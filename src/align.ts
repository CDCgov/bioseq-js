import { type AlphabetTable, alphabetSize, encode, nt5 } from "./alphabet";
import { CIGAR_DEL, CIGAR_INS, CIGAR_MATCH, CIGAR_SOFT_CLIP } from "./cigar";

/**
 * - `local`: Smith-Waterman. Best-scoring substrings of both sequences; unaligned query ends are soft-clipped.
 * - `global`: Needleman-Wunsch. Both sequences aligned end to end.
 * - `semiglobal`: the whole query is aligned to its best-matching region of the reference
 *   (free end gaps on the reference only). Typical for mapping a read, gene or primer onto a genome.
 */
export type AlignMode = "local" | "global" | "semiglobal";

/** A square substitution matrix indexed by symbol code. */
export type ScoreMatrix = readonly (readonly number[])[];

export interface AlignOptions {
  /** Default `"local"`. */
  mode?: AlignMode;
  /** `[match, mismatch]` scores or a full substitution matrix. Default `[1, -1]`. */
  matrix?: readonly [number, number] | ScoreMatrix;
  /**
   * Affine gap penalties `[open, extend]`; a gap of length k costs open + k * extend.
   * A single number means a linear gap penalty. Signs are ignored. Default `[1, 1]`.
   */
  gap?: readonly [number, number] | number;
  /**
   * Band width: only cells within this distance of the main diagonal are computed.
   * Greatly reduces time and memory for similar-length sequences. Automatically widened to at
   * least |reference.length - query.length|. Ignored for `semiglobal`. Default: unbanded.
   */
  band?: number;
  /** Character encoding table. Default `nt5`. Use `aa24` with `BLOSUM62`. */
  table?: AlphabetTable;
  /** Skip the traceback: returns only score and end coordinates, in O(query) memory. */
  scoreOnly?: boolean;
}

/** Coordinates are 0-based, end-exclusive. */
export interface AlignmentScore {
  score: number;
  refEnd: number;
  queryEnd: number;
}

export interface Alignment extends AlignmentScore {
  /** Start of the alignment on the reference (alias of `refStart`, kept for compatibility). */
  position: number;
  /** Packed BAM-style CIGAR; see `cigar2str`. */
  CIGAR: number[];
  refStart: number;
  queryStart: number;
}

type Seq = string | ArrayLike<number>;

// Large enough to never win a max(), small enough that subtracting gap penalties can't overflow int32.
const NEG_INF = -(1 << 29);
const MAX_TRACEBACK_CELLS = 2 ** 31;

function absPair(gap: AlignOptions["gap"]): [number, number] {
  if (gap === undefined) return [1, 1];
  if (typeof gap === "number") return [0, Math.abs(gap)];
  return [Math.abs(gap[0]), Math.abs(gap[1])];
}

function isPair(m: AlignOptions["matrix"]): m is readonly [number, number] {
  return m !== undefined && m.length === 2 && typeof m[0] === "number" && typeof m[1] === "number";
}

/** Build a square match/mismatch matrix; the last symbol (ambiguous/unknown) scores 0 against everything. */
export function genScoreMatrix(n: number, match: number, mismatch: number): number[][] {
  const b = mismatch > 0 ? -mismatch : mismatch;
  const m: number[][] = [];
  for (let i = 0; i < n; i++) {
    m[i] = [];
    for (let j = 0; j < n; j++) m[i][j] = i === n - 1 || j === n - 1 ? 0 : i === j ? match : b;
  }
  return m;
}

function checkCodes(s: ArrayLike<number>, n: number, what: string): void {
  for (let i = 0; i < s.length; i++) {
    if (!(s[i] >= 0 && s[i] < n)) {
      throw new RangeError(
        `${what} position ${i} encodes to symbol ${s[i]}, outside the ${n}x${n} score matrix. ` +
          "Check that the encoding table matches the matrix (e.g. use aa24 with BLOSUM62).",
      );
    }
  }
}

/**
 * Pairwise alignment with affine gap penalties.
 *
 * Returns `null` for a local alignment when no positive-scoring alignment exists.
 */
export function alignWith(
  reference: Seq,
  query: Seq,
  options: AlignOptions & { scoreOnly: true },
): AlignmentScore | null;
export function alignWith(reference: Seq, query: Seq, options?: AlignOptions): Alignment | null;
export function alignWith(reference: Seq, query: Seq, options: AlignOptions = {}): Alignment | AlignmentScore | null {
  const mode = options.mode ?? "local";
  if (mode !== "local" && mode !== "global" && mode !== "semiglobal") throw new TypeError(`Unknown mode: ${mode}`);
  const isLocal = mode === "local";
  const isSemi = mode === "semiglobal";
  const table = options.table ?? nt5;
  const scoreOnly = options.scoreOnly === true;

  const t = typeof reference === "string" ? encode(reference, table) : reference;
  const q = typeof query === "string" ? encode(query, table) : query;
  const tlen = t.length;
  const qlen = q.length;

  // Substitution matrix -> query profile qp[a * qlen + j] = S(a, q[j])
  const matrix = isPair(options.matrix)
    ? genScoreMatrix(alphabetSize(table), options.matrix[0], options.matrix[1])
    : (options.matrix ?? genScoreMatrix(alphabetSize(table), 1, -1));
  const nsym = matrix.length;
  checkCodes(t, nsym, "Reference");
  checkCodes(q, nsym, "Query");
  const qp = new Int32Array(nsym * qlen);
  for (let a = 0; a < nsym; a++) {
    const row = matrix[a];
    for (let j = 0; j < qlen; j++) qp[a * qlen + j] = row[q[j]];
  }

  // Band width
  const maxLen = qlen > tlen ? qlen : tlen;
  let w = options.band == null || options.band < 0 || isSemi ? maxLen : Math.floor(options.band);
  const lenDiff = tlen > qlen ? tlen - qlen : qlen - tlen;
  if (w < lenDiff) w = lenDiff;

  const [gapo, gape] = absPair(options.gap);
  const gapoe = gapo + gape;

  // H[j] / E[j] hold the previous row's H and the current row's E; see README "Theory".
  const H = new Int32Array(qlen + 1);
  const E = new Int32Array(qlen + 1);
  if (!isLocal) {
    H[0] = 0;
    E[0] = w >= 1 ? -(gapoe + gapoe) : NEG_INF; // opens from H(-1, 0)
    // Row -1 cell (-1, j) is inside the band iff j + 1 <= w. H[j] holds H(-1, j-1);
    // E[j] = E(0, j) opens from H(-1, j). Out-of-band predecessors are -inf.
    for (let j = 1; j <= qlen; j++) {
      H[j] = j > w ? NEG_INF : -(gapoe + gape * (j - 1));
      E[j] = j >= w ? NEG_INF : -(gapoe + gapoe + gape * j);
    }
  }

  // Traceback: one byte per in-band cell.
  const bw = qlen < 2 * w + 1 ? qlen : 2 * w + 1;
  if (!scoreOnly && tlen * bw > MAX_TRACEBACK_CELLS) {
    throw new RangeError(
      `Traceback matrix would need ${tlen * bw} cells; use the band option, scoreOnly, or shorter sequences.`,
    );
  }
  const z = scoreOnly ? new Uint8Array(0) : new Uint8Array(tlen * bw);

  let max = 0;
  let endI = -1;
  let endJ = -1;
  let semiBest = NEG_INF;

  for (let i = 0; i < tlen; i++) {
    let h1 = 0;
    let f = 0;
    let m = 0;
    let mj = -1;
    const qoff = t[i] * qlen;
    const beg = i > w ? i - w : 0;
    const end = i + w + 1 < qlen ? i + w + 1 : qlen;
    const zoff = i * bw - beg;
    if (isSemi) {
      h1 = 0; // leading reference residues are free
      f = -gapoe;
    } else if (!isLocal) {
      // H(i, -1) is inside the band iff i + 1 <= w; F(i, 0) opens from it.
      h1 = i >= w ? NEG_INF : -(gapoe + gape * i);
      f = i >= w ? NEG_INF : -(gapoe + gapoe + gape * i);
    }
    for (let j = beg; j < end; j++) {
      // h = H(i-1,j-1), e = E(i,j), f = F(i,j), h1 = H(i,j-1)
      let e = E[j];
      let h = H[j];
      let d: number;
      H[j] = h1;
      h += qp[qoff + j];
      d = h > e ? 0 : 1;
      h = h > e ? h : e;
      d = h > f ? d : 2;
      h = h > f ? h : f; // h = H(i,j)
      if (isLocal && !(h > 0)) d = 64;
      h1 = h;
      mj = m > h ? mj : j;
      m = m > h ? m : h;
      h -= gapoe;
      if (isLocal && h < 0) h = 0;
      e -= gape;
      d |= e > h ? 4 : 0;
      e = e > h ? e : h; // E(i+1,j)
      E[j] = e;
      f -= gape;
      d |= f > h ? 32 : 0;
      f = f > h ? f : h; // F(i,j+1)
      if (!scoreOnly) z[zoff + j] = d;
    }
    H[end] = h1;
    E[end] = isLocal ? 0 : NEG_INF;
    if (isLocal) {
      if (m > max) {
        max = m;
        endI = i;
        endJ = mj;
      }
    } else if (isSemi && h1 > semiBest) {
      semiBest = h1;
      endI = i;
    }
  }

  let score: number;
  let refEnd: number;
  let queryEnd: number;
  if (isLocal) {
    if (max <= 0) return null;
    score = max;
    refEnd = endI + 1;
    queryEnd = endJ + 1;
  } else if (isSemi) {
    if (qlen === 0) {
      score = 0;
      endI = -1;
    } else if (tlen === 0) {
      score = H[qlen];
    } else {
      score = semiBest;
    }
    refEnd = endI + 1;
    queryEnd = qlen;
  } else {
    score = H[qlen];
    refEnd = tlen;
    queryEnd = qlen;
  }
  if (scoreOnly) return { score, refEnd, queryEnd };

  // Traceback
  const cigar: number[] = [];
  const push = (op: number, len: number): void => {
    const last = cigar.length - 1;
    if (last < 0 || op !== (cigar[last] & 0xf)) cigar.push((len << 4) | op);
    else cigar[last] += len << 4;
  };
  let i: number;
  let k: number;
  if (isLocal) {
    i = endI;
    k = endJ;
    if (endJ !== qlen - 1) push(CIGAR_SOFT_CLIP, qlen - 1 - endJ);
  } else if (isSemi) {
    i = endI;
    k = qlen - 1;
  } else {
    i = tlen - 1;
    k = (i + w + 1 < qlen ? i + w + 1 : qlen) - 1;
  }
  let which = 0;
  while (i >= 0 && k >= 0) {
    const tmp = z[i * bw + k - (i > w ? i - w : 0)];
    which = (tmp >> (which << 1)) & 3;
    if (which === 0 && tmp >> 6) break;
    if (which === 0) which = tmp & 3;
    if (which === 0) {
      push(CIGAR_MATCH, 1);
      i--;
      k--;
    } else if (which === 1) {
      push(CIGAR_DEL, 1);
      i--;
    } else {
      push(CIGAR_INS, 1);
      k--;
    }
  }
  if (isLocal) {
    if (k >= 0) push(CIGAR_SOFT_CLIP, k + 1);
  } else {
    if (i >= 0 && !isSemi) push(CIGAR_DEL, i + 1);
    if (k >= 0) push(CIGAR_INS, k + 1);
  }
  cigar.reverse();
  const refStart = isSemi || isLocal ? i + 1 : 0;
  const queryStart = isLocal ? k + 1 : 0;
  return { score, position: refStart, CIGAR: cigar, refStart, refEnd, queryStart, queryEnd };
}
