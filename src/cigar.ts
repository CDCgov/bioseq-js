/**
 * CIGAR operations, BAM-style: each element is `(length << 4) | op`.
 * Only M, I, D and S are produced by `align`.
 */
export const CIGAR_OPS = "MIDNSHP=XB";
export const CIGAR_MATCH = 0;
export const CIGAR_INS = 1;
export const CIGAR_DEL = 2;
export const CIGAR_SOFT_CLIP = 4;

/** Minimal shape needed by the helpers below (an `Alignment` satisfies it). */
export interface CigarAlignment {
  position: number;
  CIGAR: readonly number[];
}

/** Convert a packed CIGAR array to its string form, e.g. `"3S5M2I4M"`. */
export function cigar2str(cigar: readonly number[]): string {
  let s = "";
  for (let k = 0; k < cigar.length; k++) s += (cigar[k] >>> 4).toString() + CIGAR_OPS.charAt(cigar[k] & 0xf);
  return s;
}

/** Parse a CIGAR string (e.g. `"5M2I4M"`) into the packed numeric form. */
export function str2cigar(str: string): number[] {
  const out: number[] = [];
  const re = /(\d+)([MIDNSHP=XB])/gy;
  let m: RegExpExecArray | null = re.exec(str);
  let consumed = 0;
  while (m !== null) {
    out.push((Number(m[1]) << 4) | CIGAR_OPS.indexOf(m[2]));
    consumed = re.lastIndex;
    m = re.exec(str);
  }
  if (consumed !== str.length) throw new SyntaxError(`Invalid CIGAR string: ${str}`);
  return out;
}

/**
 * Render an alignment as two gapped strings `[reference, query]`.
 *
 * With `fixed` set, insertions relative to the reference are dropped so the
 * output stays in reference coordinates (useful when aligning many queries
 * to one reference). Soft-clipped query bases are never shown.
 */
export function cigar2gaps(
  reference: string,
  query: string,
  start: number,
  cigar: readonly number[],
  fixed?: boolean,
): [string, string] {
  let oq = "";
  let ot = "";
  let lq = 0;
  let lt = start;
  for (let k = 0; k < cigar.length; k++) {
    const op = cigar[k] & 0xf;
    const len = cigar[k] >>> 4;
    if (op === CIGAR_MATCH || op === 7 || op === 8) {
      oq += query.substr(lq, len);
      ot += reference.substr(lt, len);
      lq += len;
      lt += len;
    } else if (op === CIGAR_INS) {
      if (!fixed) {
        oq += query.substr(lq, len);
        ot += "-".repeat(len);
      }
      lq += len;
    } else if (op === CIGAR_DEL) {
      oq += "-".repeat(len);
      ot += reference.substr(lt, len);
      lt += len;
    } else if (op === CIGAR_SOFT_CLIP) {
      lq += len;
    }
  }
  return [ot, oq];
}

export interface AlignmentStats {
  /** Aligned columns where the residues are identical (case-insensitive). */
  matches: number;
  /** Aligned columns where the residues differ. */
  mismatches: number;
  /** Query residues inserted relative to the reference. */
  insertions: number;
  /** Reference residues deleted from the query. */
  deletions: number;
  /** Number of distinct gap runs (I or D operations). */
  gapOpens: number;
  /** Alignment columns, excluding soft clips. */
  length: number;
  /** matches / length (0 for an empty alignment). */
  identity: number;
  /** Extended CIGAR using `=`/`X` instead of `M`, e.g. `"4=1X3=2I"`. */
  extendedCigar: string;
}

/** Compute identity, mismatch and gap statistics plus an `=`/`X` CIGAR for an alignment. */
export function alignmentStats(reference: string, query: string, aln: CigarAlignment): AlignmentStats {
  let matches = 0;
  let mismatches = 0;
  let insertions = 0;
  let deletions = 0;
  let gapOpens = 0;
  let ext = "";
  let lq = 0;
  let lt = aln.position;
  for (const c of aln.CIGAR) {
    const op = c & 0xf;
    const len = c >>> 4;
    if (op === CIGAR_MATCH || op === 7 || op === 8) {
      let runOp = "";
      let runLen = 0;
      for (let n = 0; n < len; n++) {
        const same = reference[lt + n].toUpperCase() === query[lq + n].toUpperCase();
        if (same) matches++;
        else mismatches++;
        const o = same ? "=" : "X";
        if (o === runOp) runLen++;
        else {
          if (runLen) ext += runLen + runOp;
          runOp = o;
          runLen = 1;
        }
      }
      if (runLen) ext += runLen + runOp;
      lq += len;
      lt += len;
    } else {
      ext += len + CIGAR_OPS.charAt(op);
      if (op === CIGAR_INS) {
        insertions += len;
        gapOpens++;
        lq += len;
      } else if (op === CIGAR_DEL) {
        deletions += len;
        gapOpens++;
        lt += len;
      } else if (op === CIGAR_SOFT_CLIP) {
        lq += len;
      }
    }
  }
  const length = matches + mismatches + insertions + deletions;
  return {
    matches,
    mismatches,
    insertions,
    deletions,
    gapOpens,
    length,
    identity: length ? matches / length : 0,
    extendedCigar: ext,
  };
}
