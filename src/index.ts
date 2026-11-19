import {
  type Alignment,
  type AlignmentScore,
  type AlignOptions,
  alignWith,
  genScoreMatrix,
  type ScoreMatrix,
} from "./align";
import { type AlphabetTable, encode, nt5 } from "./alphabet";

export {
  type AlignMode,
  type Alignment,
  type AlignmentScore,
  type AlignOptions,
  alignWith,
  genScoreMatrix,
  type ScoreMatrix,
} from "./align";
export {
  type AlphabetTable,
  aa24,
  aac,
  alphabetSize,
  encode,
  makeAlphabetMap,
  makeIntArray,
  nac,
  nt5,
  reverseComplement,
} from "./alphabet";
export {
  type AlignmentStats,
  alignmentStats,
  CIGAR_OPS,
  cigar2gaps,
  cigar2str,
  str2cigar,
} from "./cigar";
export { BLOSUM62 } from "./matrices";

type Seq = string | ArrayLike<number>;

/**
 * Align `query` against `reference`.
 *
 * Preferred form: `align(reference, query, { mode, matrix, gap, band, table, scoreOnly })`.
 *
 * The positional form `align(reference, query, is_local, matrix, gapsc, w, table)` from
 * bioseq 0.1 is still supported.
 *
 * Returns `null` for a local alignment when nothing scores above zero.
 */
export function align(reference: Seq, query: Seq, options: AlignOptions & { scoreOnly: true }): AlignmentScore | null;
export function align(reference: Seq, query: Seq, options?: AlignOptions): Alignment | null;
export function align(
  reference: Seq,
  query: Seq,
  is_local?: boolean,
  matrix?: readonly [number, number] | ScoreMatrix,
  gapsc?: readonly [number, number] | number,
  w?: number | null,
  table?: AlphabetTable | null,
): Alignment | null;
export function align(
  reference: Seq,
  query: Seq,
  optsOrLocal?: AlignOptions | boolean,
  matrix?: readonly [number, number] | ScoreMatrix,
  gapsc?: readonly [number, number] | number,
  w?: number | null,
  table?: AlphabetTable | null,
): Alignment | AlignmentScore | null {
  if (optsOrLocal !== null && typeof optsOrLocal === "object") return alignWith(reference, query, optsOrLocal);
  return alignWith(reference, query, {
    mode: optsOrLocal === undefined || optsOrLocal ? "local" : "global",
    matrix,
    gap: gapsc,
    band: w ?? undefined,
    table: table ?? undefined,
  });
}

// ---------------------------------------------------------------------------
// Legacy (0.1.x) helpers, kept for compatibility.
// ---------------------------------------------------------------------------

/** @deprecated Use `genScoreMatrix`. */
export const gen_score_matrix = genScoreMatrix;

/** @deprecated Use `encode`, which returns a Uint8Array. */
export function bsg_enc_seq(seq: string, table: AlphabetTable = nt5): number[] {
  return Array.from(encode(seq, table));
}

/**
 * @deprecated Internal preprocessing step; `align` builds its own profile.
 * Returns `profile[symbol][j]` = score of `symbol` against query position j.
 */
export function gen_query_profile(
  seq: Seq,
  matrix: readonly [number, number] | ScoreMatrix,
  table?: AlphabetTable | number | null,
): number[][] | null {
  const s = typeof seq === "string" ? bsg_enc_seq(seq, (typeof table === "object" && table) || nt5) : seq;
  let m: ScoreMatrix;
  if (matrix.length >= 2 && typeof matrix[0] === "number" && typeof matrix[1] === "number") {
    if (table == null) return null;
    const n = typeof table === "number" ? table : table[table.length - 1] + 1;
    m = genScoreMatrix(n, matrix[0], matrix[1]);
  } else {
    m = matrix as ScoreMatrix;
  }
  return m.map((row) => Array.from(s, (c) => row[c]));
}
