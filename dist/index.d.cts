/** A 256-entry lookup table mapping character codes to small integer symbol codes. */
type AlphabetTable = ArrayLike<number>;
/**
 * Allocate an integer typed array of the requested bit size, optionally filled.
 * Kept for backwards compatibility with bioseq <= 0.1.
 */
declare function makeIntArray(length: number, bitSize: number, fill?: number): Int8Array | Int16Array | Int32Array;
/**
 * Build an encoding table for an alphabet. The n-th character of `str` (in
 * either case) is encoded as n; every other character is encoded as `defaultVal`.
 */
declare function makeAlphabetMap(str: string, defaultVal: number): Int8Array;
/** Nucleotides: A=0 C=1 G=2 T=3, anything else 4 (scores 0 with `[match, mismatch]`). */
declare const nt5: Int8Array<ArrayBufferLike>;
/** IUPAC nucleic-acid codes, anything else 16. */
declare const nac: Int8Array<ArrayBufferLike>;
/** Amino acids (legacy ordering), anything else 32. */
declare const aac: Int8Array<ArrayBufferLike>;
/** Amino acids in NCBI matrix order (use with `BLOSUM62`); unknown characters map to X. */
declare const aa24: Int8Array<ArrayBufferLike>;
/** Encode a sequence string to symbol codes using `table`. */
declare function encode(seq: string, table?: AlphabetTable): Uint8Array;
/** Number of distinct symbol codes a table can produce (max code + 1). */
declare function alphabetSize(table: AlphabetTable): number;
/** Reverse-complement a DNA/RNA sequence (IUPAC aware, case preserving). Unknown characters are kept as-is. */
declare function reverseComplement(seq: string): string;

/**
 * - `local`: Smith-Waterman. Best-scoring substrings of both sequences; unaligned query ends are soft-clipped.
 * - `global`: Needleman-Wunsch. Both sequences aligned end to end.
 * - `semiglobal`: the whole query is aligned to its best-matching region of the reference
 *   (free end gaps on the reference only). Typical for mapping a read, gene or primer onto a genome.
 */
type AlignMode = "local" | "global" | "semiglobal";
/** A square substitution matrix indexed by symbol code. */
type ScoreMatrix = readonly (readonly number[])[];
interface AlignOptions {
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
interface AlignmentScore {
    score: number;
    refEnd: number;
    queryEnd: number;
}
interface Alignment extends AlignmentScore {
    /** Start of the alignment on the reference (alias of `refStart`, kept for compatibility). */
    position: number;
    /** Packed BAM-style CIGAR; see `cigar2str`. */
    CIGAR: number[];
    refStart: number;
    queryStart: number;
}
type Seq$1 = string | ArrayLike<number>;
/** Build a square match/mismatch matrix; the last symbol (ambiguous/unknown) scores 0 against everything. */
declare function genScoreMatrix(n: number, match: number, mismatch: number): number[][];
/**
 * Pairwise alignment with affine gap penalties.
 *
 * Returns `null` for a local alignment when no positive-scoring alignment exists.
 */
declare function alignWith(reference: Seq$1, query: Seq$1, options: AlignOptions & {
    scoreOnly: true;
}): AlignmentScore | null;
declare function alignWith(reference: Seq$1, query: Seq$1, options?: AlignOptions): Alignment | null;

/**
 * CIGAR operations, BAM-style: each element is `(length << 4) | op`.
 * Only M, I, D and S are produced by `align`.
 */
declare const CIGAR_OPS = "MIDNSHP=XB";
/** Minimal shape needed by the helpers below (an `Alignment` satisfies it). */
interface CigarAlignment {
    position: number;
    CIGAR: readonly number[];
}
/** Convert a packed CIGAR array to its string form, e.g. `"3S5M2I4M"`. */
declare function cigar2str(cigar: readonly number[]): string;
/** Parse a CIGAR string (e.g. `"5M2I4M"`) into the packed numeric form. */
declare function str2cigar(str: string): number[];
/**
 * Render an alignment as two gapped strings `[reference, query]`.
 *
 * With `fixed` set, insertions relative to the reference are dropped so the
 * output stays in reference coordinates (useful when aligning many queries
 * to one reference). Soft-clipped query bases are never shown.
 */
declare function cigar2gaps(reference: string, query: string, start: number, cigar: readonly number[], fixed?: boolean): [string, string];
interface AlignmentStats {
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
declare function alignmentStats(reference: string, query: string, aln: CigarAlignment): AlignmentStats;

declare const BLOSUM62: readonly (readonly number[])[];

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
declare function align(reference: Seq, query: Seq, options: AlignOptions & {
    scoreOnly: true;
}): AlignmentScore | null;
declare function align(reference: Seq, query: Seq, options?: AlignOptions): Alignment | null;
declare function align(reference: Seq, query: Seq, is_local?: boolean, matrix?: readonly [number, number] | ScoreMatrix, gapsc?: readonly [number, number] | number, w?: number | null, table?: AlphabetTable | null): Alignment | null;
/** @deprecated Use `genScoreMatrix`. */
declare const gen_score_matrix: typeof genScoreMatrix;
/** @deprecated Use `encode`, which returns a Uint8Array. */
declare function bsg_enc_seq(seq: string, table?: AlphabetTable): number[];
/**
 * @deprecated Internal preprocessing step; `align` builds its own profile.
 * Returns `profile[symbol][j]` = score of `symbol` against query position j.
 */
declare function gen_query_profile(seq: Seq, matrix: readonly [number, number] | ScoreMatrix, table?: AlphabetTable | number | null): number[][] | null;

export { BLOSUM62, CIGAR_OPS, aa24, aac, align, alignWith, alignmentStats, alphabetSize, bsg_enc_seq, cigar2gaps, cigar2str, encode, genScoreMatrix, gen_query_profile, gen_score_matrix, makeAlphabetMap, makeIntArray, nac, nt5, reverseComplement, str2cigar };
export type { AlignMode, AlignOptions, Alignment, AlignmentScore, AlignmentStats, AlphabetTable, ScoreMatrix };
