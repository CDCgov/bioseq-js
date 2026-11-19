/** A 256-entry lookup table mapping character codes to small integer symbol codes. */
export type AlphabetTable = ArrayLike<number>;

/**
 * Allocate an integer typed array of the requested bit size, optionally filled.
 * Kept for backwards compatibility with bioseq <= 0.1.
 */
export function makeIntArray(length: number, bitSize: number, fill?: number): Int8Array | Int16Array | Int32Array {
  const arr = bitSize <= 8 ? new Int8Array(length) : bitSize <= 16 ? new Int16Array(length) : new Int32Array(length);
  if (fill !== undefined && fill !== 0) arr.fill(fill);
  return arr;
}

/**
 * Build an encoding table for an alphabet. The n-th character of `str` (in
 * either case) is encoded as n; every other character is encoded as `defaultVal`.
 */
export function makeAlphabetMap(str: string, defaultVal: number): Int8Array {
  const aMap = makeIntArray(256, 8, defaultVal) as Int8Array;
  const lstr = str.toLowerCase();
  for (let n = 0; n < str.length; n++) {
    aMap[str.charCodeAt(n)] = n;
    aMap[lstr.charCodeAt(n)] = n;
  }
  return aMap;
}

/** Nucleotides: A=0 C=1 G=2 T=3, anything else 4 (scores 0 with `[match, mismatch]`). */
export const nt5 = makeAlphabetMap("ACGT", 4);
/** IUPAC nucleic-acid codes, anything else 16. */
export const nac = makeAlphabetMap("ACGTURYMKWSBDHVN", 16);
/** Amino acids (legacy ordering), anything else 32. */
export const aac = makeAlphabetMap("ARNDCQEGHILKMFPSTWYVBZX", 32);
/** Amino acids in NCBI matrix order (use with `BLOSUM62`); unknown characters map to X. */
export const aa24 = makeAlphabetMap("ARNDCQEGHILKMFPSTWYVBZX*", 22);

/** Encode a sequence string to symbol codes using `table`. */
export function encode(seq: string, table: AlphabetTable = nt5): Uint8Array {
  const s = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i);
    // Non-Latin-1 characters are treated like any other unknown character (table[0] is NUL's code).
    s[i] = table[c < 256 ? c : 0];
  }
  return s;
}

/** Number of distinct symbol codes a table can produce (max code + 1). */
export function alphabetSize(table: AlphabetTable): number {
  let max = 0;
  for (let i = 0; i < table.length; i++) if (table[i] > max) max = table[i];
  return max + 1;
}

const COMPLEMENT: Record<string, string> = {
  A: "T",
  C: "G",
  G: "C",
  T: "A",
  U: "A",
  R: "Y",
  Y: "R",
  S: "S",
  W: "W",
  K: "M",
  M: "K",
  B: "V",
  V: "B",
  D: "H",
  H: "D",
  N: "N",
};

/** Reverse-complement a DNA/RNA sequence (IUPAC aware, case preserving). Unknown characters are kept as-is. */
export function reverseComplement(seq: string): string {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) {
    const c = seq[i];
    const up = c.toUpperCase();
    const comp = COMPLEMENT[up];
    out += comp === undefined ? c : c === up ? comp : comp.toLowerCase();
  }
  return out;
}
