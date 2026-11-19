import fc from "fast-check";
import { describe, expect, it } from "vitest";
import {
  aa24,
  aac,
  alphabetSize,
  BLOSUM62,
  bsg_enc_seq,
  encode,
  gen_query_profile,
  gen_score_matrix,
  genScoreMatrix,
  makeAlphabetMap,
  makeIntArray,
  nac,
  nt5,
  reverseComplement,
} from "../src/index";
import { legacy } from "./helpers";

describe("encode / alphabets", () => {
  it("encodes nt5, both cases, unknown as 4", () => {
    expect(Array.from(encode("ACGTacgtNx"))).toEqual([0, 1, 2, 3, 0, 1, 2, 3, 4, 4]);
  });
  it("encodes with other tables", () => {
    expect(Array.from(encode("ACGU", nac))).toEqual([0, 1, 2, 4]);
    expect(Array.from(encode("*", aa24))).toEqual([23]);
    expect(Array.from(encode("J", aa24))).toEqual([22]);
    expect(Array.from(encode("J", aac))).toEqual([32]);
  });
  it("treats non-Latin-1 characters as unknown", () => {
    expect(Array.from(encode("A中C"))).toEqual([0, nt5[0], 1]);
  });
  it("makeAlphabetMap", () => {
    const m = makeAlphabetMap("XYZ", 9);
    expect(m.length).toBe(256);
    expect([m[88], m[121], m[90], m[65]]).toEqual([0, 1, 2, 9]);
    expect(m[120]).toBe(0);
  });
  it("alphabetSize", () => {
    expect(alphabetSize(nt5)).toBe(5);
    expect(alphabetSize(nac)).toBe(17);
    expect(alphabetSize(aa24)).toBe(24);
    expect(alphabetSize(aac)).toBe(33);
    expect(alphabetSize(makeAlphabetMap("AB", 0))).toBe(2);
  });
  it("makeIntArray picks the array width", () => {
    expect(makeIntArray(3, 8)).toBeInstanceOf(Int8Array);
    expect(makeIntArray(3, 16)).toBeInstanceOf(Int16Array);
    expect(makeIntArray(3, 32)).toBeInstanceOf(Int32Array);
    expect(Array.from(makeIntArray(3, 8, 7))).toEqual([7, 7, 7]);
  });
});

describe("reverseComplement", () => {
  it("handles plain DNA and RNA", () => {
    expect(reverseComplement("AACG")).toBe("CGTT");
    expect(reverseComplement("ACGU")).toBe("ACGT");
  });
  it("handles IUPAC codes", () => {
    expect(reverseComplement("RYKMBDHVSWN")).toBe("NWSBDHVKMRY");
    expect(reverseComplement("R")).toBe("Y");
    expect(reverseComplement("K")).toBe("M");
    expect(reverseComplement("B")).toBe("V");
    expect(reverseComplement("D")).toBe("H");
  });
  it("preserves case", () => {
    expect(reverseComplement("AcGt")).toBe("aCgT");
  });
  it("keeps unknown characters", () => {
    expect(reverseComplement("A-*1")).toBe("1*-T");
  });
  it("empty string", () => {
    expect(reverseComplement("")).toBe("");
  });
  it("is an involution for symmetric-IUPAC input", () => {
    fc.assert(
      fc.property(fc.string({ unit: fc.constantFrom(..."ACGTRYKMBDHVSWNacgtrykmbdhvswn-"), maxLength: 50 }), (s) => {
        expect(reverseComplement(reverseComplement(s))).toBe(s);
      }),
    );
  });
});

describe("BLOSUM62", () => {
  it("is 24x24 and symmetric", () => {
    expect(BLOSUM62.length).toBe(24);
    for (let i = 0; i < 24; i++) {
      expect(BLOSUM62[i].length).toBe(24);
      for (let j = 0; j < 24; j++) expect(BLOSUM62[i][j]).toBe(BLOSUM62[j][i]);
    }
  });
  it("spot checks via aa24", () => {
    const s = (a: string, b: string) => BLOSUM62[aa24[a.charCodeAt(0)]][aa24[b.charCodeAt(0)]];
    expect(s("A", "A")).toBe(4);
    expect(s("W", "W")).toBe(11);
    expect(s("*", "*")).toBe(1);
    expect(s("A", "*")).toBe(-4);
    expect(s("C", "C")).toBe(9);
    expect(s("A", "R")).toBe(-1);
    expect(s("W", "F")).toBe(1);
    expect(s("I", "V")).toBe(3);
  });
});

describe("genScoreMatrix and legacy helpers", () => {
  it("builds a square matrix with a zero last row/column", () => {
    expect(genScoreMatrix(3, 2, -3)).toEqual([
      [2, -3, 0],
      [-3, 2, 0],
      [0, 0, 0],
    ]);
  });
  it("a positive mismatch is negated", () => {
    expect(genScoreMatrix(3, 2, 3)).toEqual(genScoreMatrix(3, 2, -3));
  });
  it("gen_score_matrix matches legacy", () => {
    for (const [n, a, b] of [
      [5, 1, -1],
      [17, 2, -3],
      [3, 5, 4],
    ]) {
      expect(gen_score_matrix(n, a, b)).toEqual(legacy.gen_score_matrix(n, a, b));
    }
  });
  it("bsg_enc_seq matches legacy", () => {
    for (const s of ["", "ACGTNacgt", "AC-GT"]) expect(bsg_enc_seq(s)).toEqual(Array.from(legacy.bsg_enc_seq(s, nt5)));
    expect(bsg_enc_seq("ACGU", nac)).toEqual(Array.from(legacy.bsg_enc_seq("ACGU", legacy.nac)));
    expect(bsg_enc_seq("ARND", aac)).toEqual(Array.from(legacy.bsg_enc_seq("ARND", legacy.aac)));
  });
  it("gen_query_profile matches legacy for [match, mismatch] and full matrices", () => {
    expect(gen_query_profile("ACGTN", [1, -1], nt5)).toEqual(legacy.gen_query_profile("ACGTN", [1, -1], legacy.nt5));
    expect(gen_query_profile([0, 1, 2, 3, 4], [2, -3], 5)).toEqual(
      legacy.gen_query_profile([0, 1, 2, 3, 4], [2, -3], 5),
    );
    const m = genScoreMatrix(5, 1, -2);
    expect(gen_query_profile("ACGT", m, nt5)).toEqual(legacy.gen_query_profile("ACGT", m, legacy.nt5));
    expect(gen_query_profile([0, 1, 2], [1, -1], nt5)).toEqual(
      legacy.gen_query_profile([0, 1, 2], [1, -1], legacy.nt5),
    );
  });
  it("gen_query_profile returns null for a pair with no table", () => {
    expect(gen_query_profile("ACGT", [1, -1])).toBeNull();
  });
});
