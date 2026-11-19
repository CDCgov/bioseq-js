import fc from "fast-check";
import { describe, expect, it } from "vitest";
import { aa24, align, alignWith, BLOSUM62, cigar2str, nac } from "../src/index";
import { legacy, rescore, scoring, seqPair, subFn } from "./helpers";
import { referenceScore } from "./reference";

describe("README quick start", () => {
  it("aligns the documented example", () => {
    const r = align("ATAGCTAGCTAGCATAAGC", "AGCTAcCGCAT")!;
    expect(r.score).toBe(6);
    expect(r.position).toBe(6);
    expect(cigar2str(r.CIGAR)).toBe("5M2I4M");
  });
  it("legacy positional signature equals the options form", () => {
    const a = align("ATAGCTAGCTAGCATAAGC", "AGCTAcCGCAT", true, [2, -3], [4, 2]);
    const b = align("ATAGCTAGCTAGCATAAGC", "AGCTAcCGCAT", { mode: "local", matrix: [2, -3], gap: [4, 2] });
    expect(a).toEqual(b);
    const g1 = align("ACGTTGCA", "ACGTGCA", false, [1, -1], [1, 1], 3, null);
    const g2 = align("ACGTTGCA", "ACGTGCA", { mode: "global", band: 3 });
    expect(g1).toEqual(g2);
    expect(align("ACGT", "ACGT")).toEqual(align("ACGT", "ACGT", { mode: "local" }));
  });
  it("alignWith is the same function as the options form", () => {
    expect(alignWith("ACGT", "CG")).toEqual(align("ACGT", "CG"));
  });
});

describe("local", () => {
  it("returns null when everything mismatches", () => {
    expect(align("AAAA", "CCCC")).toBeNull();
    expect(align("AAAA", "CCCC", { scoreOnly: true })).toBeNull();
  });
  it("soft-clips unaligned query ends", () => {
    const r = align("GGGGACGTACGTGGGG", "TTTACGTACGTTTT")!;
    expect(cigar2str(r.CIGAR)).toBe("3S8M3S");
    expect(r.queryStart).toBe(3);
    expect(r.queryEnd).toBe(11);
    expect(r.refStart).toBe(4);
    expect(r.refEnd).toBe(12);
  });
  it("empty inputs give null", () => {
    expect(align("", "")).toBeNull();
    expect(align("", "ACGT")).toBeNull();
    expect(align("ACGT", "")).toBeNull();
  });
});

describe("global", () => {
  it("empty vs empty", () => {
    const r = align("", "", { mode: "global" })!;
    expect(r.score).toBe(0);
    expect(r.CIGAR).toEqual([]);
  });
  it("empty query is one deletion", () => {
    const r = align("ACG", "", { mode: "global" })!;
    expect(r.score).toBe(-4);
    expect(cigar2str(r.CIGAR)).toBe("3D");
  });
  it("aligns end to end", () => {
    const r = align("ACGT", "ACGT", { mode: "global" })!;
    expect(r).toMatchObject({ score: 4, refStart: 0, refEnd: 4, queryStart: 0, queryEnd: 4 });
    expect(cigar2str(r.CIGAR)).toBe("4M");
  });
  it("regression: insertion directly followed by deletion pays two full gap costs (E[0] fix)", () => {
    // mismatch costs 5, so the best alignment of C vs A is an insertion plus a deletion: -(1) + -(1).
    // bioseq 0.1.5 reported -1 because its initial E[0] omitted one extension penalty.
    const opts = { mode: "global", matrix: [1, -5], gap: [0, 1] } as const;
    const r = align("C", "A", opts)!;
    expect(r.score).toBe(-2);
    expect(referenceScore("C", "A", subFn(1, -5), 0, 1, "global")).toBe(-2);
    expect(rescore("C", "A", r, subFn(1, -5), 0, 1)).toBe(-2);
    expect(legacy.align("C", "A", false, [1, -5], [0, 1]).score).toBe(-1); // the legacy bug
  });
  it("the E[0] fix is found by search too (new always equals oracle)", () => {
    fc.assert(
      fc.property(
        seqPair.filter(([r, q]) => r.length > 0 && q.length > 0),
        scoring,
        ([r, q], s) => {
          const res = align(r, q, { mode: "global", matrix: [s.match, s.mismatch], gap: [s.open, s.extend] })!;
          expect(res.score).toBe(referenceScore(r, q, subFn(s.match, s.mismatch), s.open, s.extend, "global"));
        },
      ),
      { numRuns: 200 },
    );
  });
});

describe("semiglobal", () => {
  it("finds an exact substring at the right offset", () => {
    const r = align("TTTTGATTACATTTT", "GATTACA", { mode: "semiglobal" })!;
    expect(r).toMatchObject({ score: 7, refStart: 4, refEnd: 11, queryStart: 0, queryEnd: 7, position: 4 });
    expect(cigar2str(r.CIGAR)).toBe("7M");
  });
  it("empty query scores 0", () => {
    const r = align("ACGT", "", { mode: "semiglobal" })!;
    expect(r.score).toBe(0);
    expect(r.CIGAR).toEqual([]);
  });
  it("whole query is aligned even when it is poor", () => {
    const r = align("AAAA", "CCCC", { mode: "semiglobal" })!;
    expect(r.score).toBe(-4);
    expect(r.queryEnd - r.queryStart).toBe(4);
  });
  it("ignores band", () => {
    expect(align("TTTTTTTTGATTACA", "GATTACA", { mode: "semiglobal", band: 1 })).toEqual(
      align("TTTTTTTTGATTACA", "GATTACA", { mode: "semiglobal" }),
    );
  });
});

describe("band initialisation regressions (j >= w off-by-one inherited from ksw)", () => {
  it("global with empty reference and non-empty query scores one insertion", () => {
    const r = align("", "ACG", { mode: "global" })!;
    expect(r.score).toBe(-4);
    expect(cigar2str(r.CIGAR)).toBe("3I");
  });
  it("semiglobal with empty reference and non-empty query scores one insertion", () => {
    const r = align("", "ACG", { mode: "semiglobal" })!;
    expect(r.score).toBe(-4);
    expect(cigar2str(r.CIGAR)).toBe("3I");
  });
  it("a path may not leave the band through the first row (found by fast-check)", () => {
    // band 0 widens to 1; "2I3D" passes through cell (-1, 1), which is outside the band
    const r = align("AAA", "TC", { mode: "global", matrix: [1, -3], gap: [1, 1], band: 0 })!;
    expect(r.score).toBe(-8);
    expect(cigar2str(r.CIGAR)).not.toBe("2I3D"); // the only -7 path, which is out of band
  });
  it("a path running along the very edge of the band is allowed (band init off by one)", () => {
    // best path: 1 leading insertion, 9 matches on diagonal offset +1, 1 trailing deletion; all within |i-j| <= 1
    const ref = "CGTACGTTGA";
    const query = "ACGTACGTTG";
    expect(align(ref, query, { mode: "global", band: 1 })!.score).toBe(
      align(ref, query, { mode: "global", band: 2 })!.score,
    );
  });
});

describe("banding", () => {
  it("band restricts the search (legacy ignored it)", () => {
    const core = "GATTACAGATTACAGGT";
    const ref = `CCCCCCCC${core}`;
    const query = `${core}CCCCCCCC`;
    const unbanded = align(ref, query, { mode: "global" })!;
    const narrow = align(ref, query, { mode: "global", band: 2 })!;
    const wide = align(ref, query, { mode: "global", band: 20 })!;
    expect(narrow.score).toBeLessThan(unbanded.score);
    expect(wide).toEqual(unbanded);
    expect(narrow.CIGAR).not.toEqual(unbanded.CIGAR);
    // legacy 0.1.5 silently ignored the band
    const old = legacy.align(ref, query, false, [1, -1], [1, 1], 2);
    expect(old.score).toBe(unbanded.score);
  });
  it("is widened to the length difference", () => {
    const r = align("ACGTACGTACGT", "ACGT", { mode: "global", band: 0 })!;
    expect(r.score).toBe(align("ACGTACGTACGT", "ACGT", { mode: "global" })!.score);
  });
  it("negative or null band means unbanded", () => {
    expect(align("ACGTTGCA", "ACGTGCA", { band: -1 })).toEqual(align("ACGTTGCA", "ACGTGCA"));
  });
});

describe("scoreOnly", () => {
  it("returns only score and ends", () => {
    const r = align("ATAGCTAGCTAGCATAAGC", "AGCTAcCGCAT", { scoreOnly: true })!;
    expect(r).toEqual({ score: 6, refEnd: 15, queryEnd: 11 });
  });
});

describe("protein alignment", () => {
  const ref = "MKTAYIAKQRQISFVKSHFSRQ";
  const query = "MKTAYIAKQRQISFVKSHFSRQ".replace("ISF", "IF").replace("TAY", "TWY");
  const sub = (a: string, b: string) => BLOSUM62[aa24[a.charCodeAt(0)]][aa24[b.charCodeAt(0)]];
  for (const mode of ["local", "global", "semiglobal"] as const) {
    it(`${mode} matches oracle with aa24 + BLOSUM62`, () => {
      const r = align(ref, query, { mode, matrix: BLOSUM62, table: aa24, gap: [11, 1] })!;
      expect(r.score).toBe(referenceScore(ref, query, sub, 11, 1, mode));
      expect(rescore(ref, query, r, sub, 11, 1)).toBe(r.score);
    });
  }
  it("default matrix works with aa24 (24 symbols)", () => {
    expect(align("ARND", "ARND", { table: aa24 })!.score).toBe(4);
  });
  it("lower-case protein works", () => {
    expect(align("arndceq", "ARNDCEQ", { table: aa24, matrix: BLOSUM62 })).toEqual(
      align("ARNDCEQ", "ARNDCEQ", { table: aa24, matrix: BLOSUM62 }),
    );
  });
});

describe("matrix / table mismatches", () => {
  it("throws RangeError when a symbol code exceeds the matrix", () => {
    const m4 = [
      [1, -1, -1, -1],
      [-1, 1, -1, -1],
      [-1, -1, 1, -1],
      [-1, -1, -1, 1],
    ];
    expect(() => align("ACGN", "ACGT", { matrix: m4 })).toThrow(RangeError);
    expect(() => align("ACGT", "ACGN", { matrix: m4 })).toThrow(RangeError);
    expect(align("ACGT", "ACGT", { matrix: m4 })!.score).toBe(4);
  });
  it("table codes only need to fit inside the matrix", () => {
    expect(() => align("ACGT", "ACGT", { table: nac, matrix: BLOSUM62 })).not.toThrow();
    // nac maps unknown characters to 16, still < 24
    expect(() => align("ACGT-", "ACGT", { table: nac, matrix: BLOSUM62 })).not.toThrow();
    // '*' encodes to 23 with aa24, outside a 1x1 matrix
    expect(() => align("ACGT*", "ACGT", { table: aa24, matrix: [[1]] })).toThrow(RangeError);
  });
});

describe("input handling", () => {
  it("unknown mode throws TypeError", () => {
    // biome-ignore lint/suspicious/noExplicitAny: deliberately invalid
    expect(() => align("ACGT", "ACGT", { mode: "bogus" as any })).toThrow(TypeError);
  });
  it("is case-insensitive", () => {
    expect(align("acgtacgt", "ACGTACGT")).toEqual(align("ACGTACGT", "ACGTACGT"));
    expect(align("ACGTACGT", "acgtacgt")).toEqual(align("ACGTACGT", "ACGTACGT"));
  });
  it("accepts pre-encoded sequences", () => {
    expect(align([0, 1, 2, 3], [1, 2])).toEqual(align("ACGT", "CG"));
  });
  it("unknown characters (N) score zero", () => {
    expect(align("ACNGT", "ACNGT")!.score).toBe(4);
  });
  it("gap given as a single number is a linear penalty", () => {
    const a = align("ACGTTGCA", "ACGTGCA", { mode: "global", gap: 2 })!;
    const b = align("ACGTTGCA", "ACGTGCA", { mode: "global", gap: [0, 2] })!;
    expect(a).toEqual(b);
  });
  it("gap signs are ignored", () => {
    expect(align("ACGTTGCA", "ACGTGCA", { gap: [-3, -1] })).toEqual(align("ACGTTGCA", "ACGTGCA", { gap: [3, 1] }));
  });
});
