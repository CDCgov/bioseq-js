import { describe, expect, it } from "vitest";
import { align, alignmentStats, cigar2gaps, cigar2str, str2cigar } from "../src/index";

const packed = (s: string) => str2cigar(s);

describe("cigar2str / str2cigar", () => {
  it("joins without commas", () => {
    expect(cigar2str([(5 << 4) | 0, (2 << 4) | 1, (4 << 4) | 0])).toBe("5M2I4M");
    expect(cigar2str([])).toBe("");
  });
  it("round-trips all ops", () => {
    const s = "3S5M2I4M1D2N1H1P3=4X1B";
    expect(cigar2str(str2cigar(s))).toBe(s);
  });
  it("packs as (len << 4) | op", () => {
    expect(str2cigar("10M2D")).toEqual([160, 34]);
  });
  it("empty string parses to empty array", () => {
    expect(str2cigar("")).toEqual([]);
  });
  it("rejects invalid input with SyntaxError", () => {
    for (const bad of ["M", "5", "5Q", "5M,2I", "5M 2I", "abc", "5M2", "-5M"]) {
      expect(() => str2cigar(bad)).toThrow(SyntaxError);
    }
  });
});

describe("cigar2gaps", () => {
  const ref = "ATAGCTAGCTAGCATAAGC";
  const query = "AGCTAcCGCAT";
  const aln = align(ref, query)!;
  it("renders insertions as gaps in the reference", () => {
    expect(cigar2gaps(ref, query, aln.position, aln.CIGAR)).toEqual(["AGCTA--GCAT", "AGCTAcCGCAT"]);
  });
  it("fixed drops insertions to stay in reference coordinates", () => {
    expect(cigar2gaps(ref, query, aln.position, aln.CIGAR, true)).toEqual(["AGCTAGCAT", "AGCTAGCAT"]);
  });
  it("shows deletions as gaps in the query", () => {
    expect(cigar2gaps("ACGGT", "ACGT", 0, packed("3M1D1M"))).toEqual(["ACGGT", "ACG-T"]);
  });
  it("omits soft-clipped bases", () => {
    expect(cigar2gaps("GGACGTGG", "TTACGTTT", 2, packed("2S4M2S"))).toEqual(["ACGT", "ACGT"]);
  });
  it("treats = and X like M", () => {
    expect(cigar2gaps("ACGT", "AGGT", 0, packed("1=1X2="))).toEqual(["ACGT", "AGGT"]);
  });
});

describe("alignmentStats", () => {
  it("hand-checked example with a mismatch and an insertion", () => {
    // ref  ACGT--ACG
    // qry  ACCTTTACG   -> =,=,X,=,I,I,=,=,=
    const st = alignmentStats("ACGTACG", "ACCTTTACG", { position: 0, CIGAR: packed("4M2I3M") });
    expect(st).toEqual({
      matches: 6,
      mismatches: 1,
      insertions: 2,
      deletions: 0,
      gapOpens: 1,
      length: 9,
      identity: 6 / 9,
      extendedCigar: "2=1X1=2I3=",
    });
  });
  it("counts deletions and separate gap runs", () => {
    const st = alignmentStats("ACGGTAC", "ACTC", { position: 0, CIGAR: packed("2M2D1M1D1M") });
    // ACGGTAC vs AC--T-C: ref positions 0-1 match, delete GG, T vs T match, delete A, C vs C
    expect(st.matches).toBe(4);
    expect(st.mismatches).toBe(0);
    expect(st.deletions).toBe(3);
    expect(st.gapOpens).toBe(2);
    expect(st.length).toBe(7);
    expect(st.extendedCigar).toBe("2=2D1=1D1=");
  });
  it("is case-insensitive; soft clips are not counted but kept in the extended CIGAR", () => {
    const st = alignmentStats("acgt", "TTACGT", { position: 0, CIGAR: packed("2S4M") });
    expect(st).toMatchObject({ matches: 4, mismatches: 0, length: 4, identity: 1, extendedCigar: "2S4=" });
  });
  it("honours the alignment start position", () => {
    const st = alignmentStats("TTTACGT", "ACGT", { position: 3, CIGAR: packed("4M") });
    expect(st.matches).toBe(4);
  });
  it("empty alignment has identity 0", () => {
    expect(alignmentStats("", "", { position: 0, CIGAR: [] }).identity).toBe(0);
  });
  it("works on real alignment output", () => {
    const ref = "ATAGCTAGCTAGCATAAGC";
    const query = "AGCTAcCGCAT";
    const aln = align(ref, query)!;
    const st = alignmentStats(ref, query, aln);
    expect(st).toMatchObject({ matches: 9, mismatches: 0, insertions: 2, gapOpens: 1, length: 11 });
    expect(st.extendedCigar).toBe("5=2I4=");
  });
});
