import fc from "fast-check";
import { describe, expect, it } from "vitest";
import { align } from "../src/index";
import { legacy, opSum, rescore, scoring, seqPair, subFn } from "./helpers";
import { type Mode, referenceScore } from "./reference";

const runs = { numRuns: 300 };
const modes: Mode[] = ["local", "global", "semiglobal"];
const prop = (
  f: (ref: string, q: string, s: { match: number; mismatch: number; open: number; extend: number }) => void,
) =>
  fc.assert(
    fc.property(seqPair, scoring, (p, s) => f(p[0], p[1], s)),
    runs,
  );
const opts = (mode: Mode, s: { match: number; mismatch: number; open: number; extend: number }) => ({
  mode,
  matrix: [s.match, s.mismatch] as [number, number],
  gap: [s.open, s.extend] as [number, number],
});

describe("score equals naive oracle", () => {
  for (const mode of modes) {
    it(mode, () =>
      prop((ref, q, s) => {
        const expected = referenceScore(ref, q, subFn(s.match, s.mismatch), s.open, s.extend, mode);
        const r = align(ref, q, opts(mode, s));
        if (mode === "local") {
          if (expected <= 0) expect(r).toBeNull();
          else expect(r?.score).toBe(expected);
        } else expect(r?.score).toBe(expected);
      }),
    );
  }
});

describe("CIGAR consistency", () => {
  for (const mode of modes) {
    it(mode, () =>
      prop((ref, q, s) => {
        const r = align(ref, q, opts(mode, s));
        if (r === null) return;
        expect(rescore(ref, q, r, subFn(s.match, s.mismatch), s.open, s.extend)).toBe(r.score);
        expect(opSum(r, "MIS")).toBe(q.length);
        expect(opSum(r, "MD")).toBe(r.refEnd - r.refStart);
        expect(r.position).toBe(r.refStart);
        if (mode === "global") {
          expect(opSum(r, "S")).toBe(0);
          expect(r.refStart).toBe(0);
          expect(r.refEnd).toBe(ref.length);
        }
        if (mode === "semiglobal") {
          expect(opSum(r, "S")).toBe(0);
          expect(r.queryStart).toBe(0);
          expect(r.queryEnd).toBe(q.length);
        }
        if (mode === "local") expect(r.queryEnd - r.queryStart).toBe(opSum(r, "MI"));
      }),
    );
  }
});

describe("scoreOnly", () => {
  for (const mode of modes) {
    it(mode, () =>
      prop((ref, q, s) => {
        const full = align(ref, q, opts(mode, s));
        const so = align(ref, q, { ...opts(mode, s), scoreOnly: true });
        if (full === null) expect(so).toBeNull();
        else expect(so).toEqual({ score: full.score, refEnd: full.refEnd, queryEnd: full.queryEnd });
      }),
    );
  }
});

describe("band", () => {
  it("band >= max length is identical to unbanded", () =>
    prop((ref, q, s) => {
      for (const mode of ["local", "global"] as Mode[]) {
        const band = Math.max(ref.length, q.length);
        expect(align(ref, q, { ...opts(mode, s), band })).toEqual(align(ref, q, opts(mode, s)));
      }
    }));
  it("banded global score equals the banded oracle", () =>
    fc.assert(
      fc.property(seqPair, scoring, fc.integer({ min: 0, max: 8 }), ([ref, q], s, band) => {
        const r = align(ref, q, { ...opts("global", s), band });
        expect(r!.score).toBe(referenceScore(ref, q, subFn(s.match, s.mismatch), s.open, s.extend, "global", band));
      }),
      runs,
    ));
  it("banded score <= unbanded (global) and is self-consistent", () =>
    fc.assert(
      fc.property(seqPair, scoring, fc.integer({ min: 0, max: 8 }), (p, s, band) => {
        const [ref, q] = p;
        const un = align(ref, q, opts("global", s));
        const b = align(ref, q, { ...opts("global", s), band });
        expect(b!.score).toBeLessThanOrEqual(un!.score);
        expect(rescore(ref, q, b!, subFn(s.match, s.mismatch), s.open, s.extend)).toBe(b!.score);
        const l = align(ref, q, { ...opts("local", s), band });
        if (l) expect(rescore(ref, q, l, subFn(s.match, s.mismatch), s.open, s.extend)).toBe(l.score);
        const ul = align(ref, q, opts("local", s));
        expect(l?.score ?? 0).toBeLessThanOrEqual(ul?.score ?? 0);
      }),
      runs,
    ));
});

describe("local mode with empty inputs", () => {
  it("matches oracle", () =>
    fc.assert(
      fc.property(fc.constantFrom("", "A", "ACGT"), fc.constantFrom("", "A", "ACGT"), (a, b) => {
        const r = align(a, b, { mode: "local" });
        const e = referenceScore(a, b, subFn(1, -1), 1, 1, "local");
        if (e <= 0) expect(r).toBeNull();
        else expect(r?.score).toBe(e);
      }),
    ));
});

describe("legacy equivalence", () => {
  it("local mode identical to bioseq 0.1.5", () =>
    prop((ref, q, s) => {
      const old = legacy.align(ref, q, true, [s.match, s.mismatch], [s.open, s.extend]);
      const r = align(ref, q, opts("local", s));
      if (old === null || old === undefined) expect(r).toBeNull();
      else {
        expect(r!.score).toBe(old.score);
        expect(r!.position).toBe(old.position);
        expect(r!.CIGAR).toEqual(Array.from(old.CIGAR));
      }
    }));
  it("global mode identical to legacy, or equal to the oracle", () =>
    prop((ref, q, s) => {
      const old = legacy.align(ref, q, false, [s.match, s.mismatch], [s.open, s.extend]);
      const r = align(ref, q, opts("global", s))!;
      const same = r.score === old.score && JSON.stringify(r.CIGAR) === JSON.stringify(Array.from(old.CIGAR));
      if (!same) expect(r.score).toBe(referenceScore(ref, q, subFn(s.match, s.mismatch), s.open, s.extend, "global"));
    }));
});
