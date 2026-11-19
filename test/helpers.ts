import { createRequire } from "node:module";
import fc from "fast-check";
import type { Alignment } from "../src/index";

const require = createRequire(import.meta.url);
// biome-ignore lint/suspicious/noExplicitAny: legacy CommonJS oracle
export const legacy: any = require("./legacy/bioseq-0.1.5.cjs");

export const subFn = (match: number, mismatch: number) => (a: string, b: string) =>
  a.toUpperCase() === b.toUpperCase() ? match : mismatch;

/** Re-score the alignment described by a result against the sequences. Throws on inconsistency. */
export function rescore(
  ref: string,
  query: string,
  aln: Alignment,
  sub: (a: string, b: string) => number,
  open: number,
  extend: number,
): number {
  let r = aln.refStart;
  let q = aln.queryStart;
  let score = 0;
  for (const c of aln.CIGAR) {
    const op = "MIDNSHP=XB"[c & 0xf];
    const len = c >>> 4;
    if (op === "M") {
      for (let k = 0; k < len; k++) score += sub(ref[r + k], query[q + k]);
      r += len;
      q += len;
    } else if (op === "I") {
      score -= open + extend * len;
      q += len;
    } else if (op === "D") {
      score -= open + extend * len;
      r += len;
    } else if (op === "S") {
      // soft clips lie outside [queryStart, queryEnd): nothing to score
    } else throw new Error(`unexpected op ${op}`);
  }
  if (r !== aln.refEnd) throw new Error(`ref end mismatch ${r} vs ${aln.refEnd}`);
  return score;
}

export const opSum = (aln: Alignment, ops: string): number =>
  aln.CIGAR.filter((c) => ops.includes("MIDNSHP=XB"[c & 0xf])).reduce((a, c) => a + (c >>> 4), 0);

const dna = (max: number) => fc.string({ unit: fc.constantFrom(..."ACGT"), maxLength: max });

/** A pair of (reference, query): independent random strings, or a mutated copy. */
export const seqPair = fc.oneof(
  fc.tuple(dna(40), dna(40)),
  fc
    .tuple(
      dna(40),
      fc.array(fc.tuple(fc.nat(60), fc.constantFrom("sub", "ins", "del"), fc.constantFrom(..."ACGT")), {
        maxLength: 6,
      }),
      fc.nat(10),
      fc.nat(10),
    )
    .map(([base, muts, cutL, cutR]) => {
      let q = base.slice(Math.min(cutL, base.length), Math.max(0, base.length - cutR));
      for (const [pos, kind, ch] of muts) {
        const p = q.length ? pos % (q.length + 1) : 0;
        if (kind === "sub") q = q.slice(0, p) + ch + q.slice(p + 1);
        else if (kind === "ins") q = q.slice(0, p) + ch + q.slice(p);
        else q = q.slice(0, p) + q.slice(p + 1);
      }
      return [base, q] as [string, string];
    }),
);

export const scoring = fc.record({
  match: fc.integer({ min: 1, max: 5 }),
  mismatch: fc.integer({ min: -5, max: -1 }),
  open: fc.integer({ min: 0, max: 6 }),
  extend: fc.integer({ min: 1, max: 3 }),
});
