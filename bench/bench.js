// Benchmarks for bioseq: new implementation (dist/index.js) vs the legacy 0.1.5 build.
// Run with: npm run bench   (requires `npm run build` first so dist/ is current)

import { createRequire } from "node:module";
import { Bench } from "tinybench";
import * as bioseq from "../dist/index.js";

const require = createRequire(import.meta.url);
const legacy = require("../test/legacy/bioseq-0.1.5.cjs");

const SEED = 20260519;
const DNA = "ACGT";
const AMINO = "ARNDCQEGHILKMFPSTWYV";

/** Small seeded PRNG (mulberry32) so every run benchmarks identical inputs. */
function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Random sequence over `alphabet`. */
function randomSeq(rand, length, alphabet) {
  let out = "";
  for (let i = 0; i < length; i++) out += alphabet[Math.floor(rand() * alphabet.length)];
  return out;
}

/** Copy of `seq` with ~`rate` of positions substituted by a different symbol. */
function mutate(rand, seq, rate, alphabet) {
  let out = "";
  for (let i = 0; i < seq.length; i++) {
    if (rand() < rate) {
      let c;
      do {
        c = alphabet[Math.floor(rand() * alphabet.length)];
      } while (c === seq[i]);
      out += c;
    } else {
      out += seq[i];
    }
  }
  return out;
}

// Reference/query pairs, generated once so all benchmarks share inputs.
const rand = mulberry32(SEED);
const dna = {};
for (const len of [1000, 5000, 10000]) {
  const ref = randomSeq(rand, len, DNA);
  dna[len] = { ref, qry: mutate(rand, ref, 0.05, DNA) };
}
const protRef = randomSeq(rand, 500, AMINO);
const prot = { ref: protRef, qry: mutate(rand, protRef, 0.05, AMINO) };

// Scoring shared across DNA cases: match +2, mismatch -3, affine gap open 5 / extend 2.
const DNA_MATRIX = [2, -3];
const DNA_GAP = [5, 2];
const LEGACY_GAP = [-5, -2]; // legacy 0.1.x expects negative gap penalties

const bench = new Bench({ warmup: true });

// Quick-running cases: a generous time budget for stable numbers.
const FAST = { time: 1000 };
// Expensive cases: a fixed, small number of iterations to keep total runtime bounded.
const SLOW = { time: 1, iterations: 5 };

bench
  .add("legacy 0.1.5 local, 1kb", () => legacy.align(dna[1000].ref, dna[1000].qry, true, DNA_MATRIX, LEGACY_GAP), FAST)
  .add(
    "new local, 1kb",
    () => bioseq.align(dna[1000].ref, dna[1000].qry, { mode: "local", matrix: DNA_MATRIX, gap: DNA_GAP }),
    FAST,
  )
  .add("legacy 0.1.5 local, 5kb", () => legacy.align(dna[5000].ref, dna[5000].qry, true, DNA_MATRIX, LEGACY_GAP), SLOW)
  .add(
    "new local, 5kb",
    () => bioseq.align(dna[5000].ref, dna[5000].qry, { mode: "local", matrix: DNA_MATRIX, gap: DNA_GAP }),
    SLOW,
  )
  .add(
    "new global, band 50, 10kb",
    () =>
      bioseq.align(dna[10000].ref, dna[10000].qry, {
        mode: "global",
        matrix: DNA_MATRIX,
        gap: DNA_GAP,
        band: 50,
      }),
    SLOW,
  )
  .add(
    "new local scoreOnly, 5kb",
    () =>
      bioseq.align(dna[5000].ref, dna[5000].qry, {
        mode: "local",
        matrix: DNA_MATRIX,
        gap: DNA_GAP,
        scoreOnly: true,
      }),
    SLOW,
  )
  .add(
    "new protein local (aa24 + BLOSUM62, gap 10/1), 500aa",
    () =>
      bioseq.align(prot.ref, prot.qry, {
        mode: "local",
        table: bioseq.aa24,
        matrix: bioseq.BLOSUM62,
        gap: [10, 1],
      }),
    FAST,
  );

console.log(`bioseq benchmarks (seed ${SEED}, Node ${process.version})`);
await bench.run();
console.table(bench.table());
