# bioseq-js

A fast, lightweight pairwise sequence aligner: local (Smith-Waterman), global
(Needleman-Wunsch) and semi-global alignment with affine gap penalties, for DNA,
RNA and protein. Works in Node and in the browser, and ships with TypeScript types.

## Use

### Demo

[Click here to try it out.](http://cdcgov.github.io/bioseq-js/demo/)

### Installation

From npm:

```bash
npm install bioseq
```

```javascript
// ESM / TypeScript
import { align } from "bioseq";

// CommonJS
const { align } = require("bioseq");
```

In the browser, include the prebuilt script, which defines a global `bioseq`:

```html
<script src="https://unpkg.com/bioseq/dist/bioseq.min.js"></script>
<!-- or: https://cdn.jsdelivr.net/npm/bioseq/dist/bioseq.min.js -->
<!-- or a local copy of dist/bioseq.min.js -->
```

### Quick Start Example

```javascript
import { align, cigar2str, cigar2gaps, alignmentStats } from "bioseq";

const reference = "ATAGCTAGCTAGCATAAGC";
const query = "AGCTAcCGCAT";

const aln = align(reference, query, { mode: "local" });

console.log("Score:", aln.score);
console.log("CIGAR:", cigar2str(aln.CIGAR));
console.log("Reference:", aln.refStart, "-", aln.refEnd);
console.log("Query:", aln.queryStart, "-", aln.queryEnd);

const [refRow, queryRow] = cigar2gaps(reference, query, aln.refStart, aln.CIGAR);
console.log(refRow);
console.log(queryRow);

const stats = alignmentStats(reference, query, aln);
console.log(stats.matches, "matches,", stats.gapOpens, "gap,", (stats.identity * 100).toFixed(1) + "% identity");
console.log("Extended CIGAR:", stats.extendedCigar);
```

Output:

```
Score: 6
CIGAR: 5M2I4M
Reference: 6 - 15
Query: 0 - 11
AGCTA--GCAT
AGCTAcCGCAT
9 matches, 1 gap, 81.8% identity
Extended CIGAR: 5=2I4=
```

For a more comprehensive example, view the source of [demo/index.html](https://github.com/CDCgov/bioseq-js/blob/master/demo/index.html).

### API

All functions are exported by name (and as properties of the global `bioseq`
in the browser).

#### `align(reference, query, options?)`

Aligns `query` against `reference` and returns an alignment, or `null` (see below).
`reference` and `query` are strings, or arrays of already-encoded symbol codes.

```javascript
align(reference, query, { mode, matrix, gap, band, table, scoreOnly });
```

| Option      | Default     | Description |
| ----------- | ----------- | ----------- |
| `mode`      | `"local"`   | `"local"`, `"global"` or `"semiglobal"` (see below). |
| `matrix`    | `[1, -1]`   | `[match, mismatch]` scores, or a full square substitution matrix indexed by symbol code (e.g. `BLOSUM62`). |
| `gap`       | `[1, 1]`    | Affine gap penalties `[open, extend]`; a gap of length k costs `open + k * extend`. A single number is a linear gap penalty. Signs are ignored, so `[-1, -1]` works too. |
| `band`      | unbanded    | Band width: only cells within this distance of the main diagonal are computed. Automatically widened to at least the length difference of the two sequences. Ignored for `semiglobal`. |
| `table`     | `nt5`       | Encoding table mapping characters to symbol codes. Must match `matrix` (use `aa24` with `BLOSUM62`). |
| `scoreOnly` | `false`     | Skip the traceback and return only `score`, `refEnd` and `queryEnd`, using O(query) memory. |

The legacy positional signature from 0.1.x still works:
`align(reference, query, is_local, matrix, gapsc, w, table)`.

**Modes**

- `local` - Smith-Waterman. The best-scoring pair of substrings; unaligned
  query ends are soft-clipped (`S`).
- `global` - Needleman-Wunsch. Both sequences are aligned end to end.
- `semiglobal` - the whole query is aligned to its best-matching region of the
  reference; gaps at the ends of the reference are free. Typical for mapping a
  read, gene or primer onto a longer sequence.

**Return value**

```typescript
{
  score: number;
  refStart: number;   // 0-based, inclusive
  refEnd: number;     // 0-based, exclusive
  queryStart: number; // 0-based, inclusive
  queryEnd: number;   // 0-based, exclusive
  position: number;   // alias of refStart, kept for compatibility
  CIGAR: number[];    // packed CIGAR, see below
}
```

With `scoreOnly: true` the result contains only `score`, `refEnd` and `queryEnd`.

A **local** alignment returns `null` when no positive-scoring alignment exists
(for example `align("AAAA", "TTTT")`), so check the result before using it.

A `RangeError` is thrown if a sequence contains a symbol outside the score
matrix (usually an encoding table / matrix mismatch), or if the traceback would
need too much memory (use `band` or `scoreOnly`).

#### CIGAR

A CIGAR describes how the query lines up with the reference. `align` returns it
as an array of numbers in the BAM way: each element is `(length << 4) | op`,
where `op` indexes the string `"MIDNSHP=XB"`. Use `cigar2str` to get the
familiar text form.

| Op  | Meaning |
| --- | ------- |
| `M` | alignment match (the residues may be identical or different) |
| `I` | insertion to the reference (extra residues in the query) |
| `D` | deletion from the reference (residues missing from the query) |
| `N` | skipped region of the reference |
| `S` | soft clipping (query residues present but not aligned) |
| `H` | hard clipping (query residues absent) |
| `P` | padding |
| `=` | sequence match |
| `X` | sequence mismatch |

`align` itself only produces `M`, `I`, `D` and `S`; the others are understood by
the helpers below (`alignmentStats` produces `=` and `X`).

#### `cigar2str(cigar)` and `str2cigar(str)`

Convert between the packed form and the string form.

```javascript
cigar2str(aln.CIGAR);     // "5M2I4M"
str2cigar("5M2I4M");      // [80, 33, 64]
```

(In 0.1.x `cigar2str` returned a comma-separated string like `"5M,2I,4M"`.)

#### `cigar2gaps(reference, query, start, cigar, fixed?)`

Returns `[referenceRow, queryRow]`: the two sequences laid out according to the
CIGAR, with `-` for gaps. Soft-clipped query bases are not shown.

- _start_ - the `refStart` of the alignment
- _cigar_ - the `CIGAR` of the alignment
- _fixed_ - if `true`, insertions relative to the reference are dropped so the
  output stays in reference coordinates; useful when aligning many queries to
  the same reference.

#### `alignmentStats(reference, query, aln)`

Takes the same sequences (as strings) and the result of `align`, and returns:

```javascript
{
  matches, mismatches,     // aligned columns that are identical / different
  insertions, deletions,   // number of residues inserted / deleted
  gapOpens,                // number of gap runs
  length,                  // alignment columns (excluding soft clips)
  identity,                // matches / length, between 0 and 1
  extendedCigar,           // CIGAR using = and X instead of M, e.g. "5=2I4="
}
```

#### Encoding tables and matrices

Characters are converted to small integer codes using an encoding table. Each
table is an array indexed by character code.

| Table  | Alphabet |
| ------ | -------- |
| `nt5`  | `ACGT` plus one catch-all symbol for anything else (default) |
| `nac`  | IUPAC nucleotide codes `ACGTURYMKWSBDHVN` |
| `aac`  | amino acids `ARNDCQEGHILKMFPSTWYVBZX` |
| `aa24` | amino acids plus stop: `ARNDCQEGHILKMFPSTWYVBZX*`; matches the `BLOSUM62` matrix |

`BLOSUM62` is the standard 24x24 protein substitution matrix. A protein example:

```javascript
import { align, aa24, BLOSUM62, cigar2str, cigar2gaps } from "bioseq";

const ref = "MKTAYIAKQRQISFVK";
const query = "MKAYIAKQRISFVK";
const aln = align(ref, query, {
  mode: "global",
  table: aa24,
  matrix: BLOSUM62,
  gap: [11, 1],
});

console.log(aln.score, cigar2str(aln.CIGAR)); // 43 2M1D7M1D5M
console.log(cigar2gaps(ref, query, aln.refStart, aln.CIGAR));
// [ 'MKTAYIAKQRQISFVK', 'MK-AYIAKQR-ISFVK' ]
```

Other helpers: `encode(seq, table)` returns a `Uint8Array` of symbol codes;
`alphabetSize(table)` returns the number of symbols in a table;
`genScoreMatrix(n, match, mismatch)` builds an n x n match/mismatch matrix (the
last symbol scores 0 against everything); `makeAlphabetMap(str, defaultVal)`
builds your own table; `reverseComplement(seq)` returns the reverse complement
of a nucleotide string (IUPAC codes aware, case preserved).

```javascript
reverseComplement("ATGCN"); // "NGCAT"
```

#### Deprecated

These 0.1.x helpers are still exported but will be removed in a future
release: `gen_score_matrix` (use `genScoreMatrix`), `bsg_enc_seq` (use
`encode`), `gen_query_profile` (not needed; `align` builds its own profile) and
`makeIntArray`.

### Performance

Alignment takes O(n * m) time for sequences of length n and m. The traceback
needs one byte per cell: n * m bytes, or about n * (2 * band + 1) when `band`
is used. Practical tips:

- For sequences of similar length, set `band` (e.g. a global alignment of two
  10 kb sequences with `band: 50` takes about 23 ms instead of about 440 ms).
- If you only need the score, pass `scoreOnly: true`; no traceback is stored.
- Run `npm run bench` to benchmark on your machine.

## Development

```bash
npm install
npm test          # vitest
npm run lint      # biome
npm run typecheck # tsc --noEmit
npm run build     # tsup -> dist/
npm run bench
```

Source lives in `src/`. The built `dist/` directory is committed because the
GitHub Pages demo loads `dist/bioseq.js`. See [CHANGELOG.md](CHANGELOG.md) for
release notes.

## Theory

`bioseq` implements local and global pairwise alignment with affine gap
penalties. There are two formulations: the Durbin formulation as is
described in his book and the Green formulation as is implemented in [phrap](http://www.phrap.org/phredphrap/phrap.html).
The Durbin formulation is easier to understand, while the Green formulation
is simpler to code and probably faster in practice.

The Durbin formulation is:

    M(i,j) = max{M(i-1,j-1)+S(i,j), E(i-1,j-1), F(i-1,j-1)}

    E(i,j) = max{M(i-1,j)-q-r, F(i-1,j)-q-r, E(i-1,j)-r}

    F(i,j) = max{M(i,j-1)-q-r, F(i,j-1)-r, E(i,j-1)-q-r}

where q is the gap open penalty, r the gap extension penalty and S(i,j) is
the score between the i-th residue in the row sequence and the j-th residue
in the column sequence. Note that the original Durbin formulation disallows
transitions between between E and F states, but we allow them here.

In the Green formulation, we introduce:

    H(i,j) = max{M(i,j), E(i,j), F(i,j)}

The recursion becomes:

    H(i,j) = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}

    E(i,j) = max{H(i-1,j)-q, E(i-1,j)} - r

    F(i,j) = max{H(i,j-1)-q, F(i,j-1)} - r

It is in fact equivalent to the Durbin formulation. In implementation, we
calculate the scores in a different order:

    H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}

    E(i+1,j) = max{H(i,j)-q, E(i,j)} - r

    F(i,j+1) = max{H(i,j)-q, F(i,j)} - r

i.e. at cell (i,j), we compute E for the next row and F for the next column.

## History

In an online conversation, Istvan Albert was complaining he could not find a
good Smith-Waterman implementation in Javascript.
[Li Heng](https://github.com/lh3) thought he could write one over night by
porting [ksw.c](https://github.com/lh3/bwa/blob/master/ksw.c) to javascript. It
took longer than he planned because he found a couple of subtle bugs in ksw.c.
And while he was porting the C code to javascript, he realized that it is not
that difficult to merge local and banded global alignments in one function.
Achieving that also took extra time.

The end product is a fast and lightweight javascript library for affine-gap
local and banded global pairwise alignment. With a modern Javascript engine, it
is not much slower than a non-SSE C implementation.

This repository adapted Li's original algorithm for
[release on NPM](https://www.npmjs.com/package/bioseq). Li did all the hard
work, and we owe him a large debt of gratitude for it.

## Public Domain

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License

The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
Surveillance Platform [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/privacy.html](http://www.cdc.gov/privacy.html).

## Contributing

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page are subject to the [Presidential Records Act](http://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
