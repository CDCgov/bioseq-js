# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-11-19

### Added

- Options-object API: `align(ref, query, { mode, matrix, gap, band, table, scoreOnly })`.
  The positional signature from 0.1.x still works.
- `semiglobal` mode: the whole query is aligned to the best region of the
  reference, with free end gaps on the reference.
- Richer results: `refStart`, `refEnd`, `queryStart` and `queryEnd` (0-based,
  end-exclusive). `position` is kept as an alias of `refStart`.
- `scoreOnly` option: skips the traceback and uses O(query) memory.
- Protein support: the `aa24` encoding table and the `BLOSUM62` matrix.
- Helpers: `alignmentStats` (matches, mismatches, gaps, identity and an
  `=`/`X` extended CIGAR), `str2cigar`, `reverseComplement`, `encode`,
  `alphabetSize` and `genScoreMatrix`.
- Clear `RangeError`s when the encoding table and score matrix do not match, or
  when the traceback would need too much memory.
- ESM and CommonJS builds, TypeScript type declarations, and a minified browser build.
- Tests (vitest), linting (Biome), a benchmark script and GitHub Actions CI.

### Changed

- Rewritten in TypeScript (`src/`, built with tsup) using typed arrays: roughly
  2.5-3x faster and far less memory (a 10 kb x 10 kb alignment used about 1.6 GB
  of heap before).
- Local alignments still return `null` when no positive-scoring alignment
  exists; this is now documented.
- `gen_score_matrix`, `gen_query_profile`, `bsg_enc_seq` and `makeIntArray`
  are deprecated (still exported).
- README rewritten with a full API reference; demo page updated with mode and
  band controls, alignment statistics and `null` handling.

### Fixed

- The band width parameter `w` was silently ignored in 0.1.x (a typo,
  `t.reference` instead of `t.length`, made it `NaN`). Banding now works: a 10 kb
  global alignment with a band of 50 takes about 23 ms instead of about 440 ms.
- `cigar2str` returned a comma-separated string (`"5M,2I,4M"`); it now returns
  the standard `"5M2I4M"`.
- Global alignment under-penalized an alignment that starts with an insertion
  directly followed by a deletion.
- Band boundaries in global alignment (inherited from `ksw.c`): the edge diagonal
  of the band was excluded, and a leading gap could pass through cells outside
  the band. Banded global scores now exactly match a reference implementation
  restricted to the same band.
- Global and semi-global alignment of an empty reference against a non-empty
  query returned a huge negative score instead of an all-insertion alignment.

### Breaking/behavior changes

- `cigar2str` output format changed from `"5M,2I,4M"` to `"5M2I4M"`.
- The band width is now honored, so results from calls that set `w` can differ
  from 0.1.x.
- Global scores can differ from 0.1.x in the rare case of an insertion
  immediately followed by a deletion at the start of the alignment.
- The package now has an `exports` map and an ESM build, and `main` moved to
  `dist/index.cjs`. Deep imports of other files in the package are no longer
  possible.

## [0.1.x]

Single-file JavaScript implementation ported from Heng Li's `ksw.c`.

