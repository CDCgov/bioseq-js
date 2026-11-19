/**
 * Naive O(n*m) Gotoh oracle. Computes ONLY the optimal score.
 *
 * Three states per cell (i = reference residues consumed, j = query residues consumed):
 *   M[i][j]  alignment ends with ref[i-1] aligned to query[j-1]
 *   X[i][j]  ends with a deletion (reference residue against a gap)
 *   Y[i][j]  ends with an insertion (query residue against a gap)
 * A gap of length k costs open + k * extend. A gap may be opened directly after any
 * state, including the opposite gap type (a deletion may sit next to an insertion).
 *
 * `pred(i, j)` is the best score of anything that may precede the next column when
 * we stand at cell (i, j): any of the three states, or a free start where the mode allows it.
 */
export type Mode = "local" | "global" | "semiglobal";

const NEG = Number.NEGATIVE_INFINITY;

export function referenceScore(
  ref: string,
  query: string,
  sub: (a: string, b: string) => number,
  open: number,
  extend: number,
  mode: Mode,
  band?: number,
): number {
  const n = ref.length;
  const m = query.length;
  const grid = () => Array.from({ length: n + 1 }, () => new Array<number>(m + 1).fill(NEG));
  const M = grid();
  const X = grid();
  const Y = grid();

  const freeStart = (i: number, j: number): boolean =>
    mode === "local" || (mode === "global" && i === 0 && j === 0) || (mode === "semiglobal" && j === 0);
  const pred = (i: number, j: number): number => Math.max(M[i][j], X[i][j], Y[i][j], freeStart(i, j) ? 0 : NEG);

  // Banding (as in the library): only cells with |i - j| <= w are reachable, w >= |n - m|.
  const w = band === undefined ? Number.POSITIVE_INFINITY : Math.max(band, Math.abs(n - m));
  for (let i = 0; i <= n; i++) {
    for (let j = 0; j <= m; j++) {
      if (Math.abs(i - j) > w) continue;
      if (i > 0 && j > 0) M[i][j] = pred(i - 1, j - 1) + sub(ref[i - 1], query[j - 1]);
      if (i > 0) X[i][j] = Math.max(X[i - 1][j] - extend, pred(i - 1, j) - open - extend);
      if (j > 0) Y[i][j] = Math.max(Y[i][j - 1] - extend, pred(i, j - 1) - open - extend);
    }
  }

  if (mode === "global") return pred(n, m);
  if (mode === "semiglobal") {
    let s = NEG;
    for (let i = 0; i <= n; i++) s = Math.max(s, pred(i, m));
    return s;
  }
  // local: best cell anywhere; 0 means "no positive alignment"
  let s = 0;
  for (let i = 0; i <= n; i++) for (let j = 0; j <= m; j++) s = Math.max(s, M[i][j], X[i][j], Y[i][j]);
  return s;
}
