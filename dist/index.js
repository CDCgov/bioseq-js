// src/alphabet.ts
function makeIntArray(length, bitSize, fill) {
  const arr = bitSize <= 8 ? new Int8Array(length) : bitSize <= 16 ? new Int16Array(length) : new Int32Array(length);
  if (fill !== void 0 && fill !== 0) arr.fill(fill);
  return arr;
}
function makeAlphabetMap(str, defaultVal) {
  const aMap = makeIntArray(256, 8, defaultVal);
  const lstr = str.toLowerCase();
  for (let n = 0; n < str.length; n++) {
    aMap[str.charCodeAt(n)] = n;
    aMap[lstr.charCodeAt(n)] = n;
  }
  return aMap;
}
var nt5 = makeAlphabetMap("ACGT", 4);
var nac = makeAlphabetMap("ACGTURYMKWSBDHVN", 16);
var aac = makeAlphabetMap("ARNDCQEGHILKMFPSTWYVBZX", 32);
var aa24 = makeAlphabetMap("ARNDCQEGHILKMFPSTWYVBZX*", 22);
function encode(seq, table = nt5) {
  const s = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i);
    s[i] = table[c < 256 ? c : 0];
  }
  return s;
}
function alphabetSize(table) {
  let max = 0;
  for (let i = 0; i < table.length; i++) if (table[i] > max) max = table[i];
  return max + 1;
}
var COMPLEMENT = {
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
  N: "N"
};
function reverseComplement(seq) {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) {
    const c = seq[i];
    const up = c.toUpperCase();
    const comp = COMPLEMENT[up];
    out += comp === void 0 ? c : c === up ? comp : comp.toLowerCase();
  }
  return out;
}

// src/cigar.ts
var CIGAR_OPS = "MIDNSHP=XB";
var CIGAR_MATCH = 0;
var CIGAR_INS = 1;
var CIGAR_DEL = 2;
var CIGAR_SOFT_CLIP = 4;
function cigar2str(cigar) {
  let s = "";
  for (let k = 0; k < cigar.length; k++) s += (cigar[k] >>> 4).toString() + CIGAR_OPS.charAt(cigar[k] & 15);
  return s;
}
function str2cigar(str) {
  const out = [];
  const re = /(\d+)([MIDNSHP=XB])/gy;
  let m = re.exec(str);
  let consumed = 0;
  while (m !== null) {
    out.push(Number(m[1]) << 4 | CIGAR_OPS.indexOf(m[2]));
    consumed = re.lastIndex;
    m = re.exec(str);
  }
  if (consumed !== str.length) throw new SyntaxError(`Invalid CIGAR string: ${str}`);
  return out;
}
function cigar2gaps(reference, query, start, cigar, fixed) {
  let oq = "";
  let ot = "";
  let lq = 0;
  let lt = start;
  for (let k = 0; k < cigar.length; k++) {
    const op = cigar[k] & 15;
    const len = cigar[k] >>> 4;
    if (op === CIGAR_MATCH || op === 7 || op === 8) {
      oq += query.substr(lq, len);
      ot += reference.substr(lt, len);
      lq += len;
      lt += len;
    } else if (op === CIGAR_INS) {
      if (!fixed) {
        oq += query.substr(lq, len);
        ot += "-".repeat(len);
      }
      lq += len;
    } else if (op === CIGAR_DEL) {
      oq += "-".repeat(len);
      ot += reference.substr(lt, len);
      lt += len;
    } else if (op === CIGAR_SOFT_CLIP) {
      lq += len;
    }
  }
  return [ot, oq];
}
function alignmentStats(reference, query, aln) {
  let matches = 0;
  let mismatches = 0;
  let insertions = 0;
  let deletions = 0;
  let gapOpens = 0;
  let ext = "";
  let lq = 0;
  let lt = aln.position;
  for (const c of aln.CIGAR) {
    const op = c & 15;
    const len = c >>> 4;
    if (op === CIGAR_MATCH || op === 7 || op === 8) {
      let runOp = "";
      let runLen = 0;
      for (let n = 0; n < len; n++) {
        const same = reference[lt + n].toUpperCase() === query[lq + n].toUpperCase();
        if (same) matches++;
        else mismatches++;
        const o = same ? "=" : "X";
        if (o === runOp) runLen++;
        else {
          if (runLen) ext += runLen + runOp;
          runOp = o;
          runLen = 1;
        }
      }
      if (runLen) ext += runLen + runOp;
      lq += len;
      lt += len;
    } else {
      ext += len + CIGAR_OPS.charAt(op);
      if (op === CIGAR_INS) {
        insertions += len;
        gapOpens++;
        lq += len;
      } else if (op === CIGAR_DEL) {
        deletions += len;
        gapOpens++;
        lt += len;
      } else if (op === CIGAR_SOFT_CLIP) {
        lq += len;
      }
    }
  }
  const length = matches + mismatches + insertions + deletions;
  return {
    matches,
    mismatches,
    insertions,
    deletions,
    gapOpens,
    length,
    identity: length ? matches / length : 0,
    extendedCigar: ext
  };
}

// src/align.ts
var NEG_INF = -(1 << 29);
var MAX_TRACEBACK_CELLS = 2 ** 31;
function absPair(gap) {
  if (gap === void 0) return [1, 1];
  if (typeof gap === "number") return [0, Math.abs(gap)];
  return [Math.abs(gap[0]), Math.abs(gap[1])];
}
function isPair(m) {
  return m !== void 0 && m.length === 2 && typeof m[0] === "number" && typeof m[1] === "number";
}
function genScoreMatrix(n, match, mismatch) {
  const b = mismatch > 0 ? -mismatch : mismatch;
  const m = [];
  for (let i = 0; i < n; i++) {
    m[i] = [];
    for (let j = 0; j < n; j++) m[i][j] = i === n - 1 || j === n - 1 ? 0 : i === j ? match : b;
  }
  return m;
}
function checkCodes(s, n, what) {
  for (let i = 0; i < s.length; i++) {
    if (!(s[i] >= 0 && s[i] < n)) {
      throw new RangeError(
        `${what} position ${i} encodes to symbol ${s[i]}, outside the ${n}x${n} score matrix. Check that the encoding table matches the matrix (e.g. use aa24 with BLOSUM62).`
      );
    }
  }
}
function alignWith(reference, query, options = {}) {
  const mode = options.mode ?? "local";
  if (mode !== "local" && mode !== "global" && mode !== "semiglobal") throw new TypeError(`Unknown mode: ${mode}`);
  const isLocal = mode === "local";
  const isSemi = mode === "semiglobal";
  const table = options.table ?? nt5;
  const scoreOnly = options.scoreOnly === true;
  const t = typeof reference === "string" ? encode(reference, table) : reference;
  const q = typeof query === "string" ? encode(query, table) : query;
  const tlen = t.length;
  const qlen = q.length;
  const matrix = isPair(options.matrix) ? genScoreMatrix(alphabetSize(table), options.matrix[0], options.matrix[1]) : options.matrix ?? genScoreMatrix(alphabetSize(table), 1, -1);
  const nsym = matrix.length;
  checkCodes(t, nsym, "Reference");
  checkCodes(q, nsym, "Query");
  const qp = new Int32Array(nsym * qlen);
  for (let a = 0; a < nsym; a++) {
    const row = matrix[a];
    for (let j = 0; j < qlen; j++) qp[a * qlen + j] = row[q[j]];
  }
  const maxLen = qlen > tlen ? qlen : tlen;
  let w = options.band == null || options.band < 0 || isSemi ? maxLen : Math.floor(options.band);
  const lenDiff = tlen > qlen ? tlen - qlen : qlen - tlen;
  if (w < lenDiff) w = lenDiff;
  const [gapo, gape] = absPair(options.gap);
  const gapoe = gapo + gape;
  const H = new Int32Array(qlen + 1);
  const E = new Int32Array(qlen + 1);
  if (!isLocal) {
    H[0] = 0;
    E[0] = w >= 1 ? -(gapoe + gapoe) : NEG_INF;
    for (let j = 1; j <= qlen; j++) {
      H[j] = j > w ? NEG_INF : -(gapoe + gape * (j - 1));
      E[j] = j >= w ? NEG_INF : -(gapoe + gapoe + gape * j);
    }
  }
  const bw = qlen < 2 * w + 1 ? qlen : 2 * w + 1;
  if (!scoreOnly && tlen * bw > MAX_TRACEBACK_CELLS) {
    throw new RangeError(
      `Traceback matrix would need ${tlen * bw} cells; use the band option, scoreOnly, or shorter sequences.`
    );
  }
  const z = scoreOnly ? new Uint8Array(0) : new Uint8Array(tlen * bw);
  let max = 0;
  let endI = -1;
  let endJ = -1;
  let semiBest = NEG_INF;
  for (let i2 = 0; i2 < tlen; i2++) {
    let h1 = 0;
    let f = 0;
    let m = 0;
    let mj = -1;
    const qoff = t[i2] * qlen;
    const beg = i2 > w ? i2 - w : 0;
    const end = i2 + w + 1 < qlen ? i2 + w + 1 : qlen;
    const zoff = i2 * bw - beg;
    if (isSemi) {
      h1 = 0;
      f = -gapoe;
    } else if (!isLocal) {
      h1 = i2 >= w ? NEG_INF : -(gapoe + gape * i2);
      f = i2 >= w ? NEG_INF : -(gapoe + gapoe + gape * i2);
    }
    for (let j = beg; j < end; j++) {
      let e = E[j];
      let h = H[j];
      let d;
      H[j] = h1;
      h += qp[qoff + j];
      d = h > e ? 0 : 1;
      h = h > e ? h : e;
      d = h > f ? d : 2;
      h = h > f ? h : f;
      if (isLocal && !(h > 0)) d = 64;
      h1 = h;
      mj = m > h ? mj : j;
      m = m > h ? m : h;
      h -= gapoe;
      if (isLocal && h < 0) h = 0;
      e -= gape;
      d |= e > h ? 4 : 0;
      e = e > h ? e : h;
      E[j] = e;
      f -= gape;
      d |= f > h ? 32 : 0;
      f = f > h ? f : h;
      if (!scoreOnly) z[zoff + j] = d;
    }
    H[end] = h1;
    E[end] = isLocal ? 0 : NEG_INF;
    if (isLocal) {
      if (m > max) {
        max = m;
        endI = i2;
        endJ = mj;
      }
    } else if (isSemi && h1 > semiBest) {
      semiBest = h1;
      endI = i2;
    }
  }
  let score;
  let refEnd;
  let queryEnd;
  if (isLocal) {
    if (max <= 0) return null;
    score = max;
    refEnd = endI + 1;
    queryEnd = endJ + 1;
  } else if (isSemi) {
    if (qlen === 0) {
      score = 0;
      endI = -1;
    } else if (tlen === 0) {
      score = H[qlen];
    } else {
      score = semiBest;
    }
    refEnd = endI + 1;
    queryEnd = qlen;
  } else {
    score = H[qlen];
    refEnd = tlen;
    queryEnd = qlen;
  }
  if (scoreOnly) return { score, refEnd, queryEnd };
  const cigar = [];
  const push = (op, len) => {
    const last = cigar.length - 1;
    if (last < 0 || op !== (cigar[last] & 15)) cigar.push(len << 4 | op);
    else cigar[last] += len << 4;
  };
  let i;
  let k;
  if (isLocal) {
    i = endI;
    k = endJ;
    if (endJ !== qlen - 1) push(CIGAR_SOFT_CLIP, qlen - 1 - endJ);
  } else if (isSemi) {
    i = endI;
    k = qlen - 1;
  } else {
    i = tlen - 1;
    k = (i + w + 1 < qlen ? i + w + 1 : qlen) - 1;
  }
  let which = 0;
  while (i >= 0 && k >= 0) {
    const tmp = z[i * bw + k - (i > w ? i - w : 0)];
    which = tmp >> (which << 1) & 3;
    if (which === 0 && tmp >> 6) break;
    if (which === 0) which = tmp & 3;
    if (which === 0) {
      push(CIGAR_MATCH, 1);
      i--;
      k--;
    } else if (which === 1) {
      push(CIGAR_DEL, 1);
      i--;
    } else {
      push(CIGAR_INS, 1);
      k--;
    }
  }
  if (isLocal) {
    if (k >= 0) push(CIGAR_SOFT_CLIP, k + 1);
  } else {
    if (i >= 0 && !isSemi) push(CIGAR_DEL, i + 1);
    if (k >= 0) push(CIGAR_INS, k + 1);
  }
  cigar.reverse();
  const refStart = isSemi || isLocal ? i + 1 : 0;
  const queryStart = isLocal ? k + 1 : 0;
  return { score, position: refStart, CIGAR: cigar, refStart, refEnd, queryStart, queryEnd };
}

// src/matrices.ts
var BLOSUM62 = [
  [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
  // A
  [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
  // R
  [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
  // N
  [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
  // D
  [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
  // C
  [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
  // Q
  [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
  // E
  [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
  // G
  [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
  // H
  [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
  // I
  [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
  // L
  [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
  // K
  [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
  // M
  [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
  // F
  [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
  // P
  [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
  // S
  [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
  // T
  [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
  // W
  [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
  // Y
  [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
  // V
  [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4],
  // B
  [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
  // Z
  [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
  // X
  [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]
  // *
];

// src/index.ts
function align(reference, query, optsOrLocal, matrix, gapsc, w, table) {
  if (optsOrLocal !== null && typeof optsOrLocal === "object") return alignWith(reference, query, optsOrLocal);
  return alignWith(reference, query, {
    mode: optsOrLocal === void 0 || optsOrLocal ? "local" : "global",
    matrix,
    gap: gapsc,
    band: w ?? void 0,
    table: table ?? void 0
  });
}
var gen_score_matrix = genScoreMatrix;
function bsg_enc_seq(seq, table = nt5) {
  return Array.from(encode(seq, table));
}
function gen_query_profile(seq, matrix, table) {
  const s = typeof seq === "string" ? bsg_enc_seq(seq, typeof table === "object" && table || nt5) : seq;
  let m;
  if (matrix.length >= 2 && typeof matrix[0] === "number" && typeof matrix[1] === "number") {
    if (table == null) return null;
    const n = typeof table === "number" ? table : table[table.length - 1] + 1;
    m = genScoreMatrix(n, matrix[0], matrix[1]);
  } else {
    m = matrix;
  }
  return m.map((row) => Array.from(s, (c) => row[c]));
}
export {
  BLOSUM62,
  CIGAR_OPS,
  aa24,
  aac,
  align,
  alignWith,
  alignmentStats,
  alphabetSize,
  bsg_enc_seq,
  cigar2gaps,
  cigar2str,
  encode,
  genScoreMatrix,
  gen_query_profile,
  gen_score_matrix,
  makeAlphabetMap,
  makeIntArray,
  nac,
  nt5,
  reverseComplement,
  str2cigar
};
//# sourceMappingURL=index.js.map