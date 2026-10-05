// Score explorer: ribokit's codon lead score (ribokit.score) on one simulated ORF.
(function () {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const C = { f0: "#2a9d8f", off: "#9aa5b1", cds: "#4c78a8", text: "#263238", muted: "#78909c", axis: "#b0bec5" };
  const CUT = 0.0026;   // the p cut of q < 0.05 in orfs.md's 6 real libraries

  // ln k!, for the binomial.
  const LF = [0];
  function lnf(k) {
    while (LF.length <= k) LF.push(LF[LF.length - 1] + Math.log(LF.length));
    return LF[k];
  }
  function binom(k, n, p) {
    if (k < 0 || k > n) return 0;
    if (p <= 0) return k === 0 ? 1 : 0;
    if (p >= 1) return k === n ? 1 : 0;
    return Math.exp(lnf(n) - lnf(k) - lnf(n - k) + k * Math.log(p) + (n - k) * Math.log(1 - p));
  }
  // As score.lead_probability: P(c0 > max(c1, c2)) for n reads, multinomial over the frames with probabilities q.
  function leadProbability(n, q) {
    const r = q[1] + q[2] > 0 ? q[1] / (q[1] + q[2]) : 0;
    let s = 0;
    for (let k0 = 1; k0 <= n; k0++) {
      const m = n - k0, lo = Math.max(m - k0 + 1, 0), hi = Math.min(m, k0 - 1);
      let ahead = 0;
      for (let c1 = lo; c1 <= hi; c1++) ahead += binom(c1, m, r);
      s += binom(k0, n, q[0]) * ahead;
    }
    return s;
  }
  // As score.upper_tail: P(X >= k), X the sum of independent Bernoulli(p).
  function upperTail(p, k) {
    let dist = [1];
    for (const x of p) {
      const next = new Array(dist.length + 1).fill(0);
      dist.forEach((d, i) => { next[i] += d * (1 - x); next[i + 1] += d * x; });
      dist = next;
    }
    return Math.min(dist.slice(k).reduce((a, b) => a + b, 0), 1);
  }

  // Seeded PRNG (mulberry32): the same seed gives the same draw.
  function prng(seed) {
    return function () {
      seed = (seed + 0x6d2b79f5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  const piOf = (p) => [p.pi0, (1 - p.pi0) / 2, (1 - p.pi0) / 2];
  // ORF frame f is CDS frame (f + h) mod 3 for an ORF in frame h of its CDS.
  const inCds = (p, j) => p.overlap > 0 && j >= p.codons - p.overlap;

  // Reads per voting codon and ORF frame. Reads come in clumps of `clump` reads on one nucleotide: one codon, one frame.
  function simulate(p, rand) {
    const counts = Array.from({ length: p.codons }, () => [0, 0, 0]);
    const pick = (q) => { const u = rand(); return u < q[0] ? 0 : u < q[0] + q[1] ? 1 : 2; };
    function drop(total, first, span, q) {
      for (let left = total; left > 0; left -= p.clump) counts[first + Math.floor(rand() * span)][pick(q)] += Math.min(p.clump, left);
    }
    const pi = piOf(p), even = [1 / 3, 1 / 3, 1 / 3];
    drop(p.reads, 0, p.codons, p.translated ? pi : even);
    drop(p.background, 0, p.codons, even);
    if (p.overlap > 0)
      drop(Math.round(p.cdsDensity * p.overlap), p.codons - p.overlap, p.overlap, [0, 1, 2].map((f) => pi[(f + p.cdsFrame) % 3]));
    return counts;
  }

  // The null, as in score.py: each position weighted by the ORF's own density (no frame) plus, inside the host CDS,
  // the CDS's density times 3 pi(that position's frame in the CDS). Densities in reads per nt.
  function nullQ(p, j) {
    const pi = piOf(p), dOrf = (p.reads + p.background) / (3 * p.codons), dCds = p.cdsDensity / 3;
    const w = [0, 1, 2].map((f) => dOrf + (inCds(p, j) ? dCds * 3 * pi[(f + p.cdsFrame) % 3] : 0));
    const s = w[0] + w[1] + w[2];
    return s > 0 ? w.map((x) => x / s) : [1 / 3, 1 / 3, 1 / 3];
  }

  function score(counts, p, independent) {
    const rows = [];
    counts.forEach((c, j) => {
      const n = c[0] + c[1] + c[2];
      if (!n) return;
      const q = nullQ(p, j), multi = leadProbability(n, q);
      rows.push({ j, n, lead: c[0] > Math.max(c[1], c[2]), p: independent ? multi : Math.max(q[0], multi) });
    });
    const ps = rows.map((r) => r.p), leads = rows.filter((r) => r.lead).length;
    const expected = ps.reduce((a, b) => a + b, 0), v = ps.reduce((a, b) => a + b * (1 - b), 0);
    return { rows, leads, expected, z: v > 0 ? (leads - expected) / Math.sqrt(v) : NaN,
             p: upperTail(ps, leads), minP: ps.reduce((a, b) => a * b, 1) };
  }

  function el(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs || {})) e.setAttribute(k, v);
    if (parent) parent.appendChild(e);
    return e;
  }
  function label(parent, x, y, s, attrs) {
    const t = el("text", Object.assign({ x, y, "font-size": 12, fill: C.text, "text-anchor": "middle" }, attrs), parent);
    t.textContent = s;
    return t;
  }
  // Drawing width in SVG units: 720, scaled to the column; on a narrow column (a phone), the column's width in px
  // (at least 300), so that labels keep their size.
  const layoutWidth = (box) => (box.clientWidth < 600 ? Math.max(300, Math.round(box.clientWidth)) : 720);
  const fmtP = (x) => (x >= 0.01 ? x.toFixed(2) : x.toExponential(1));

  function init() {
    const root = document.getElementById("score-explorer");
    if (!root || root.dataset.ready) return;
    root.dataset.ready = "1";
    root.innerHTML = "";

    const controls = document.createElement("div");
    controls.className = "rk-controls";
    const range = (k, text, min, max, step, value, unit) =>
      `<label>${text} <input data-k="${k}" type="range" min="${min}" max="${max}" step="${step}" value="${value}"> ` +
      `<output>${value}</output>${unit ? " " + unit : ""}</label>`;
    controls.innerHTML =
      range("codons", "voting codons", 3, 60, 1, 15) +
      range("reads", "reads from the ORF", 0, 200, 1, 20) +
      '<label><input data-k="translated" type="checkbox" checked> translated</label>' +
      range("pi0", "frame-0 share π<sub>0</sub>", 0.34, 0.96, 0.02, 0.9) +
      range("clump", "reads per clump", 1, 4, 1, 2) +
      range("background", "background reads, no frame", 0, 100, 1, 0) +
      range("overlap", "codons inside a host CDS", 0, 30, 1, 0) +
      range("cdsDensity", "CDS reads per codon", 0, 20, 0.5, 5) +
      '<label>ORF in CDS frame <select data-k="cdsFrame"><option value="1">1</option><option value="2">2</option></select></label>' +
      '<label><input data-k="independent" type="checkbox"> null: independent reads</label>' +
      '<span><button type="button" data-act="draw">new draw</button></span>';
    root.appendChild(controls);
    const q = (k) => controls.querySelector(`[data-k="${k}"]`);

    let W = layoutWidth(controls);
    const svg = el("svg", { role: "img",
                            "aria-label": "Reads per codon and frame of a simulated ORF, and the codons that lead" }, root);
    const readout = document.createElement("p");
    readout.className = "rk-readout";
    root.appendChild(readout);

    let seed = 1;

    function params() {
      const p = {};
      for (const k of ["codons", "reads", "pi0", "clump", "background", "overlap", "cdsDensity"]) {
        p[k] = +q(k).value;
        q(k).nextElementSibling.textContent = k === "pi0" ? p[k].toFixed(2) : String(p[k]);
      }
      p.overlap = Math.min(p.overlap, p.codons);
      p.translated = q("translated").checked;
      p.cdsFrame = +q("cdsFrame").value;
      return p;
    }

    function draw() {
      const p = params(), independent = q("independent").checked;
      const counts = simulate(p, prng(seed));
      const s = score(counts, p, independent);
      root.rk = { params: p, counts, score: s };
      while (svg.firstChild) svg.removeChild(svg.firstChild);

      // Legend: one row; narrow, "✓ leads" on a second row, and the chart 16 lower.
      const narrow = W < 720, [l0, l1, lx, ly] = narrow ? [0, 140, 0, 35] : [40, 190, 330, 19], dy = narrow ? 16 : 0;
      svg.setAttribute("viewBox", `0 0 ${W} ${250 + dy}`);
      el("rect", { x: l0, y: 10, width: 14, height: 10, fill: C.f0, rx: 1 }, svg);
      label(svg, l0 + 20, 19, "frame 0 (the ORF's)", { "text-anchor": "start", "font-size": 11, fill: C.muted });
      el("rect", { x: l1, y: 10, width: 14, height: 10, fill: C.off, rx: 1 }, svg);
      label(svg, l1 + 20, 19, "frame 1 / frame 2", { "text-anchor": "start", "font-size": 11, fill: C.muted });
      label(svg, lx, ly, "✓ leads", { "text-anchor": "start", "font-size": 11, fill: C.f0 });

      const x0 = 40, cw = (W - x0 - 20) / p.codons, bw = Math.min(cw / 3.4, 10), base = 190 + dy, ph = 140;
      const ymax = Math.max(4, ...counts.map((c) => Math.max(...c)));
      if (p.overlap > 0) {
        const xs = x0 + (p.codons - p.overlap) * cw;
        el("rect", { x: xs, y: base - ph - 6, width: p.overlap * cw, height: ph + 6, fill: C.cds, opacity: 0.12 }, svg);
        label(svg, Math.min(xs + (p.overlap * cw) / 2, W - 60), base - ph - 10, "inside the host CDS",
              { "font-size": 11, fill: C.muted });
      }
      el("line", { x1: x0, x2: x0 + p.codons * cw, y1: base, y2: base, stroke: C.axis }, svg);
      counts.forEach((c, j) => {
        const cx = x0 + j * cw + cw / 2;
        c.forEach((n, f) => {
          const h = (n / ymax) * ph;
          el("rect", { x: cx + (f - 1.5) * bw, y: base - h, width: bw - 0.6, height: h, fill: f === 0 ? C.f0 : C.off }, svg);
        });
        if (c[0] > Math.max(c[1], c[2])) label(svg, cx, base + 14, "✓", { "font-size": Math.min(11, cw * 1.25), fill: C.f0 });
        if ((j + 1) % 5 === 0 || p.codons <= 15) label(svg, cx, base + 28, String(j + 1), { "font-size": 10, fill: C.muted });
      });
      label(svg, x0 - 6, base + 28, "codon", { "font-size": 10, fill: C.muted, "text-anchor": "end" });
      label(svg, x0 - 6, base - ph + 4, String(ymax), { "font-size": 10, fill: C.muted, "text-anchor": "end" });
      label(svg, x0 - 6, base, "0", { "font-size": 10, fill: C.muted, "text-anchor": "end" });

      let t = `<b>${s.rows.length}</b> of ${p.codons} codons have reads and vote · <b>${s.leads}</b> lead ` +
              `(${s.expected.toFixed(1)} expected) · z ${isNaN(s.z) ? "–" : s.z.toFixed(1)} · ` +
              `<b>p ${fmtP(s.p)}</b> · min_p ${fmtP(s.minP)}.`;
      if (s.minP > CUT)
        t += ` min_p is above ${CUT}, the p cut of 6 real libraries: with these codons with reads, the ORF could not` +
             " pass even if every vote led.";
      else if (s.p <= CUT) t += ` p is at or below ${CUT}, the p cut of 6 real libraries.`;
      else t += ` min_p is below ${CUT}, the p cut of 6 real libraries, but too few votes led to reach it.`;
      if (independent)
        t += " The independent-reads null takes a codon's reads as separate draws, so a clump that leads looks like" +
             " a rare event: compare p with the box unticked.";
      readout.innerHTML = t;
    }

    controls.addEventListener("input", draw);
    controls.addEventListener("click", (e) => {
      if (e.target.dataset && e.target.dataset.act === "draw") { seed += 1; draw(); }
    });
    draw();
    if (typeof ResizeObserver !== "undefined")
      new ResizeObserver(() => { const w = layoutWidth(controls); if (w !== W) { W = w; draw(); } }).observe(root);
  }

  if (typeof document$ !== "undefined") document$.subscribe(init);
  else document.addEventListener("DOMContentLoaded", init);
})();
