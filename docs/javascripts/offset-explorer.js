// Offset explorer: ribokit's CDS-inclusion window search on simulated reads of one length.
(function () {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const C = { cds: "#4c78a8", utr: "#9aa5b1", stop: "#c0392b", in: "#2a9d8f", out: "#e76f51",
              bar: "#cfd8dc", cur: "#f58518", text: "#263238", muted: "#78909c" };

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

  function mulberry32(a) {
    return function () {
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Toy transcript, 0-based: CDS is [CS, CE) without the stop codon; the last sense codon starts at CE - 3.
  const TX_LEN = 420, CS = 60, CE = 360, LEN = 29, MIN_FLANK = 3;
  const PLANTED = [[12, 0.6], [11, 0.15], [13, 0.25]];   // true offset per phase, and its share of reads
  const CENTRES = [];
  for (let c = MIN_FLANK + 1; c < LEN - 3 - MIN_FLANK; c++) CENTRES.push(c);
  const mod3 = (x) => ((x % 3) + 3) % 3;

  function simulate() {
    const rnd = mulberry32(20261002);
    const offset = () => {
      let r = rnd();
      for (const [d, p] of PLANTED) if ((r -= p) < 0) return d;
      return PLANTED[0][0];
    };
    const psites = [];
    const codons = (CE - CS) / 3;
    for (let i = 0; i < 300; i++) psites.push(CS + 3 * Math.floor(rnd() * codons));
    for (let i = 0; i < 20; i++) psites.push(CS, CE - 3);
    for (let i = 0; i < 30; i++) psites.push(rnd() < 0.5 ? Math.floor(rnd() * CS) : CE + Math.floor(rnd() * (TX_LEN - CE)));
    const reads = [];
    for (const p of psites) {
      const pos5 = p - offset();
      if (pos5 >= 0 && pos5 + LEN <= TX_LEN) reads.push({ pos5, phase: mod3(pos5 - CS) });
    }
    return reads;
  }

  // Window centred on c: each phase takes the one of c-1, c, c+1 that puts its P-site on a codon start.
  const windowOffsets = (c) => [0, 1, 2].map((ph) => [c - 1, c, c + 1].find((d) => mod3(d + ph) === 0));
  const inCds = (r, d) => r.pos5 + d >= CS && r.pos5 + d <= CE - 3;
  const spansStart = (r) => CS - r.pos5 >= 0 && CS - r.pos5 < LEN;
  const spansStop = (r) => CE - 3 - r.pos5 >= 0 && CE - 3 - r.pos5 < LEN;
  const countIn = (reads, c) => { const ds = windowOffsets(c); return reads.filter((r) => inCds(r, ds[r.phase])).length; };

  function init() {
    const root = document.getElementById("offset-explorer");
    if (!root || root.dataset.ready) return;
    root.dataset.ready = "1";

    const reads = simulate();
    const startReads = reads.filter(spansStart).sort((a, b) => a.pos5 - b.pos5);
    const stopReads = reads.filter(spansStop).sort((a, b) => a.pos5 - b.pos5);
    const support = startReads.length + stopReads.length;
    const scores = CENTRES.map((c) => countIn(reads, c));
    const supScores = CENTRES.map((c) => countIn(startReads, c) + countIn(stopReads, c));
    const order = CENTRES.map((_, i) => i).sort((i, j) => scores[j] - scores[i] || CENTRES[i] - CENTRES[j]);
    const best = order[0];
    const margin = (scores[order[0]] - scores[order[1]]) / support;
    root.rk = { reads, centres: CENTRES, scores, best: CENTRES[best], margin, offsets: windowOffsets(CENTRES[best]) };

    root.innerHTML = "";
    const controls = document.createElement("div");
    controls.className = "rk-controls";
    controls.innerHTML =
      '<label>window centre <i>c</i> <input type="range" min="' + CENTRES[0] + '" max="' + CENTRES[CENTRES.length - 1] +
      '" value="8" step="1"></label> <button type="button" data-act="prev">◀</button>' +
      '<button type="button" data-act="next">▶</button> <button type="button" data-act="best">best window</button>';
    root.appendChild(controls);
    const slider = controls.querySelector("input");

    const W = 720, panelW = 340, gap = 40, span = 70, ppn = panelW / span;
    const rowH = 5, top = 50;
    const rows = Math.max(startReads.length, stopReads.length);
    const chartTop = top + rows * rowH + 46, chartH = 120;
    const H = chartTop + chartH + 44;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": "Footprints around the start and stop codon, and the score of each offset window" }, root);
    const readout = document.createElement("p");
    readout.className = "rk-readout";
    root.appendChild(readout);

    const panels = [
      { x0: 0, from: CS - 36, reads: startReads, title: "reads spanning the start codon" },
      { x0: panelW + gap, from: CE - 3 - 34, reads: stopReads, title: "reads spanning the last sense codon" },
    ];
    const dyn = [];
    panels.forEach((p, k) => {
      const g = el("g", {}, svg);
      const clip = el("clipPath", { id: `rk-off-clip-${k}` }, g);
      el("rect", { x: p.x0, y: 0, width: panelW, height: H }, clip);
      const X = (pos) => p.x0 + (pos - p.from) * ppn;
      label(g, p.x0 + panelW / 2, 12, p.title, { fill: C.muted });
      const body = el("g", { "clip-path": `url(#rk-off-clip-${k})` }, g);
      el("rect", { x: X(0), y: 26, width: (CS - 0) * ppn, height: 4, fill: C.utr }, body);
      el("rect", { x: X(CE + 3), y: 26, width: (TX_LEN - CE - 3) * ppn, height: 4, fill: C.utr }, body);
      el("rect", { x: X(CS), y: 22, width: (CE - CS) * ppn, height: 12, fill: C.cds }, body);
      el("rect", { x: X(CE), y: 22, width: 3 * ppn, height: 12, fill: C.stop }, body);
      if (k === 0) label(body, X(CS) + 1.5 * ppn, 44, "start", { "font-size": 11, fill: C.muted });
      else label(body, X(CE) + 1.5 * ppn, 44, "stop", { "font-size": 11, fill: C.stop });
      const edge = k === 0 ? CS : CE;
      el("line", { x1: X(edge), x2: X(edge), y1: 20, y2: top + rows * rowH, stroke: "#90a4ae", "stroke-dasharray": "3,3" }, body);
      p.reads.forEach((r, i) => {
        const y = top + i * rowH;
        const bar = el("rect", { x: X(r.pos5), y, width: LEN * ppn, height: rowH - 1, opacity: 0.35 }, body);
        const ps = el("rect", { y, width: 3 * ppn, height: rowH - 1 }, body);
        dyn.push({ r, bar, ps, X });
      });
    });

    // Score chart: support reads whose P-site lands in the CDS, per window centre.
    const cw = W / CENTRES.length;
    const yS = (v) => chartTop + chartH - (v / support) * chartH;
    label(svg, 0, chartTop - 14, "reads spanning the start or last codon whose P-site lands inside the CDS",
          { "text-anchor": "start", fill: C.muted });
    el("line", { x1: 0, x2: W, y1: chartTop + chartH, y2: chartTop + chartH, stroke: "#b0bec5" }, svg);
    el("line", { x1: 0, x2: W, y1: yS(support), y2: yS(support), stroke: "#cfd8dc", "stroke-dasharray": "3,3" }, svg);
    label(svg, W, yS(support) - 4, `all ${support}`, { "text-anchor": "end", "font-size": 11, fill: C.muted });
    const bars = CENTRES.map((c, i) => {
      const b = el("rect", { x: i * cw + 4, width: cw - 8, y: yS(supScores[i]), height: chartTop + chartH - yS(supScores[i]),
                             rx: 2, style: "cursor:pointer" }, svg);
      b.addEventListener("click", () => set(c));
      label(svg, i * cw + cw / 2, chartTop + chartH + 16, String(c), { "font-size": 11, fill: C.muted });
      return b;
    });
    label(svg, best * cw + cw / 2, yS(supScores[best]) - 6, "★", { fill: C.cur, "font-size": 14 });
    label(svg, W / 2, chartTop + chartH + 36, "window centre c", { fill: C.muted });

    function set(c) {
      c = Math.max(CENTRES[0], Math.min(CENTRES[CENTRES.length - 1], c));
      slider.value = c;
      const ds = windowOffsets(c);
      for (const { r, bar, ps, X } of dyn) {
        const d = ds[r.phase], col = inCds(r, d) ? C.in : C.out;
        bar.setAttribute("fill", col);
        ps.setAttribute("x", X(r.pos5 + d));
        ps.setAttribute("fill", col);
      }
      const i = CENTRES.indexOf(c);
      bars.forEach((b, j) => b.setAttribute("fill", j === i ? C.cur : C.bar));
      const isBest = i === best;
      readout.innerHTML =
        `<b>c = ${c}</b> → offsets: phase 0 → <b>${ds[0]}</b>, phase 1 → <b>${ds[1]}</b>, phase 2 → <b>${ds[2]}</b> nt. ` +
        `<span style="color:${C.in}">${supScores[i]}</span> of ${support} reads that span the start or last codon get their ` +
        `P-site inside the CDS (<span style="color:${C.out}">${support - supScores[i]}</span> outside). ` +
        (isBest ? `<b>This is the best window</b> (★); margin ${margin.toFixed(2)}.`
                : `The best window is c = ${CENTRES[best]} (★).`);
    }
    slider.addEventListener("input", () => set(+slider.value));
    controls.addEventListener("click", (e) => {
      const act = e.target.dataset && e.target.dataset.act;
      if (act === "prev") set(+slider.value - 1);
      if (act === "next") set(+slider.value + 1);
      if (act === "best") set(CENTRES[best]);
    });
    set(+slider.value);
  }

  if (typeof document$ !== "undefined") document$.subscribe(init);
  else document.addEventListener("DOMContentLoaded", init);
})();
