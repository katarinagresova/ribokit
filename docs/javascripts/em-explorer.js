// EM explorer: ribokit's EM on two components that share reads. #em-explorer (method.md): two CDSs, with and
// without the 1/L term. #em-explorer-frame (orfs.md): a uoORF over its CDS's start, with and without the frame term.
(function () {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const C = { A: "#4c78a8", B: "#f58518", AB: "#8e6bbf", uoORF: "#e76f51", CDS: "#4c78a8", uoORFCDS: "#9a6f7d",
              text: "#263238", muted: "#78909c", grid: "#eceff1" };
  const TOL = 1e-3, MAX_ITER = 100000;
  const SCENARIOS = {
    nested: { name: "B lies inside A", regions: [{ len: 300, members: [0] }, { len: 600, members: [0, 1] }] },
    ends: { name: "A and B each have a part of their own",
            regions: [{ len: 300, members: [0] }, { len: 600, members: [0, 1] }, { len: 150, members: [1] }] },
    tie: { name: "A and B are identical", regions: [{ len: 600, members: [0, 1] }] },
  };
  // A uoORF in frame 1 of its CDS: its own part, the overlap, the CDS's own part.
  const FRAME = { regions: [{ len: 90, members: [0] }, { len: 90, members: [0, 1] }, { len: 180, members: [1] }] };

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
  const fmt = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(1));

  // Expected reads per region = length x summed density of the CDSs that contain it.
  function model(sc, rho) {
    const L = [0, 0], unique = [0, 0], classes = [];
    for (const r of sc.regions) {
      for (const t of r.members) L[t] += r.len;
      const n = r.len * r.members.reduce((s, t) => s + rho[t], 0);
      if (r.members.length === 1) unique[r.members[0]] += n;
      else classes.push({ members: r.members, count: n });
    }
    return { L, unique, classes, truth: [rho[0] * L[0], rho[1] * L[1]] };
  }

  // The uoORF case. The overlap's reads form 3 classes by CDS frame g; g is the uoORF's frame g - 1. Each ORF's
  // weight in a class is 3 pi(its frame), as in orfs.py's weighted classes (1 without the frame term). pi = (pi0,
  // the rest split evenly). `peak` reads of the CDS's start peak fall in the uoORF's own part (the U4 case).
  function frameModel(rho, pi0, peak, useFrame) {
    const pi = [pi0, (1 - pi0) / 2, (1 - pi0) / 2];
    const [own, ov, cds] = FRAME.regions.map((r) => r.len);
    const classes = [0, 1, 2].map((g) => {
      const u = (g + 2) % 3;
      return { members: [0, 1], frame: g, count: ov * (rho[0] * pi[u] + rho[1] * pi[g]),
               f: useFrame ? [3 * pi[u], 3 * pi[g]] : [1, 1] };
    });
    return { L: [own + ov, ov + cds], unique: [rho[0] * own + peak, rho[1] * cds], classes,
             truth: [rho[0] * (own + ov), rho[1] * (ov + cds) + peak] };
  }

  // Plain EM step, as in ribokit.quant.em: even start, then alpha_t = u_t + sum_c n_c w_t / sum_{s in c} w_s, w = alpha / L
  // (times the class weight f, if the class has one).
  // `ribokit quant` accelerates this with SQUAREM (method.md, "Acceleration"); shown here unaccelerated, step by step.
  function em(m, useLength) {
    const len = useLength ? m.L : [1, 1];
    let alpha = m.unique.slice();
    for (const c of m.classes) for (const t of c.members) alpha[t] += c.count / c.members.length;
    const trace = [alpha], changes = [Infinity];
    for (let it = 1; it <= MAX_ITER; it++) {
      const w = alpha.map((a, t) => a / len[t]);
      const next = m.unique.slice();
      for (const c of m.classes) {
        const f = c.f || c.members.map(() => 1);
        const denom = c.members.reduce((s, t, i) => s + w[t] * f[i], 0);
        if (denom > 0) c.members.forEach((t, i) => { next[t] += (c.count * w[t] * f[i]) / denom; });
      }
      const change = Math.max(...next.map((a, t) => Math.abs(a - alpha[t])));
      alpha = next;
      trace.push(alpha);
      changes.push(change);
      if (change < TOL) break;
    }
    return { trace, changes };
  }

  function setup(root, frame) {
    if (root.dataset.ready) return;
    root.dataset.ready = "1";
    root.innerHTML = "";
    const NAMES = frame ? ["uoORF", "CDS"] : ["A", "B"];

    const controls = document.createElement("div");
    controls.className = "rk-controls";
    const buttons = '<span><button type="button" data-act="reset">reset</button> <button type="button" data-act="step">step</button> ' +
      '<button type="button" data-act="run">run</button></span>';
    controls.innerHTML = frame ?
      '<label>density of the uoORF <input data-k="rA" type="range" min="0" max="3" step="0.1" value="0"> <output>0.0</output> reads/nt</label>' +
      '<label>density of the CDS <input data-k="rB" type="range" min="0" max="3" step="0.1" value="2"> <output>2.0</output> reads/nt</label>' +
      '<label>frame-0 share π<sub>0</sub> <input data-k="pi" type="range" min="0.34" max="0.96" step="0.02" value="0.9"> <output>0.90</output></label>' +
      '<label>CDS start-peak reads in the uoORF\'s own part <input data-k="peak" type="range" min="0" max="120" step="5" value="40"> <output>40</output></label>' +
      '<label><input data-k="frame" type="checkbox" checked> frame term 3π</label>' + buttons :
      '<label>case <select data-k="sc">' +
      Object.entries(SCENARIOS).map(([k, s]) => `<option value="${k}">${s.name}</option>`).join("") + "</select></label>" +
      '<label>density of A <input data-k="rA" type="range" min="0" max="3" step="0.1" value="1"> <output>1.0</output> reads/nt</label>' +
      '<label>density of B <input data-k="rB" type="range" min="0" max="3" step="0.1" value="2"> <output>2.0</output> reads/nt</label>' +
      '<label><input data-k="len" type="checkbox" checked> length term 1/L</label>' + buttons;
    root.appendChild(controls);
    const q = (k) => controls.querySelector(`[data-k="${k}"]`);

    let W = layoutWidth(controls);
    const svg = el("svg", { role: "img",
                            "aria-label": frame ? "A uoORF over its CDS's start, their reads, and the EM estimate over steps" :
                              "Two CDSs, their reads, and the EM estimate over steps" }, root);
    const readout = document.createElement("p");
    readout.className = "rk-readout";
    root.appendChild(readout);

    let m, res, k = 0, timer = null;

    function recompute() {
      const rho = [+q("rA").value, +q("rB").value];
      q("rA").nextElementSibling.textContent = rho[0].toFixed(1);
      q("rB").nextElementSibling.textContent = rho[1].toFixed(1);
      if (frame) {
        q("pi").nextElementSibling.textContent = (+q("pi").value).toFixed(2);
        q("peak").nextElementSibling.textContent = q("peak").value;
        m = frameModel(rho, +q("pi").value, +q("peak").value, q("frame").checked);
        res = em(m, true);
      } else {
        m = model(SCENARIOS[q("sc").value], rho);
        res = em(m, q("len").checked);
      }
      root.rk = { model: m, steps: res.trace.length - 1, alpha: res.trace[res.trace.length - 1] };
    }

    function draw() {
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const sc = frame ? FRAME : SCENARIOS[q("sc").value];
      const K = res.trace.length - 1, alpha = res.trace[k];
      const narrow = W < 720;

      // Structure: the CDSs as rows, regions coloured by which CDSs contain them.
      const total = sc.regions.reduce((s, r) => s + r.len, 0), x0 = frame ? 64 : 40;
      const sx = ((narrow ? W - 10 : 680) - x0) / total;
      let x = x0, line = 0;
      // Wide: each label under its region. Narrow: one line per label, in region order.
      const note = (xm, y, s, attrs, indent) => narrow ?
        label(svg, indent ? 14 : 2, 92 + 16 * line++, s, Object.assign({ "text-anchor": "start" }, attrs)) :
        label(svg, xm, y, s, attrs);
      NAMES.forEach((n, t) => label(svg, frame ? 30 : 14, 31 + 30 * t, n, { "font-weight": "bold", fill: C[n] }));
      sc.regions.forEach((r, i) => {
        const key = r.members.map((t) => NAMES[t]).join("");
        const w = r.len * sx;
        for (const t of r.members) el("rect", { x, y: 18 + 30 * t, width: w - 2, height: 18, rx: 3, fill: C[key] }, svg);
        const n = frame ? [m.unique[0], m.classes.reduce((s, c) => s + c.count, 0), m.unique[1]][i] :
          r.len * r.members.reduce((s, t) => s + [+q("rA").value, +q("rB").value][t], 0);
        label(svg, x + w / 2, 31 + 30 * (r.members.length === 1 ? r.members[0] : 0), `${r.len} nt`,
              { fill: "white", "font-size": 11 });
        const to = frame ? (r.members.length > 1 ? "both" : NAMES[r.members[0]]) :
          `{${r.members.map((t) => NAMES[t]).join(", ")}}`;
        note(x + w / 2, 92, `${fmt(n)} reads → ${to}`, { fill: C[key], "font-size": 12 });
        if (frame && i === 0 && +q("peak").value > 0)
          note(x + w / 2, 108, `incl. ${q("peak").value} start-peak reads`, { fill: C.muted, "font-size": 11 }, true);
        if (frame && i === 1)
          note(x + w / 2, 124, `CDS frame 0 / 1 / 2: ${m.classes.map((c) => fmt(c.count)).join(" / ")}`,
               { fill: C.muted, "font-size": 11 }, true);
        x += w;
      });

      // Bars: truth (outline) vs current estimate (filled). Wide: the trace to their right. Narrow: the trace below,
      // under room for the most label lines (5 with the frame term, 3 without).
      const top = narrow ? 92 + 16 * (frame ? 4 : 2) + 36 : frame ? 150 : 128, ph = narrow ? 120 : 150, base = top + ph;
      const bx0 = narrow ? (W - 250) / 2 : 0;
      const ymax = res.trace.reduce((mx, a) => Math.max(mx, a[0], a[1]), Math.max(1, ...m.truth)) * 1.1;
      const Y = (v, b = base) => b - (v / ymax) * ph;
      el("line", { x1: bx0, x2: bx0 + 250, y1: base, y2: base, stroke: "#b0bec5" }, svg);
      NAMES.forEach((n, t) => {
        const bx = bx0 + 30 + t * 115;
        el("rect", { x: bx, y: Y(alpha[t]), width: 40, height: base - Y(alpha[t]), fill: C[n], rx: 2 }, svg);
        el("rect", { x: bx + 44, y: Y(m.truth[t]), width: 40, height: base - Y(m.truth[t]), fill: "none",
                     stroke: C[n], "stroke-dasharray": "4,3", "stroke-width": 1.5, rx: 2 }, svg);
        label(svg, bx + 20, Y(alpha[t]) - 5, fmt(alpha[t]), { "font-size": 11 });
        label(svg, bx + 64, Y(m.truth[t]) - 5, fmt(m.truth[t]), { "font-size": 11, fill: C.muted });
        label(svg, bx + 20, base + 15, frame ? "EM" : `${n} EM`, { "font-size": 11 });
        label(svg, bx + 64, base + 15, frame ? "true" : `${n} true`, { "font-size": 11, fill: C.muted });
        if (frame) label(svg, bx + 42, base + 30, n, { "font-size": 11, "font-weight": "bold", fill: C[n] });
      });
      label(svg, bx0 + 125, top - 10, "expected reads", { fill: C.muted });

      // Trace: estimate per step, truth dashed.
      const ttop = narrow ? base + 66 : top, tbase = ttop + ph;
      const px0 = narrow ? 10 : 300, pw = W - px0 - 10, X = (i) => px0 + (i / Math.max(K, 10)) * pw;
      svg.setAttribute("viewBox", `0 0 ${W} ${tbase + 52}`);
      el("line", { x1: px0, x2: px0 + pw, y1: tbase, y2: tbase, stroke: "#b0bec5" }, svg);
      el("line", { x1: px0, x2: px0, y1: ttop, y2: tbase, stroke: "#b0bec5" }, svg);
      NAMES.forEach((n, t) => {
        el("line", { x1: px0, x2: px0 + pw, y1: Y(m.truth[t], tbase), y2: Y(m.truth[t], tbase), stroke: C[n],
                     "stroke-dasharray": "4,3", opacity: 0.6 }, svg);
        const pts = res.trace.slice(0, k + 1).map((a, i) => `${X(i).toFixed(1)},${Y(a[t], tbase).toFixed(1)}`).join(" ");
        el("polyline", { points: pts, fill: "none", stroke: C[n], "stroke-width": 2 }, svg);
        el("circle", { cx: X(k), cy: Y(alpha[t], tbase), r: 3.5, fill: C[n] }, svg);
      });
      label(svg, px0 + pw / 2, ttop - 10, "estimate per step (dashed: truth)", { fill: C.muted });
      label(svg, px0, tbase + 15, "0", { "font-size": 11, fill: C.muted });
      label(svg, px0 + pw, tbase + 15, String(Math.max(K, 10)), { "font-size": 11, fill: C.muted, "text-anchor": "end" });
      label(svg, px0 + pw / 2, tbase + 30, "EM step", { fill: C.muted, "font-size": 11 });

      // Readout.
      const done = k === K;
      let s = `<b>EM step ${k}</b>` + (k === 0 ? " (even start: each shared read split equally)" :
        `, largest change ${res.changes[k] < 1e-6 ? "0" : res.changes[k].toPrecision(2)} reads`) +
        (done ? ` · <b>converged</b> after ${K} EM steps (no estimate moves by ${TOL} reads or more).` : ".");
      if (done && frame) {
        const peak = +q("peak").value;
        s += q("frame").checked ?
          " With the frame term, each frame class of the overlap is split by density times 3π of its frame in each ORF." :
          " Without the frame term, the overlap is split by density alone.";
        if (peak === 0) s += " With the reads spread evenly, density alone already splits it right: the frame term" +
                            " matters where density misleads.";
        else if (q("frame").checked) s += " The start-peak reads still count for the uoORF, and so do most of the CDS's" +
                                         " frame-1 reads, which are in the uoORF's frame 0: a known limit.";
        else s += " The start-peak reads make the uoORF look denser than it is, so it also takes a share of the CDS's" +
                  " reads in the overlap.";
      } else if (done) {
        if (q("sc").value === "tie") {
          s += " No read tells A and B apart, so the EM keeps its even starting split whatever the true densities:" +
               " ribokit lists such CDSs in <code>ties.tsv</code>.";
        } else if (!q("len").checked) {
          s += " Without 1/L, the shared reads end up split in proportion to each CDS's own reads, however long" +
               " its own part is, not by density.";
        } else {
          s += " With 1/L, each shared read is split by the CDSs' densities (reads per nt), which recovers the truth.";
        }
      }
      readout.innerHTML = s;
    }

    function stop() { if (timer) { cancelAnimationFrame(timer); timer = null; } }
    function run() {
      stop();
      const K = res.trace.length - 1;
      if (k >= K) k = 0;
      const stepSize = Math.max(1, Math.ceil((K - k) / 60));
      const tick = () => {
        k = Math.min(K, k + stepSize);
        draw();
        timer = k < K ? requestAnimationFrame(tick) : null;
      };
      timer = requestAnimationFrame(tick);
    }

    controls.addEventListener("input", () => { stop(); recompute(); k = 0; run(); });
    controls.addEventListener("click", (e) => {
      const act = e.target.dataset && e.target.dataset.act;
      if (!act) return;
      stop();
      if (act === "reset") k = 0;
      if (act === "step") k = Math.min(res.trace.length - 1, k + 1);
      if (act === "run") return run();
      draw();
    });
    recompute();
    k = res.trace.length - 1;
    draw();
    if (typeof ResizeObserver !== "undefined")
      new ResizeObserver(() => { const w = layoutWidth(controls); if (w !== W) { W = w; draw(); } }).observe(root);
  }

  function init() {
    const a = document.getElementById("em-explorer"), b = document.getElementById("em-explorer-frame");
    if (a) setup(a, false);
    if (b) setup(b, true);
  }

  if (typeof document$ !== "undefined") document$.subscribe(init);
  else document.addEventListener("DOMContentLoaded", init);
})();
