// EM explorer: ribokit's EM on two CDSs that share reads, with and without the 1/L term.
(function () {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const C = { A: "#4c78a8", B: "#f58518", AB: "#8e6bbf", text: "#263238", muted: "#78909c", grid: "#eceff1" };
  const TOL = 1e-3, MAX_ITER = 100000;
  const NAMES = ["A", "B"];

  const SCENARIOS = {
    nested: { name: "B lies inside A", regions: [{ len: 300, members: [0] }, { len: 600, members: [0, 1] }] },
    ends: { name: "A and B each have a part of their own",
            regions: [{ len: 300, members: [0] }, { len: 600, members: [0, 1] }, { len: 150, members: [1] }] },
    tie: { name: "A and B are identical", regions: [{ len: 600, members: [0, 1] }] },
  };

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

  // Same steps as ribokit.quant.em: even start, then alpha_t = u_t + sum_c n_c w_t / sum_{s in c} w_s, w = alpha / L.
  function em(m, useLength) {
    const len = useLength ? m.L : [1, 1];
    let alpha = m.unique.slice();
    for (const c of m.classes) for (const t of c.members) alpha[t] += c.count / c.members.length;
    const trace = [alpha], changes = [Infinity];
    for (let it = 1; it <= MAX_ITER; it++) {
      const w = alpha.map((a, t) => a / len[t]);
      const next = m.unique.slice();
      for (const c of m.classes) {
        const denom = c.members.reduce((s, t) => s + w[t], 0);
        if (denom > 0) for (const t of c.members) next[t] += (c.count * w[t]) / denom;
      }
      const change = Math.max(...next.map((a, t) => Math.abs(a - alpha[t])));
      alpha = next;
      trace.push(alpha);
      changes.push(change);
      if (change < TOL) break;
    }
    return { trace, changes };
  }

  function init() {
    const root = document.getElementById("em-explorer");
    if (!root || root.dataset.ready) return;
    root.dataset.ready = "1";
    root.innerHTML = "";

    const controls = document.createElement("div");
    controls.className = "rk-controls";
    controls.innerHTML =
      '<label>case <select data-k="sc">' +
      Object.entries(SCENARIOS).map(([k, s]) => `<option value="${k}">${s.name}</option>`).join("") + "</select></label>" +
      '<label>density of A <input data-k="rA" type="range" min="0" max="3" step="0.1" value="1"> <output>1.0</output> reads/nt</label>' +
      '<label>density of B <input data-k="rB" type="range" min="0" max="3" step="0.1" value="2"> <output>2.0</output> reads/nt</label>' +
      '<label><input data-k="len" type="checkbox" checked> length term 1/L</label>' +
      '<span><button type="button" data-act="reset">reset</button> <button type="button" data-act="step">step</button> ' +
      '<button type="button" data-act="run">run</button></span>';
    root.appendChild(controls);
    const q = (k) => controls.querySelector(`[data-k="${k}"]`);

    const W = 720, H = 330;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
                            "aria-label": "Two CDSs, their reads, and the EM estimate over iterations" }, root);
    const readout = document.createElement("p");
    readout.className = "rk-readout";
    root.appendChild(readout);

    let m, res, k = 0, timer = null;

    function recompute() {
      const rho = [+q("rA").value, +q("rB").value];
      q("rA").nextElementSibling.textContent = rho[0].toFixed(1);
      q("rB").nextElementSibling.textContent = rho[1].toFixed(1);
      m = model(SCENARIOS[q("sc").value], rho);
      res = em(m, q("len").checked);
      root.rk = { model: m, iterations: res.trace.length - 1, alpha: res.trace[res.trace.length - 1] };
    }

    function draw() {
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const sc = SCENARIOS[q("sc").value];
      const K = res.trace.length - 1, alpha = res.trace[k];

      // Structure: the CDSs as rows, regions coloured by which CDSs contain them.
      const total = sc.regions.reduce((s, r) => s + r.len, 0), x0 = 40, sx = 640 / total;
      let x = x0;
      NAMES.forEach((n, t) => label(svg, 14, 31 + 30 * t, n, { "font-weight": "bold", fill: C[n] }));
      for (const r of sc.regions) {
        const key = r.members.map((t) => NAMES[t]).join("");
        const w = r.len * sx;
        for (const t of r.members) el("rect", { x, y: 18 + 30 * t, width: w - 2, height: 18, rx: 3, fill: C[key] }, svg);
        const n = r.len * r.members.reduce((s, t) => s + [+q("rA").value, +q("rB").value][t], 0);
        label(svg, x + w / 2, 31 + 30 * (r.members.length === 1 ? r.members[0] : 0), `${r.len} nt`,
              { fill: "white", "font-size": 11 });
        label(svg, x + w / 2, 92, `${fmt(n)} reads → {${r.members.map((t) => NAMES[t]).join(", ")}}`,
              { fill: C[key], "font-size": 12 });
        x += w;
      }

      // Bars: truth (outline) vs current estimate (filled).
      const top = 128, ph = 150, base = top + ph;
      const ymax = res.trace.reduce((mx, a) => Math.max(mx, a[0], a[1]), Math.max(1, ...m.truth)) * 1.1;
      const Y = (v) => base - (v / ymax) * ph;
      el("line", { x1: 0, x2: 250, y1: base, y2: base, stroke: "#b0bec5" }, svg);
      NAMES.forEach((n, t) => {
        const bx = 30 + t * 115;
        el("rect", { x: bx, y: Y(alpha[t]), width: 40, height: base - Y(alpha[t]), fill: C[n], rx: 2 }, svg);
        el("rect", { x: bx + 44, y: Y(m.truth[t]), width: 40, height: base - Y(m.truth[t]), fill: "none",
                     stroke: C[n], "stroke-dasharray": "4,3", "stroke-width": 1.5, rx: 2 }, svg);
        label(svg, bx + 20, Y(alpha[t]) - 5, fmt(alpha[t]), { "font-size": 11 });
        label(svg, bx + 64, Y(m.truth[t]) - 5, fmt(m.truth[t]), { "font-size": 11, fill: C.muted });
        label(svg, bx + 20, base + 15, `${n} EM`, { "font-size": 11 });
        label(svg, bx + 64, base + 15, `${n} true`, { "font-size": 11, fill: C.muted });
      });
      label(svg, 125, top - 10, "expected reads", { fill: C.muted });

      // Trace: estimate per iteration, truth dashed.
      const px0 = 300, pw = W - px0 - 10, X = (i) => px0 + (i / Math.max(K, 10)) * pw;
      el("line", { x1: px0, x2: px0 + pw, y1: base, y2: base, stroke: "#b0bec5" }, svg);
      el("line", { x1: px0, x2: px0, y1: top, y2: base, stroke: "#b0bec5" }, svg);
      NAMES.forEach((n, t) => {
        el("line", { x1: px0, x2: px0 + pw, y1: Y(m.truth[t]), y2: Y(m.truth[t]), stroke: C[n],
                     "stroke-dasharray": "4,3", opacity: 0.6 }, svg);
        const pts = res.trace.slice(0, k + 1).map((a, i) => `${X(i).toFixed(1)},${Y(a[t]).toFixed(1)}`).join(" ");
        el("polyline", { points: pts, fill: "none", stroke: C[n], "stroke-width": 2 }, svg);
        el("circle", { cx: X(k), cy: Y(alpha[t]), r: 3.5, fill: C[n] }, svg);
      });
      label(svg, px0 + pw / 2, top - 10, "estimate per iteration (dashed: truth)", { fill: C.muted });
      label(svg, px0, base + 15, "0", { "font-size": 11, fill: C.muted });
      label(svg, px0 + pw, base + 15, String(Math.max(K, 10)), { "font-size": 11, fill: C.muted, "text-anchor": "end" });
      label(svg, px0 + pw / 2, base + 30, "iteration", { fill: C.muted, "font-size": 11 });

      // Readout.
      const done = k === K;
      let s = `<b>Iteration ${k}</b>` + (k === 0 ? " (even start: each shared read split equally)" :
        `, largest change ${res.changes[k] < 1e-6 ? "0" : res.changes[k].toPrecision(2)} reads`) +
        (done ? ` · <b>converged</b> after ${K} iterations (no estimate moves by ${TOL} reads or more).` : ".");
      if (done) {
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
  }

  if (typeof document$ !== "undefined") document$.subscribe(init);
  else document.addEventListener("DOMContentLoaded", init);
})();
