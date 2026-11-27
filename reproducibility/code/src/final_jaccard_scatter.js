
/* Offline permutation scatter for saved module-pair statistics. */
(() => {
  const SVG_NS = "http://www.w3.org/2000/svg";
  const bounds = { left: 76, right: 875, top: 108, bottom: 520 };
  const width = 900;
  const height = 580;
  const maxLogP = 4.2;
  const maxJaccard = .25;
  const perPage = 50;
  const highlightCount = 25;
  const outlierLabel = dataset => dataset.afterOutlier ? "After outlier removal" : "Before outlier removal";
  const xPosition = value => bounds.left + value / maxJaccard * (bounds.right - bounds.left);
  const yPosition = value => bounds.bottom - value / maxLogP * (bounds.bottom - bounds.top);
  const pText = value => value < .001 ? value.toExponential(2) : String(Number(value.toPrecision(3)));
  const cellLabel = value => value === "MAIT_gDT_cell" ? "MAIT / γδ T cell" : value.replaceAll("_", " ");
  const esc = value => String(value).replace(/[&<>"']/g, char => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[char]
  ));

  function svgNode(tag, attrs, value) {
    const node = document.createElementNS(SVG_NS, tag);
    for (const [key, item] of Object.entries(attrs || {})) node.setAttribute(key, item);
    if (value !== undefined) node.textContent = value;
    return node;
  }

  const svgStyle = `
    text{font:12px Arial,Helvetica,sans-serif;fill:#c0c6d0}
    .scatter-title{font-weight:bold;font-size:18px;fill:#f4f6f9}
    .scatter-subtitle{font-size:12px}
    .scatter-axis-title{font-size:14px;fill:#f4f6f9}
    .scatter-grid{stroke:#363b42;stroke-width:1}
    .scatter-axis{stroke:#79808a;stroke-width:1}
    .scatter-point{fill:#9ba1ab;fill-opacity:.28;cursor:pointer}
    .scatter-point.included{fill:#58aae8;fill-opacity:.9}
    .scatter-point:hover{fill:#eef7ff;fill-opacity:1}
    .scatter-cut{stroke:#a6b7c8;stroke-dasharray:5 5;stroke-width:1.2}
    .scatter-cut-text{fill:#e5effa;paint-order:stroke;stroke:#1b1d20;stroke-width:5}
    .scatter-cut-hit{stroke:transparent;stroke-width:14;cursor:ns-resize;touch-action:none}
    .scatter-rank{fill:#deefff;font-size:12px;font-weight:bold;paint-order:stroke;stroke:#1b1d20;stroke-width:4}
    .scatter-table-head{font-weight:bold}
  `;

  function createPlot(points, dataset, href) {
    const svg = svgNode("svg", {
      viewBox: `0 0 ${width} ${height}`,
      "aria-label": "kME-weighted Jaccard versus permutation p scatter plot",
    });
    svg.append(svgNode("style", {}, svgStyle),
      svgNode("rect", { width, height, fill: "#1b1d20" }),
      svgNode("text", { x: 23, y: 26, class: "scatter-title" },
        `${dataset.label} · ${dataset.cells} cells · HVG 9,000`),
      svgNode("text", { x: 23, y: 46, class: "scatter-subtitle" },
        `CSF RNA × CSF protein · ${dataset.afterOutlier ? "After" : "Before"} common outlier removal · Other excluded`),
      svgNode("circle", { cx: 27, cy: 70, r: 4, fill: "#9ba1ab", opacity: .6 }),
      svgNode("text", { x: 39, y: 74 }, "Above p cutoff"),
      svgNode("circle", { cx: 165, cy: 70, r: 4, fill: "#58aae8" }),
      svgNode("text", { x: 177, y: 74 }, "At or below cutoff"),
      svgNode("circle", { cx: 321, cy: 70, r: 5.4, fill: "none", stroke: "#daeeff" }),
      svgNode("text", { x: 335, y: 74 }, "Top 25 by Jaccard")
    );

    for (let tick = 0; tick <= 5; tick++) {
      const number = tick * .05, x = xPosition(number);
      svg.append(svgNode("line", { x1: x, x2: x, y1: bounds.top, y2: bounds.bottom, class: "scatter-grid" }),
        svgNode("text", { x, y: bounds.bottom + 19, "text-anchor": "middle" }, number.toFixed(2)));
    }
    for (let tick = 0; tick <= 4; tick++) {
      const y = yPosition(tick);
      svg.append(svgNode("line", { x1: bounds.left, x2: bounds.right, y1: y, y2: y, class: "scatter-grid" }),
        svgNode("text", { x: bounds.left - 13, y: y + 4, "text-anchor": "end" }, String(tick)));
    }
    svg.append(
      svgNode("path", { d: `M${bounds.left},${bounds.top}V${bounds.bottom}H${bounds.right}`, fill: "none", class: "scatter-axis" }),
      svgNode("text", { x: (bounds.left + bounds.right) / 2, y: 565, "text-anchor": "middle", class: "scatter-axis-title" },
        "kME-weighted Jaccard similarity"),
      svgNode("text", { transform: "translate(25 314) rotate(-90)", "text-anchor": "middle", class: "scatter-axis-title" },
        "−log₁₀(permutation p)")
    );

    const cutLine = svgNode("line", { x1: bounds.left, x2: bounds.right, class: "scatter-cut" });
    const drag = svgNode("line", { x1: bounds.left, x2: bounds.right, class: "scatter-cut-hit", "data-ui": 1 });
    const cutoffLabel = svgNode("text", { x: bounds.right - 7, "text-anchor": "end", class: "scatter-cut-text", "pointer-events": "none" });
    const unselected = svgNode("g"), selected = svgNode("g");
    const ranks = svgNode("g", { "pointer-events": "none", "aria-hidden": "true" });
    const plotted = points.map(pair => {
      const x = xPosition(pair.jaccard), y = yPosition(-Math.log10(pair.p));
      const link = svgNode("a", { href: href(pair), target: "_top", tabindex: "-1" });
      const circle = svgNode("circle", { cx: x, cy: y, r: 2.5, class: "scatter-point" });
      link.append(svgNode("title", {},
        `${cellLabel(pair.cell)} · RNA ${pair.transcript_module} × protein ${pair.proteomics_module}\nJaccard ${pair.jaccard.toFixed(4)} · p ${pText(pair.p)} · q ${pText(pair.q)} · ${pair.shared} shared genes`), circle);
      unselected.append(link);
      return { ...pair, x, y, link, circle };
    });
    svg.append(cutLine, drag, unselected, selected, ranks, cutoffLabel);
    return { svg, cutLine, drag, cutoffLabel, unselected, selected, plotted, ranks, active: new Set() };
  }

  function drawRanks(layer, ranked) {
    layer.replaceChildren();
    const occupied = [];
    for (const [index, point] of ranked.slice(0, highlightCount).entries()) {
      let label;
      for (const radius of [14, 20, 27, 36, 47, 63, 81]) {
        for (const degrees of [-90, -45, -135, 0, 180, 45, 135, 90]) {
          const radians = degrees * Math.PI / 180;
          const x = point.x + radius * Math.cos(radians), y = point.y + radius * Math.sin(radians);
          if (x < bounds.left + 9 || x > bounds.right - 9 || y < bounds.top + 8 || y > bounds.bottom - 8) continue;
          if (occupied.some(([ox, oy]) => Math.abs(ox - x) < 22 && Math.abs(oy - y) < 16)) continue;
          if (ranked.slice(0, highlightCount).some(p => Math.hypot(p.x - x, p.y - y) < 11)) continue;
          label = [x, y];
          break;
        }
        if (label) break;
      }
      layer.append(svgNode("circle", { cx: point.x, cy: point.y, r: 5.2, fill: "none", stroke: "#daeeff" }));
      if (!label) continue;
      const [x, y] = label;
      occupied.push(label);
      layer.append(svgNode("line", { x1: point.x, x2: x, y1: point.y, y2: y, stroke: "#9ec3de", "stroke-width": .7 }),
        svgNode("text", { x, y: y + 4, "text-anchor": "middle", class: "scatter-rank" }, String(index + 1)));
    }
  }

  function exportPlot(plot, eligible, cutoff, dataset) {
    const copy = plot.svg.cloneNode(true);
    const exportedWidth = width + 438;
    copy.setAttribute("xmlns", SVG_NS);
    copy.setAttribute("viewBox", `0 0 ${exportedWidth} ${height}`);
    copy.setAttribute("width", exportedWidth);
    copy.setAttribute("height", height);
    copy.querySelectorAll("[data-ui]").forEach(node => node.remove());
    copy.querySelectorAll("a").forEach(link => link.replaceWith(...link.childNodes));
    copy.querySelector("rect").setAttribute("width", exportedWidth);
    const start = width + 11;
    const top = eligible.slice(0, highlightCount);
    copy.append(svgNode("line", { x1: width, x2: width, y1: 9, y2: height - 12, class: "scatter-grid" }),
      svgNode("text", { x: start, y: 26, class: "scatter-title" }, "Ranked module pairs"),
      svgNode("text", { x: start, y: 47 },
        `${top.length} shown of ${eligible.length.toLocaleString()} pairs with p ≤ ${pText(cutoff)}`));
    const columns = [[0,"#"], [29,"RNA cell type"], [189,"RNA module"], [287,"Protein module"]];
    for (const [x, name] of columns) copy.append(svgNode("text", { x: start + x, y: 74, class: "scatter-table-head" }, name));
    for (const [index, pair] of top.entries()) {
      const y = 95 + index * 18;
      if (index % 2 === 0) copy.append(svgNode("rect", { x: start - 4, y: y - 13, width: 405, height: 18, fill: "#242a30" }));
      const values = [String(index + 1), cellLabel(pair.cell), pair.transcript_module, pair.proteomics_module];
      values.forEach((text, col) => copy.append(svgNode("text", { x: start + columns[col][0], y }, text)));
    }
    const blob = new Blob([new XMLSerializer().serializeToString(copy)], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob), link = document.createElement("a");
    link.href = url;
    const branch = dataset.afterOutlier ? "after_outlier" : "before_outlier";
    link.download = `${dataset.stem}_hvg9000_${branch}_permutation_scatter.svg`;
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function create(root, datasets, { href, onCorrection }) {
    root.innerHTML = `
      <div class="scatter-controls">
        <div><span class="control-label">CSF analysis condition</span><div class="button-row" data-conditions></div></div>
        <div class="actions scatter-downloads"><button type="button" data-svg>Save SVG</button>
          <a class="action" download data-csv>Download source CSV</a></div>
      </div>
      <div class="scatter-threshold">
        <label for="scatterPThreshold">Permutation p cutoff <output id="scatterPValue" for="scatterPThreshold"></output></label>
        <input id="scatterPThreshold" type="range" min="0" max="4.10" step=".01" value="2">
        <button type="button" data-reset>Reset to 0.01</button>
        <span id="scatterCount" role="status"></span>
      </div>
      <div class="scatter-layout">
        <div class="scatter-plot" id="scatterPlot"></div>
        <div class="scatter-pairs">
          <h2>Pairs at or below the cutoff</h2>
          <div class="table-wrap scatter-table"><table class="ranking">
            <thead><tr><th>#</th><th>RNA cell type</th><th>RNA module</th><th>Protein module</th><th>Jaccard</th><th>p</th><th>BH-FDR q</th><th>Genes</th></tr></thead>
            <tbody id="scatterPairs"></tbody></table></div>
          <div class="scatter-pages"><button type="button" data-prev>Previous</button><span data-page></span><button type="button" data-next>Next</button></div>
        </div>
      </div>`;
    const get = selector => root.querySelector(selector);
    const slider = get("#scatterPThreshold"), pValue = get("#scatterPValue");
    const plotContainer = get("#scatterPlot"), conditions = get("[data-conditions]");
    const prev = get("[data-prev]"), next = get("[data-next]");
    let dataset, plot, conditionKey, cutoff, eligible = [], page = 0, frame = null;
    const entries = Object.entries(datasets);
    const previousCondition = new Map();
    conditions.innerHTML = entries.map(([key, data]) =>
      `<button type="button" data-condition="${esc(key)}">${esc(data.label)} · ${data.cells} cells · ${outlierLabel(data)}</button>`).join("");
    const buttons = [...conditions.querySelectorAll("button")];
    buttons.forEach(button => {
      button.onclick = () => {
        setCondition(button.dataset.condition);
        onCorrection(dataset.correction);
      };
    });

    function renderPairs() {
      const start = page * perPage, end = Math.min(start + perPage, eligible.length);
      get("#scatterPairs").innerHTML = eligible.slice(start, end).map((pair, index) => `
        <tr><td>${start + index + 1}</td><td>${esc(cellLabel(pair.cell))}</td><td>${esc(pair.transcript_module)}</td>
        <td>${esc(pair.proteomics_module)}</td><td class="num">${pair.jaccard.toFixed(4)}</td>
        <td class="num">${pText(pair.p)}</td><td class="num">${pText(pair.q)}</td>
        <td><a class="action" target="_top" href="${esc(href(pair))}">Inspect</a></td></tr>`).join("")
        || '<tr><td colspan="8" class="empty">No pairs meet this cutoff.</td></tr>';
      prev.disabled = page === 0;
      next.disabled = end >= eligible.length;
      get("[data-page]").textContent = eligible.length ? `${start + 1}–${end} of ${eligible.length.toLocaleString()}` : "0 pairs";
      get(".scatter-table").scrollTop = 0;
    }

    function applyCutoff() {
      if (frame !== null) cancelAnimationFrame(frame);
      frame = null;
      if (!plot) return;
      const logP = Number(slider.value), y = yPosition(logP);
      cutoff = 10 ** -logP;
      eligible = plot.plotted.filter(pair => pair.p <= cutoff);
      const now = new Set(eligible);
      for (const point of plot.active) {
        if (now.has(point)) continue;
        point.circle.classList.remove("included");
        point.circle.setAttribute("r", 2.5);
        plot.unselected.append(point.link);
      }
      for (const point of now) {
        if (plot.active.has(point)) continue;
        point.circle.classList.add("included");
        point.circle.setAttribute("r", 3.4);
        plot.selected.append(point.link);
      }
      plot.active = now;
      for (const line of [plot.cutLine, plot.drag]) {
        line.setAttribute("y1", y);
        line.setAttribute("y2", y);
      }
      plot.cutoffLabel.setAttribute("y", Math.max(bounds.top + 13, y - 8));
      plot.cutoffLabel.textContent = `p = ${pText(cutoff)}`;
      pValue.textContent = pText(cutoff);
      slider.setAttribute("aria-valuetext", `p ${pText(cutoff)}`);
      get("#scatterCount").textContent = `${eligible.length.toLocaleString()} of ${plot.plotted.length.toLocaleString()} module pairs`;
      drawRanks(plot.ranks, eligible);
      page = 0;
      renderPairs();
    }

    function scheduleCutoff() {
      if (frame === null) frame = requestAnimationFrame(applyCutoff);
    }

    function setCondition(key) {
      if (conditionKey === key) return;
      dataset = datasets[key];
      conditionKey = key;
      previousCondition.set(dataset.correction, key);
      buttons.forEach(button => {
        const selected = button.dataset.condition === key;
        button.classList.toggle("active", selected);
        button.setAttribute("aria-pressed", selected);
      });
      const points = dataset.points.map(([record, transcript_module, proteomics_module, jaccard, p, q, shared]) => {
        const [recordId, cell] = dataset.records[record];
        return { recordId, cell, transcript_module, proteomics_module, jaccard, p, q, shared };
      }).sort((a, b) => b.jaccard - a.jaccard || a.p - b.p || a.cell.localeCompare(b.cell) ||
        a.transcript_module.localeCompare(b.transcript_module) || a.proteomics_module.localeCompare(b.proteomics_module));
      plot = createPlot(points, dataset, href);
      plotContainer.replaceChildren(plot.svg);
      plot.drag.onpointerdown = event => {
        event.preventDefault();
        plot.drag.setPointerCapture(event.pointerId);
      };
      plot.drag.onpointermove = event => {
        if (!plot.drag.hasPointerCapture(event.pointerId)) return;
        const cursor = new DOMPoint(event.clientX, event.clientY).matrixTransform(plot.svg.getScreenCTM().inverse());
        const logP = (bounds.bottom - cursor.y) / (bounds.bottom - bounds.top) * maxLogP;
        slider.value = Math.max(0, Math.min(4.1, logP));
        scheduleCutoff();
      };
      get("[data-csv]").href = dataset.csv;
      applyCutoff();
    }

    function setCorrection(correction) {
      const key = previousCondition.get(correction) || entries.find(([, data]) => data.correction === correction)?.[0];
      if (key) {
        setCondition(key);
      } else if (!conditionKey && entries.length) {
        setCondition(entries[0][0]);
        onCorrection(dataset.correction);
      }
    }

    get("#scatterPairs").onclick = event => {
      if (!event.target.closest("a")) event.target.closest("tr")?.querySelector("a")?.click();
    };
    slider.oninput = scheduleCutoff;
    get("[data-reset]").onclick = () => { slider.value = 2; applyCutoff(); };
    get("[data-svg]").onclick = () => {
      if (frame !== null) applyCutoff();
      exportPlot(plot, eligible, cutoff, dataset);
    };
    prev.onclick = () => { page--; renderPairs(); };
    next.onclick = () => { page++; renderPairs(); };

    return { setCorrection };
  }

  window.FinalJaccardScatter = { create };
})();
