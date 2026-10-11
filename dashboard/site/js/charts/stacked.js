// The additive decomposition of a disease's prediction: one bar per protein family, one segment per
// contributing phenolog (the paper's case study figure, interactive). rows = additive.parquet rows of
// one disease. Everything is re-sliced client side from the per (family, phenolog) rows.
import { esc, fmt, speciesOf, route } from "../format.js";
import { showTip, moveTip, hideTip } from "./tooltip.js";

const PALETTE = [...d3.schemeTableau10, ...d3.schemeSet3.filter(c => c !== "#d9d9d9" && c !== "#ffffb3")];

export const ORDERS = {
  total: ["Summed contribution", f => -f.total],
  count: ["Number of phenologs", f => -f.segments.length],
  max: ["Strongest single phenolog", f => -Math.max(...f.segments.map(s => s.score))],
  known: ["Families with known genes first", f => [f.known ? 0 : 1, -f.total]],
  name: ["Family id", f => f.family],
};

export function slice(rows, { k = 20, species = null, topFamilies = 20, order = "total" } = {}) {
  const active = rows.filter(r => r.phenolog_rank <= k && (!species || species.has(r.species)));
  const fams = new Map();
  for (const r of active) {
    if (!fams.has(r.family)) fams.set(r.family, { family: r.family, segments: [], total: 0, known: r.known_symbols, other: r.other_symbols });
    const f = fams.get(r.family);
    f.segments.push(r);
    f.total += r.score;
  }
  const cmp = ORDERS[order][1];
  const list = [...fams.values()];
  list.forEach(f => f.segments.sort((a, b) => a.phenolog_rank - b.phenolog_rank));
  list.sort((a, b) => {
    const x = cmp(a), y = cmp(b);
    return Array.isArray(x) ? (x[0] - y[0]) || (x[1] - y[1]) : (x < y ? -1 : x > y ? 1 : 0);
  });
  const phenologs = [...new Map(active.map(r => [r.phenolog_rank, r])).values()].sort((a, b) => a.phenolog_rank - b.phenolog_rank);
  return { families: list.slice(0, topFamilies), phenologs };
}

export function additiveChart(el, rows, opts = {}) {
  const { families, phenologs } = slice(rows, opts);
  el.innerHTML = "";
  if (!families.length) { el.innerHTML = `<div class="empty">No contributing protein families under these settings</div>`; return; }

  const color = d3.scaleOrdinal().domain(phenologs.map(p => p.phenolog_rank)).range(PALETTE);
  const width = Math.max(680, el.clientWidth || 900);
  const rowH = 22, mLeft = 112, mRight = 300, mTop = 8, mBottom = 30;
  const height = mTop + families.length * rowH + mBottom;
  const x = d3.scaleLinear().domain([0, d3.max(families, f => f.total) || 1]).nice().range([0, width - mLeft - mRight]);

  const svg = d3.select(el).append("svg").attr("width", width).attr("height", height).attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${mLeft},${mTop})`);

  g.append("g").attr("class", "axis").attr("transform", `translate(0,${families.length * rowH + 4})`)
    .call(d3.axisBottom(x).ticks(6)).call(a => a.select(".domain").remove());
  g.append("g").selectAll("line").data(x.ticks(6)).join("line").attr("class", "gridline")
    .attr("x1", d => x(d)).attr("x2", d => x(d)).attr("y1", 0).attr("y2", families.length * rowH);
  svg.append("text").attr("x", mLeft + x.range()[1] / 2).attr("y", height - 2).attr("text-anchor", "middle")
    .attr("font-size", 11).attr("fill", "#8b97a8").text("additive prediction score  Σ (overlap / disease orthogroups) × (1 − p)");

  const fam = g.selectAll("g.fam").data(families).join("g").attr("class", "fam")
    .attr("transform", (f, i) => `translate(0,${i * rowH})`);

  fam.append("a").attr("href", f => route.family(f.family)).append("text")
    .attr("x", -8).attr("y", rowH / 2 + 4).attr("text-anchor", "end").attr("font-size", 11.5).attr("font-family", "JetBrains Mono, monospace")
    .attr("fill", "#4fc3f7").text(f => f.family);

  fam.each(function (f) {
    let acc = 0;
    d3.select(this).selectAll("rect").data(f.segments.map(s => { const d = { ...s, x0: acc }; acc += s.score; return d; })).join("rect")
      .attr("x", d => x(d.x0)).attr("y", 3).attr("height", rowH - 6)
      .attr("width", d => Math.max(0, x(d.x0 + d.score) - x(d.x0)))
      .attr("fill", d => color(d.phenolog_rank)).attr("stroke", "#0d1117").attr("stroke-width", 0.6)
      .on("mousemove", (e, d) => {
        showTip(e, `<b>${esc(d.phenotype_name)}</b><br><i>${esc(speciesOf(d.species).binomial)}</i> · phenolog #${d.phenolog_rank}<br>` +
                   `contribution ${fmt.score(d.score)} of ${fmt.score(f.total)}`);
        moveTip(e);
      })
      .on("mouseleave", hideTip);
  });

  // gene labels: known genes first, in the accent colour
  fam.append("text").attr("x", f => x(f.total) + 8).attr("y", rowH / 2 + 4).attr("font-size", 11.5).each(function (f) {
    const t = d3.select(this);
    const known = (f.known || "").split(", ").filter(Boolean);
    const other = (f.other || "").split(", ").filter(Boolean);
    const shown = [...known.map(s => [s, true]), ...other.map(s => [s, false])];
    const limit = 4;
    shown.slice(0, limit).forEach(([s, isKnown], i) => {
      t.append("tspan").attr("fill", isKnown ? "#f0b84a" : "#e6ebf2").attr("font-weight", isKnown ? 700 : 400)
        .text((i ? ", " : "") + s);
    });
    if (shown.length > limit) t.append("tspan").attr("fill", "#5c687a").text(`  +${shown.length - limit}`);
  });

  // legend: phenologs in rank order
  const legend = document.createElement("div");
  legend.className = "legend";
  legend.innerHTML = phenologs.map(p => `
    <span class="item"><span class="sw" style="background:${color(p.phenolog_rank)}"></span>
      <span class="faint">#${p.phenolog_rank}</span>
      <a href="${route.phenotype(p.phenotype_id)}">${esc(p.phenotype_name)}</a>
      <i class="muted">${esc(speciesOf(p.species).binomial)}</i></span>`).join("");
  el.appendChild(legend);
}
