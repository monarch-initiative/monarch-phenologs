// The species tree of a release: a dated cladogram (divergence from the base species, Mya) with a
// metric drawn beside each leaf. tree = nested {name, taxon_id, time, children}; species = rows of species.parquet
import { esc, fmt, speciesOf, route } from "./format.js";
import { showTip, moveTip, hideTip } from "./tooltip.js";

export function speciesTree(el, tree, species, { metric = "n_phenologs", metricLabel = "phenologs", width = 760 } = {}) {
  el.innerHTML = "";
  if (!tree) { el.innerHTML = `<div class="empty">No tree for this release</div>`; return; }
  const byName = new Map(species.map(s => [s.species, s]));
  const root = d3.hierarchy(tree);
  const leaves = root.leaves();
  const rowH = 34, labelW = 250, barW = 170, mLeft = 40, mTop = 12, mBottom = 38;
  const treeW = width - labelW - barW - mLeft - 20;
  const height = leaves.length * rowH + mTop + mBottom;

  const maxTime = d3.max(root.descendants(), d => d.data.time) || 1;
  const x = d3.scaleSqrt().domain([0, maxTime]).range([treeW, 0]);
  leaves.forEach((d, i) => { d.y = mTop + i * rowH + rowH / 2; });
  root.eachAfter(d => { if (d.children) d.y = d3.mean(d.children, c => c.y); });
  root.each(d => { d.x = x(d.data.time); });

  const svg = d3.select(el).append("svg").attr("width", width).attr("height", height).attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${mLeft},0)`);

  // time axis
  const axis = d3.axisBottom(x).tickValues([0, 100, 300, 600, 1000].filter(t => t <= maxTime)).tickFormat(d => d === 0 ? "now" : d);
  g.append("g").attr("class", "axis").attr("transform", `translate(0,${height - mBottom + 6})`).call(axis)
    .call(a => a.select(".domain").remove());
  g.append("text").attr("x", 0).attr("y", height - 3).attr("text-anchor", "start").attr("font-size", 10.5).attr("fill", "#5c687a")
    .text("divergence from " + (speciesOf(species.find(s => s.is_base)?.species || "").binomial) + ", million years ago (TimeTree, approximate)");

  // elbow links
  g.selectAll("path.link").data(root.links()).join("path").attr("class", "link")
    .attr("d", d => `M${d.source.x},${d.source.y}V${d.target.y}H${d.target.x}`)
    .attr("fill", "none").attr("stroke", "#3a4759").attr("stroke-width", 1.5);

  // inner node times
  g.selectAll("text.t").data(root.descendants().filter(d => d.children)).join("text").attr("class", "t")
    .attr("x", d => d.x + 4).attr("y", d => d.y - 5).attr("text-anchor", "start").attr("font-size", 9.5).attr("fill", "#5c687a")
    .text(d => d.data.time);

  // leaves
  const maxMetric = Math.max(1, ...species.filter(s => !s.is_base).map(s => s[metric] || 0));
  const bar = d3.scaleLinear().domain([0, maxMetric]).range([0, barW - 60]);
  const leaf = g.selectAll("g.leaf").data(leaves).join("g").attr("class", "leaf")
    .attr("transform", d => `translate(${d.x},${d.y})`).style("cursor", "pointer")
    .on("click", (e, d) => { location.hash = route.species(d.data.name); });
  leaf.append("circle").attr("r", 6).attr("fill", d => speciesOf(d.data.name).color).attr("stroke", "#0d1117").attr("stroke-width", 1.5);
  leaf.append("text").attr("x", 12).attr("y", -2).attr("font-size", 12.5).attr("font-style", "italic").attr("fill", "#e6ebf2")
    .text(d => speciesOf(d.data.name).binomial);
  leaf.append("text").attr("x", 12).attr("y", 11).attr("font-size", 10.5).attr("fill", "#8b97a8")
    .text(d => { const s = byName.get(d.data.name); return s ? (s.is_base ? "base species · " + s.ontology : `${s.common}${s.ontology ? " · " + s.ontology : ""}`) : ""; });

  const barX = treeW + labelW;
  leaf.filter(d => !byName.get(d.data.name)?.is_base).each(function (d) {
    const s = byName.get(d.data.name) || {};
    const v = s[metric] || 0;
    const gg = d3.select(this).append("g").attr("transform", `translate(${barX - d.x},0)`);
    gg.append("rect").attr("x", 0).attr("y", -7).attr("height", 14).attr("width", Math.max(1.5, bar(v))).attr("rx", 3)
      .attr("fill", speciesOf(d.data.name).color).attr("fill-opacity", 0.85)
      .on("mousemove", e => { showTip(e, `<b>${esc(speciesOf(d.data.name).binomial)}</b><br>${fmt.int(v)} ${esc(metricLabel)}`); moveTip(e); })
      .on("mouseleave", hideTip);
    gg.append("text").attr("x", Math.max(1.5, bar(v)) + 6).attr("y", 4).attr("font-size", 11.5).attr("fill", "#e6ebf2").text(fmt.int(v));
  });
  leaf.filter(d => byName.get(d.data.name)?.is_base).append("text").attr("x", d => barX - d.x).attr("y", 4)
    .attr("font-size", 11).attr("fill", "#5c687a").text("—");
  g.append("text").attr("x", barX).attr("y", mTop - 1).attr("font-size", 10.5).attr("fill", "#8b97a8").text(metricLabel);
}
