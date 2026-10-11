// A radar (player card) chart. axes: [{ label, value (0..1), text }]
import { showTip, moveTip, hideTip } from "./tooltip.js";

export function radar(el, axes, { size = 150, color = "#f0b84a", levels = 4, labels = false } = {}) {
  const pad = labels ? 64 : 6;
  const r = (size - 2 * pad) / 2;
  const cx = size / 2, cy = size / 2;
  const n = axes.length;
  const angle = i => -Math.PI / 2 + (2 * Math.PI * i) / n;
  const pt = (i, v) => [cx + r * v * Math.cos(angle(i)), cy + r * v * Math.sin(angle(i))];

  const svg = d3.select(el).append("svg").attr("width", size).attr("height", size).attr("viewBox", `0 0 ${size} ${size}`).style("overflow", "visible");

  for (let l = 1; l <= levels; l++) {
    const v = l / levels;
    svg.append("polygon")
      .attr("points", axes.map((_, i) => pt(i, v).join(",")).join(" "))
      .attr("fill", l === levels ? "rgba(255,255,255,0.02)" : "none")
      .attr("stroke", "#273243").attr("stroke-width", 1);
  }
  axes.forEach((_, i) => {
    const [x, y] = pt(i, 1);
    svg.append("line").attr("x1", cx).attr("y1", cy).attr("x2", x).attr("y2", y).attr("stroke", "#273243");
  });

  const poly = axes.map((a, i) => pt(i, Math.max(0, Math.min(1, a.value || 0))));
  svg.append("polygon").attr("points", poly.map(p => p.join(",")).join(" "))
    .attr("fill", color).attr("fill-opacity", 0.28).attr("stroke", color).attr("stroke-width", 2).attr("stroke-linejoin", "round");

  svg.selectAll("circle.v").data(axes).join("circle").attr("class", "v")
    .attr("cx", (a, i) => poly[i][0]).attr("cy", (a, i) => poly[i][1]).attr("r", 3)
    .attr("fill", "#fff").attr("stroke", color).attr("stroke-width", 1.5)
    .on("mousemove", (e, a) => { showTip(e, `<b>${a.label}</b><br>${a.text ?? ""}`); moveTip(e); })
    .on("mouseleave", hideTip);

  if (labels) {
    axes.forEach((a, i) => {
      const [x, y] = pt(i, 1.17);
      svg.append("text").attr("x", x).attr("y", y).attr("text-anchor", "middle").attr("dominant-baseline", "middle")
        .attr("font-size", 10).attr("fill", "#8b97a8").text(a.label);
      svg.append("text").attr("x", x).attr("y", y + 12).attr("text-anchor", "middle").attr("dominant-baseline", "middle")
        .attr("font-size", 10.5).attr("font-weight", 600).attr("fill", color).text(a.text ?? "");
    });
  }
  return svg;
}

// Normalise one species' stats against the field (log scale, so a species with 83 phenologs still
// shows beside one with 50,000)
export function speciesAxes(s, all, defs) {
  return defs.map(([key, label, format]) => {
    const max = Math.max(1, ...all.map(x => x[key] || 0));
    const v = s[key] || 0;
    return { label, value: Math.log10(1 + v) / Math.log10(1 + max), text: format ? format(v) : v.toLocaleString() };
  });
}
