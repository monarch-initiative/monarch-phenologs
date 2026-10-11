// A ranked list with proportional bars. items: [{label, value, href, color, sub}]
import { esc, fmt } from "../format.js";

export function barList(el, items, { max = null, color = null } = {}) {
  const top = max ?? Math.max(1, ...items.map(d => d.value));
  el.innerHTML = `<div class="bar-list">` + items.map(d => `
    <div class="row">
      <span class="lab" title="${esc(d.label)}">${d.href ? `<a href="${d.href}">${esc(d.label)}</a>` : esc(d.label)}
        ${d.sub ? `<span class="faint"> ${d.sub}</span>` : ""}</span>
      <span><span class="bar" style="width:${(100 * d.value / top).toFixed(1)}%; background:${d.color || color || "var(--accent)"}"></span></span>
      <span class="n">${fmt.int(d.value)}</span>
    </div>`).join("") + `</div>`;
}
