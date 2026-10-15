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

// One bar per species, in the order given and in species colours, with the number beside it.
// items: [{species, label, value, color, href, text}]. log scales bars by log10(1 + value), so the
// small species stay visible next to the large ones
export function speciesBars(el, items, { log = false } = {}) {
  const f = log ? v => Math.log10(1 + v) : v => v;
  const top = Math.max(...items.map(d => f(d.value)), Number.MIN_VALUE);
  el.innerHTML = `<div class="sp-bars">` + items.map(d => `
    <a class="row" data-sp="${esc(d.species)}" href="${d.href}">
      <span class="lab"><i>${esc(d.label)}</i></span>
      <span class="track"><span class="bar" style="width:${d.value > 0 ? Math.max(1, 100 * f(d.value) / top).toFixed(1) : 0}%; background:${d.color}"></span></span>
      <span class="n">${d.text ?? fmt.int(d.value)}</span>
    </a>`).join("") + `</div>`;
}

// A funnel: each step keeps part of the one before. steps: [{label, n}]; bars and percentages are of the first step
export function funnel(el, steps, { color = "var(--accent)" } = {}) {
  const first = Math.max(1, steps[0]?.n || 0);
  el.innerHTML = `<div class="funnel">` + steps.map(st => `
    <div class="step">
      <span class="lab">${st.label}</span>
      <span class="track"><span class="bar" style="width:${(100 * st.n / first).toFixed(1)}%; background:${color}"></span></span>
      <span class="n">${fmt.int(st.n)}</span>
      <span class="pct">${fmt.pct(st.n / first)}</span>
    </div>`).join("") + `</div>`;
}
