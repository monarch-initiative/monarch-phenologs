// Small HTML builders shared by the entity pages
import { esc, fmt } from "../format.js";

export function pageHead({ kicker, title, id, idHref, links = [], after = "" }) {
  return `
    <div class="page-head">
      <div class="titles">
        <div class="kicker">${esc(kicker)}</div>
        <h1>${title}</h1>
        ${id ? `<div class="id muted" style="margin-top:4px">${idHref ? `<a href="${idHref}" target="_blank" rel="noopener">${esc(id)}</a>` : esc(id)}</div>` : ""}
        ${links.length ? `<div class="links">${links.filter(l => l && l.href).map(l => `<a href="${l.href}" target="_blank" rel="noopener">${esc(l.label)} ↗</a>`).join("")}</div>` : ""}
      </div>
      ${after}
    </div>`;
}

export function statRow(stats) {
  return `<div class="stat-row">${stats.map(([n, l, cls]) => `<div class="stat ${cls || ""}"><div class="n">${typeof n === "number" ? fmt.int(n) : n}</div><div class="l">${esc(l)}</div></div>`).join("")}</div>`;
}

export function card(title, bodyId, { sub = "", tools = "", cls = "" } = {}) {
  return `
    <div class="card ${cls}">
      <div class="card-head"><h2>${esc(title)}</h2>${tools || (sub ? `<span class="sub">${sub}</span>` : "")}</div>
      <div id="${bodyId}"></div>
    </div>`;
}

export function notFound(kind, id) {
  return `<div class="card"><h3>No ${esc(kind)} found for <span class="id">${esc(id)}</span></h3>
          <p class="muted">It may not be part of this release. Try the search box.</p></div>`;
}
