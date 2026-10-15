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

export function statRow(stats, cls = "") {
  return `<div class="stat-row ${cls}">${stats.map(([n, l, cls]) => `<div class="stat ${cls || ""}"><div class="n">${typeof n === "number" ? fmt.int(n) : n}</div><div class="l">${esc(l)}</div></div>`).join("")}</div>`;
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

// A small colour key for a card head: [[colour, label], ...]
export function legend(items) {
  return `<span class="key">${items.map(([color, label]) =>
    `<span class="key-item"><span class="sw" style="background:${color}"></span>${esc(label)}</span>`).join("")}</span>`;
}

// The two kinds of gene-disease pair every ranking shows
export const PAIR_KEY = [["var(--accent)", "known gene–disease pair"], ["var(--accent-2)", "candidate pair"]];

// Names in left-aligned columns, as many as fit: [{ label, href, cls, note }]
export function nameGrid(items) {
  return `<div class="name-grid">${items.map(i =>
    `<a class="${i.cls || ""}" href="${i.href}" title="${esc(i.title || i.label)}">${esc(i.label)}${i.note ? ` <span class="faint">${esc(i.note)}</span>` : ""}</a>`).join("")}</div>`;
}
