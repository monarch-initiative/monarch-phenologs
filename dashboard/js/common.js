// Small HTML builders shared by the entity pages
import { state, esc, fmt } from "./format.js";

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

// Whether a known pair was recovered, for the known term and known gene grids (each page adds a third key)
export const recoveryKey = () => [["var(--accent)", `known ${state.net.pair} pair recovered`],
                                  ["var(--text)", `${state.net.pair} pair not recovered by associated phenologs`]];

// The two kinds of gene-term pair every ranking shows
export const pairKey = () => [["var(--accent)", `known ${state.net.pair} pair`], ["var(--accent-2)", "novel candidate pair"]];

// The column every candidate table shares, and its header tooltip below
export const supportLabel = () => `Species phenologs supporting this ${state.net.pair} pair`;

// Header tooltips for the columns every candidate table shares
export const HELP = {
  get score() {
    return "The candidate ranking's hypergeometric p-value: how unlikely this many phenotype annotations of the gene's orthologs are " +
           `among the ${state.net.term}'s nearest phenologs. Smaller ranks higher. A ranking score, not a corrected significance level`;
  },
  get support() { return `x/y: x of the ${state.net.term}'s y phenologs in that species are phenotypes of an ortholog of this gene`; },
  get rank() { return `Place in the ${state.net.term}'s ranking by score (ties share a place), out of every gene it scored, known genes included`; },
};

// Header tooltips for the FDR tables (release and species pages)
export const FDR_HELP = {
  cutoff: "Largest p-value still called a phenolog at this FDR",
  get phenologs() { return `Real ${state.net.phenolog} comparisons at or below the cutoff`; },
  expected: "Randomised comparisons at or below the cutoff, per trial (mean over the trials)",
  realised: "Expected false positives / phenologs",
};

// Header tooltips for the phenolog tables (term, phenotype and species pages)
export const PHENOLOG_HELP = {
  get rank() { return `Place among the ${state.net.term}'s phenologs, by p-value`; },
  get overlap() { return `Protein families carrying both the ${state.net.term} and the ${state.net.model}`; },
  get n_term_orthogroups() { return `Protein families of the ${state.net.term}'s genes, counting only families the two species share`; },
  get n_phenotype_orthogroups() { return `Protein families of the ${state.net.model}'s genes, counting only families the two species share`; },
  n_common: "Protein families present in both species: the background of the test",
  pval: "Hypergeometric p-value of the overlap; a phenolog is significant at or below its species' FDR cutoff",
};

// Names in left-aligned columns, as many as fit: [{ label, href, cls, note, title, ext }]. wide: columns for longer
// names. ext: the link leaves the site (opens a new tab)
export function nameGrid(items, { wide = false } = {}) {
  return `<div class="name-grid ${wide ? "wide" : ""}">${items.map(i =>
    `<a class="${i.cls || ""}" href="${i.href}" title="${esc(i.title || i.label)}" ${i.ext ? `target="_blank" rel="noopener"` : ""}>` +
    `${esc(i.label)}${i.note ? ` <span class="faint">${esc(i.note)}</span>` : ""}</a>`).join("")}</div>`;
}
