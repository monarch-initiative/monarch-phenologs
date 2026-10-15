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

// Whether a known pair was recovered, for the known disease and known gene grids (each page adds a third key)
export const RECOVERY_KEY = [["var(--accent)", "known gene–disease pair recovered"], ["var(--text)", "gene–disease pair not recovered by associated phenologs"]];

// The two kinds of gene-disease pair every ranking shows
export const PAIR_KEY = [["var(--accent)", "known gene–disease pair"], ["var(--accent-2)", "novel candidate pair"]];

// The column every candidate table shares, and its header tooltip below
export const SUPPORT_LABEL = "Species phenologs supporting this gene–disease pair";

// Header tooltips for the columns every candidate table shares
export const HELP = {
  score: "The candidate ranking's hypergeometric p-value: how unlikely this many phenotype annotations of the gene's orthologs are " +
         "among the disease's nearest phenologs. Smaller ranks higher. A ranking score, not a corrected significance level",
  support: "x/y: x of the disease's y phenologs in that species are phenotypes of an ortholog of this gene",
  rank: "Place in the disease's ranking by score (ties share a place), out of every gene it scored, known genes included",
};

// Header tooltips for the FDR tables (release and species pages)
export const FDR_HELP = {
  cutoff: "Largest p-value still called a phenolog at this FDR",
  phenologs: "Real disease–phenotype comparisons at or below the cutoff",
  expected: "Randomised comparisons at or below the cutoff, per trial (mean over the trials)",
  realised: "Expected false positives / phenologs",
};

// Header tooltips for the phenolog tables (disease, phenotype and species pages)
export const PHENOLOG_HELP = {
  rank: "Place among the disease's phenologs, by p-value",
  overlap: "Protein families carrying both the disease and the phenotype",
  n_disease_orthogroups: "Protein families of the disease's genes, counting only families the two species share",
  n_phenotype_orthogroups: "Protein families of the phenotype's genes, counting only families the two species share",
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
