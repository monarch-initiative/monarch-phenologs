// Formatting, species lookups, the open network's words and outbound links shared by every view.

// A release holds one prediction network, named for the human side of its phenologs. Its tables call those terms
// either way (term_id, terms); the pages use these words and routes
export const NETWORKS = {
  disease: { key: "disease", term: "disease", terms: "diseases", Term: "Disease", Terms: "Diseases", model: "phenotype",
             human: "Human disease", humans: "Human diseases", pair: "gene–disease", phenolog: "disease–phenotype",
             title: "Disease phenologs", ontology: "Mondo", route: "disease", browse: "diseases" },
  phenotype: { key: "phenotype", term: "human phenotype", terms: "human phenotypes", Term: "Human phenotype", Terms: "Human phenotypes",
               model: "model phenotype", human: "Human phenotype", humans: "Human phenotypes", pair: "gene–phenotype",
               phenolog: "human–model phenotype", title: "Phenotype phenologs", ontology: "HPO", route: "human-phenotype",
               browse: "human-phenotypes" },
};

// The same page in another network's release (a reload: the network is part of the address, ?network=)
export function networkHref(network, hash = location.hash) {
  const p = new URLSearchParams(location.search);
  p.set("network", network);
  return `${location.pathname}?${p}${hash}`;
}

export const state = { overview: null, species: new Map(), net: NETWORKS.disease };

export function setOverview(ov) {
  state.overview = ov;
  state.species = new Map(ov.species.map(s => [s.species, s]));
  state.net = NETWORKS[ov.release.prediction_network] || NETWORKS.disease;
}

export const fmt = {
  int: n => (n == null ? "–" : Number(n).toLocaleString("en-US")),
  pct: x => (x == null ? "–" : (100 * x).toFixed(1) + "%"),
  p: p => {
    if (p == null) return "–";
    if (p === 0) return "0";
    if (p < 1e-3) return Number(p).toExponential(1).replace("e-", "e−");
    return Number(p).toPrecision(2);
  },
  score: s => (s == null ? "–" : Number(s).toFixed(3)),
  dec: x => (x == null ? "–" : Number(x).toLocaleString("en-US", { minimumFractionDigits: 1, maximumFractionDigits: 1 })),
};

export function speciesOf(name) {
  return state.species.get(name) || { species: name, binomial: name.replace(/-/g, " "), common: "", color: "#7f7f7f" };
}

export function speciesColor(name) {
  return speciesOf(name).color;
}

export function binomial(name) {
  return speciesOf(name).binomial;
}

// "Mus-musculus" --> "M. musculus", for tight chart labels
export function shortBinomial(name) {
  const [genus, ...rest] = speciesOf(name).binomial.split(" ");
  return rest.length ? `${genus[0]}. ${rest.join(" ")}` : genus;
}

// Outbound links. dismech() is the hook for the DisMech knowledge base once we know which
// disorders it covers; it returns null until then and the views skip it.
export const links = {
  monarch: id => `https://monarchinitiative.org/${id}`,
  panther: id => `https://www.pantherdb.org/panther/family.do?clsAccession=${id}`,
  dismech: () => null,
};

export const route = {
  term: id => `#/${state.net.route}/${id}`,
  gene: id => `#/gene/${id}`,
  phenotype: id => `#/phenotype/${id}`,
  family: id => `#/family/${id}`,
  species: name => `#/species/${name}`,
};

export function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

export function speciesChip(name, { count, href } = {}) {
  const s = speciesOf(name);
  const inner = `<span class="dot" style="background:${s.color}"></span><i>${esc(s.binomial)}</i>${count != null ? ` <span class="muted">${fmt.int(count)}</span>` : ""}`;
  return `<a class="chip" href="${href || route.species(name)}">${inner}</a>`;
}

// A candidate's support, "species:phenologs:of the term's:orthologs:hits,..." (most phenologs first)
export function parseSupport(s) {
  if (!s) return [];
  return s.split(",").map(item => {
    const [species, phenologs, used, orthologs, hits] = item.split(":");
    return { species, phenologs: +phenologs, used: +used, orthologs: +orthologs, hits: +hits };
  });
}

// One chip per species: how many of the term's phenologs with that species are phenotypes of an ortholog
// of the gene, out of how many the ranking used. The tooltip adds how many orthologs that takes
export function supportChips(s, gene) {
  return parseSupport(s).map(x => {
    const sp = speciesOf(x.species);
    const title = `${fmt.int(x.phenologs)} of ${fmt.int(x.used)} ${shortBinomial(x.species)} phenologs, ` +
                  `through ${fmt.int(x.orthologs)} ortholog${x.orthologs === 1 ? "" : "s"} of ${gene}`;
    return `<a class="chip" href="${route.species(x.species)}" title="${esc(title)}"><span class="dot" style="background:${sp.color}"></span>` +
           `<i>${esc(shortBinomial(x.species))}</i> <span>${fmt.int(x.phenologs)}<span class="faint">/${fmt.int(x.used)}</span></span></a>`;
  }).join(" ");
}
