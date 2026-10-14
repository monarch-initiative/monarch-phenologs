// Formatting, species lookups and outbound links shared by every view.
export const state = { overview: null, species: new Map() };

export function setOverview(ov) {
  state.overview = ov;
  state.species = new Map(ov.species.map(s => [s.species, s]));
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
  disease: id => `#/disease/${id}`,
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

// "Mus-musculus:30,Danio-rerio:5" --> [["Mus-musculus", 30], ...]
export function parseOccurrence(s) {
  if (!s) return [];
  return s.split(",").filter(Boolean).map(item => { const [sp, n] = item.split(":"); return [sp, Number(n)]; });
}

export function occurrenceChips(s) {
  return parseOccurrence(s).map(([sp, n]) => speciesChip(sp, { count: n })).join(" ");
}
