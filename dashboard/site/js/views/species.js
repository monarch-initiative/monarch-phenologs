// One species: what the graph holds for it, its player card, its FDR control, where its phenologs land
import { ensure, query, lit } from "../db.js";
import { state, fmt, esc, route, links } from "../format.js";
import { pageHead, statRow, card, notFound } from "./common.js";
import { radar, speciesAxes } from "../charts/radar.js";
import { barList } from "../charts/bars.js";
import { table } from "../charts/table.js";

export const RADAR_AXES = [
  ["n_phenologs", "Phenologs"],
  ["n_diseases", "Diseases"],
  ["n_phenotypes_matched", "Phenotypes"],
  ["n_candidates", "Candidates"],
  ["n_common_with_base", "Shared families"],
  ["n_phenotyped_genes", "Genes w/ phenotypes"],
];

// The species' numbers, from the graph in to the phenologs out; null ones (the base species) are left out
export function speciesStats(s) {
  return [[s.n_genes, "genes in graph"], [s.n_ortholog_genes, "protein-coding genes with an ortholog (any species)"],
          [s.n_human_ortholog_genes, "protein-coding genes with a human ortholog"], [s.n_phenotyped_genes, "genes with phenotypes"],
          [s.n_phenotypes_annotated, "annotated phenotypes"], [s.n_orthogroups, "protein families"],
          [s.n_common_with_base, "families shared with human"], [s.n_phenologs, "phenologs"], [s.n_diseases, "human diseases"],
          [s.n_candidates, "candidate pairs supported"]]
    .filter(([n]) => n != null && !(s.is_base && n === 0));
}

export async function render(el, name) {
  const ov = state.overview;
  const s = ov.species.find(x => x.species === name);
  if (!s) { el.innerHTML = notFound("species", name); return; }
  const model = ov.species.filter(x => !x.is_base);
  await ensure("phenologs", "phenotypes", "diseases");

  const [topDiseases, topPhenotypes, phen] = s.is_base ? [[], [], []] : await Promise.all([
    query(`SELECT p.disease_id, any_value(p.disease_name) AS name, count(*) AS n, min(p.pval) AS best
           FROM phenologs p WHERE species = ${lit(name)} GROUP BY p.disease_id ORDER BY n DESC, best LIMIT 15`),
    query(`SELECT phenotype_id, name, n_phenologs, n_diseases, n_orthogroups FROM phenotypes
           WHERE species = ${lit(name)} AND n_phenologs > 0 ORDER BY n_phenologs DESC, name LIMIT 15`),
    query(`SELECT rank, disease_id, disease_name, phenotype_id, phenotype_name, overlap, n_disease_orthogroups, n_phenotype_orthogroups, n_common, pval
           FROM phenologs WHERE species = ${lit(name)} ORDER BY pval`),
  ]);
  const fdr = ov.fdr.filter(f => f.species === name).sort((a, b) => b.fdr - a.fdr);

  el.innerHTML = `
    ${pageHead({ kicker: "Species", title: `<i style="color:${s.color}">${esc(s.binomial)}</i>`, id: s.taxon_id, idHref: links.monarch(s.taxon_id),
                 links: [{ label: "Monarch", href: links.monarch(s.taxon_id) }],
                 after: `<div class="muted" style="align-self:center">${esc(s.common)}${s.ontology ? ` · phenotypes from <b>${esc(s.ontology)}</b>` : ""}${s.is_base ? " · base species" : ""}</div>` })}
    ${statRow(speciesStats(s), "wide")}
    ${s.is_base ? `<div class="card section"><p class="muted">The base species: its diseases are what every model species is compared against. Browse <a href="#/browse/diseases">diseases</a> or <a href="#/browse/genes">genes</a>.</p></div>` : `
    <div class="grid species-top">
      <div class="card radar-card" id="radar"></div>
      <div class="card"><div class="card-head"><h2>False discovery control</h2><span class="sub">${fmt.int(fdr[0]?.n_random_trials)} randomised trials</span></div>
        <div class="tbl-wrap"><table class="tbl"><thead><tr><th class="num">fdr</th><th class="num">p cutoff</th><th class="num">phenologs</th><th class="num">expected false</th><th class="num">realised</th></tr></thead>
        <tbody>${fdr.map(f => `<tr class="${f.fdr === ov.release.fdr ? "known-mark" : ""}"><td class="num">${f.fdr}</td><td class="num mono">${fmt.p(f.pvalue_cutoff)}</td><td class="num">${fmt.int(f.phenolog_count)}</td><td class="num">${Number(f.expected_false_positives).toFixed(1)}</td><td class="num">${fmt.pct(f.realized_fdr)}</td></tr>`).join("")}</tbody></table></div>
      </div>
    </div>
    <div class="grid cols-2 section">
      ${card("Diseases with the most phenologs here", "top-d")}
      ${card("Phenotypes matching the most diseases", "top-p")}
    </div>
    <div class="section">${card("All phenologs", "phens", { sub: `${fmt.int(phen.length)} significant phenologs with <i>${esc(s.binomial)}</i>` })}</div>`}`;

  if (s.is_base) return;
  radar(el.querySelector("#radar"), speciesAxes(s, model, RADAR_AXES), { size: 310, color: s.color, labels: true });
  barList(el.querySelector("#top-d"), topDiseases.map(d => ({ label: d.name, value: d.n, href: route.disease(d.disease_id), color: s.color })));
  barList(el.querySelector("#top-p"), topPhenotypes.map(p => ({ label: p.name, value: p.n_phenologs, href: route.phenotype(p.phenotype_id), color: s.color, sub: `${p.n_orthogroups} fam.` })));
  table(el.querySelector("#phens"), phen, [
    { key: "disease_name", label: "Human disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}">${esc(v)}</a>` },
    { key: "phenotype_name", label: "Phenotype", format: (v, r) => `<a href="${route.phenotype(r.phenotype_id)}">${esc(v)}</a>` },
    { key: "rank", label: "Rank in disease", num: true },
    { key: "overlap", label: "Overlap", num: true },
    { key: "n_disease_orthogroups", label: "Disease fam.", num: true },
    { key: "n_phenotype_orthogroups", label: "Phenotype fam.", num: true },
    { key: "pval", label: "p", num: true, format: v => `<span class="mono">${fmt.p(v)}</span>` },
  ], { sortKey: "pval", sortDir: "asc", search: ["disease_name", "phenotype_name"], placeholder: "Filter by disease or phenotype…" });
}
