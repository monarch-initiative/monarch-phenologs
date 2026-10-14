// Browse pages: every disease / gene / phenotype / family / species as a sortable table
import { ensure, query, lit } from "../db.js";
import { state, fmt, esc, speciesChip, speciesOf, route } from "../format.js";
import { table } from "../charts/table.js";

const num = (key, label, extra = {}) => ({ key, label, num: true, ...extra });
const p = v => `<span class="mono">${fmt.p(v)}</span>`;

const KINDS = {
  diseases: {
    title: "Diseases", sub: "Every human disease with at least one significant phenolog",
    sql: () => `SELECT * FROM diseases`,
    columns: [
      { key: "name", label: "Disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}">${esc(v)}</a>` },
      { key: "disease_id", label: "Id", format: v => `<span class="id muted">${esc(v)}</span>` },
      num("n_known_genes", "Known genes"), num("n_phenologs", "Phenologs"), num("n_species", "Species"),
      num("n_known_recovered", "Recovered"), num("n_candidates", "Ranked genes"),
      { key: "top_species", label: "Top species", format: v => v ? speciesChip(v) : "" },
      num("best_pval", "Best p", { format: p }),
    ],
    sort: "n_phenologs", search: ["name", "disease_id"], placeholder: "Filter by disease name or Mondo id…",
  },
  genes: {
    title: "Genes", sub: "Genes with phenotypes, a curated disease, or a candidate",
    sql: () => `SELECT * FROM genes WHERE n_phenotypes > 0 OR n_known_diseases > 0 OR n_predicted_diseases > 0`,
    columns: [
      { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}">${esc(v)}</a>` },
      { key: "gene_id", label: "Id", format: v => `<span class="id muted">${esc(v)}</span>` },
      { key: "species", label: "Species", format: v => speciesChip(v) },
      { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : "" },
      num("n_ortholog_species", "Species with orthologs", { help: "Other species of this release the gene has an ortholog in" }),
      num("n_phenotypes", "Phenotypes"), num("n_known_diseases", "Known diseases"), num("n_predicted_diseases", "Ranked for"),
      num("n_novel_diseases", "Novel"), num("best_rank", "Best rank"),
    ],
    sort: "n_predicted_diseases", search: ["symbol", "gene_id", "family"], placeholder: "Filter by gene name or id…",
    filters: true,
  },
  phenotypes: {
    title: "Phenotypes", sub: "All species phenotypes",
    sql: () => `SELECT * FROM phenotypes`,
    columns: [
      { key: "name", label: "Phenotype", format: (v, r) => `<a href="${route.phenotype(r.phenotype_id)}">${esc(v)}</a>` },
      { key: "phenotype_id", label: "Id", format: v => `<span class="id muted">${esc(v)}</span>` },
      { key: "species", label: "Species", format: v => speciesChip(v) },
      num("n_genes", "Genes"), num("n_orthogroups", "Families"),
      num("n_known_diseases", "Known diseases", { help: "Human diseases with a curated gene in one of this phenotype's protein families" }),
      num("n_diseases", "Phenolog disease associations", { help: "Human diseases this phenotype is a significant phenolog of" }),
    ],
    sort: "n_diseases", search: ["name", "phenotype_id"], placeholder: "Filter by phenotype name or id…",
    filters: true,
  },
  families: {
    title: "Protein families", sub: "PANTHER families with at least one gene in the release",
    sql: () => `SELECT * FROM families`,
    columns: [
      { key: "family", label: "Family", format: v => `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` },
      { key: "base_symbols", label: "Human genes", ellipsis: true, format: v => esc(v || "") },
      num("n_genes", "Genes"), num("n_species", "Species"), num("n_known_diseases", "Known diseases"),
      num("n_predicted_diseases", "Ranked for"), num("n_novel_diseases", "Novel"), num("best_score", "Best score", { format: p }),
    ],
    sort: "n_predicted_diseases", search: ["family", "base_symbols"], placeholder: "Filter by family or human gene…",
  },
  species: {
    title: "Species", sub: "The species of this release",
    sql: () => `SELECT * FROM species ORDER BY plot_order`,
    columns: [
      { key: "binomial", label: "Species", format: (v, r) => `<a href="${route.species(r.species)}"><span class="dot" style="background:${r.color}"></span> <i>${esc(v)}</i></a> <span class="muted">${esc(r.common)}</span>` },
      { key: "taxon_id", label: "Taxon", format: v => `<span class="id muted">${esc(v)}</span>` },
      { key: "ontology", label: "Ontology" },
      num("n_genes", "Genes"), num("n_phenotyped_genes", "With phenotypes"), num("n_phenotypes_annotated", "Phenotypes"),
      num("n_orthogroups", "Families"), num("n_common_with_base", "Shared w/ human"), num("n_phenologs", "Phenologs"),
      num("n_diseases", "Diseases"), num("n_candidates", "Candidate pairs"), num("fdr_cutoff", "p cutoff", { format: p }),
    ],
    sort: "n_phenologs", search: ["binomial", "common"], placeholder: "Filter by species…",
  },
};

export async function render(el, kind) {
  const cfg = KINDS[kind];
  if (!cfg) { el.innerHTML = `<div class="card">Unknown list</div>`; return; }
  const tableName = kind;
  await ensure(tableName);
  const rows = await query(cfg.sql());

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Browse</div><h1>${esc(cfg.title)}</h1><div class="muted">${esc(cfg.sub)}</div></div></div>
    <div class="card" id="list"></div>`;

  // species filter chips for the big lists
  let active = null;
  const speciesPresent = cfg.filters ? state.overview.species.map(s => s.species).filter(s => rows.some(r => r.species === s)) : [];
  const tools = speciesPresent.length ? `<span class="chips" id="sp-filter">${speciesPresent.map(s => `<span class="chip toggle" data-sp="${s}"><span class="dot" style="background:${speciesOf(s).color}"></span><i>${esc(speciesOf(s).binomial)}</i></span>`).join("")}</span>` : "";

  const t = table(el.querySelector("#list"), rows, cfg.columns, {
    sortKey: cfg.sort, sortDir: "desc", search: cfg.search, placeholder: cfg.placeholder, pageSize: 50, tools,
    onDraw: root => root.querySelectorAll("#sp-filter .chip").forEach(chip => {
      chip.classList.toggle("off", active != null && active !== chip.dataset.sp);
      chip.addEventListener("click", () => {
        active = active === chip.dataset.sp ? null : chip.dataset.sp;
        t.setRows(active ? rows.filter(r => r.species === active) : rows);
      });
    }),
  });
}
