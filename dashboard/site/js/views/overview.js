// The release page: what this month's graph gave us, species by species
import { state, fmt, esc, speciesOf, shortBinomial, route } from "../format.js";
import { speciesTree } from "../charts/tree.js";
import { barList, speciesBars, funnel } from "../charts/bars.js";

// One card per species statistic: [key, title, what it counts]
const METRICS = [
  ["n_phenologs", "Phenologs", "significant disease–phenotype pairs"],
  ["n_diseases", "Diseases", "human diseases with at least one phenolog"],
  ["n_phenotypes_matched", "Phenotypes with diseases", "model phenotypes in at least one phenolog"],
  ["n_candidates", "Candidate pairs", "gene–disease candidates its phenologs support"],
  ["n_phenotyped_genes", "Genes with phenotypes", "annotated genes the test starts from"],
  ["n_common_with_base", "Families shared with human", "PANTHER families both species carry"],
];

// The species table's counts, as columns of the graph inputs table
const INPUT_COLUMNS = [
  ["n_genes", "Genes"],
  ["n_phenotyped_genes", "With phenotypes"],
  ["n_ortholog_genes", "With orthologs"],
  ["n_both", "With both"],
  ["n_phenotypes", "Phenotypes"],
  ["n_orthogroups", "Protein families"],
  ["n_g2p_edges", "Phenotype annotations"],
];

const TREE_METRICS = {
  n_phenologs: "phenologs",
  n_diseases: "diseases with a phenolog",
  n_phenotyped_genes: "genes with phenotypes",
  n_phenotypes_annotated: "annotated phenotypes",
  n_common_with_base: "protein families shared with human",
  n_candidates: "candidate pairs supported",
};

export async function render(el) {
  const ov = state.overview, r = ov.release, t = ov.tiles;
  const base = ov.species.find(s => s.is_base);
  const model = ov.species.filter(s => !s.is_base);
  const trials = ov.fdr.length ? ov.fdr[0].n_random_trials : "";

  el.innerHTML = `
    <section class="hero">
      <div class="card">
        <div class="kicker">Release</div>
        <h1>Phenologs from the Monarch KG <span class="mono" style="color:var(--accent)">${esc(r.graph_version)}</span></h1>
        <p class="lede">Orthologous phenotypes between <i>${esc(base.binomial)}</i> diseases and the phenotypes of
          ${model.length} model species, and the gene candidates they imply. Overlap is measured on PANTHER protein
          families, significance comes from ${fmt.int(trials)} randomised trials at an FDR of ${r.fdr}, and candidates are
          ranked by the ${r.rank_metric === "hg" ? "hypergeometric" : "naive Bayes"} score over the k=${fmt.int(r.kneighbs)} nearest phenologs.</p>
        <div class="meta">
          <span>graph <b>${esc(r.graph_version)}</b></span>
          ${r.phenio_version ? `<span>phenio <b>${esc(r.phenio_version)}</b></span>` : ""}
          <span>fdr <b>${r.fdr}</b></span><span>k <b>${r.kneighbs}</b></span><span>metric <b>${esc(r.rank_metric)}</b></span>
          <span>built <b>${esc(r.built)}</b></span>
        </div>
      </div>
      <div class="tiles">
        ${tile(t.diseases, "diseases with phenologs", "accent")}
        ${tile(t.phenologs, "significant phenologs", "accent")}
        ${tile(t.model_phenotypes, "model phenotypes with diseases")}
        ${tile(t.model_species, "model species")}
        ${tile(t.known_recovered, `known gene–disease pairs recovered (${fmt.pct(t.known_recovered / t.known_in_world)} of ${fmt.int(t.known_in_world)} in reach)`, "", "pipeline")}
        ${tile(t.candidate_genes, "genes with one or more candidate disease associations")}
        ${tile(t.families_implicated, "protein families implicated")}
        ${tile(t.candidate_pairs, "ranked gene–disease pairs")}
      </div>
    </section>

    <section class="section">
      <div class="card-head"><h2>Model species</h2>
        <span class="sub">hover a species to follow it across the cards · click for its page
          <span class="chips" id="metric-scale"><span class="chip toggle" data-scale="linear">linear</span><span class="chip toggle off" data-scale="log">log</span></span></span></div>
      <div class="grid cols-3 metrics" id="metrics">
        ${METRICS.map(([key, title, sub]) => metricCard(title, sub, "m-" + key)).join("")}
      </div>
    </section>

    <section class="section grid cols-2">
      <div class="card">
        <div class="card-head"><h2>Species tree</h2>
          <select id="tree-metric">${Object.entries(TREE_METRICS).map(([k, v]) => `<option value="${k}">${v}</option>`).join("")}</select></div>
        <div id="tree" style="overflow:auto"></div>
      </div>
      <div class="card">
        <div class="card-head"><h2>False discovery control</h2><span class="sub">p-value cutoff per species pair at FDR ${r.fdr}</span></div>
        <div id="fdr"></div>
      </div>
    </section>

    <section class="section grid cols-2">
      ${listCard("Diseases with the most phenologs", "top-diseases", "#/browse/diseases")}
      ${listCard("Model phenotypes matching the most diseases", "top-phenotypes", "#/browse/phenotypes", "phenotypes carried by at least 5 protein families")}
      ${listCard("Human genes predicted for the most diseases", "top-genes", "#/browse/genes")}
      ${listCard("Protein families behind the most predictions", "top-families", "#/browse/families")}
    </section>

    <section class="section" id="pipeline">
      <div class="card-head"><h2>From graph to phenologs</h2><span class="sub">what the knowledge graph gave the pipeline, and how far its curated disease genes get</span></div>
      <div class="card"><div class="card-head"><h2>Species in the graph</h2><span class="sub">the pipeline's species table · a species takes part when it has genes with both a phenotype and an ortholog</span></div>
        <div id="inputs"></div></div>
      <div class="grid cols-2" style="margin-top:16px">
        <div class="card"><div class="card-head"><h2>Known gene–disease pairs</h2><span class="sub">curated for ${esc(base.binomial)}</span></div><div id="funnel-pairs"></div></div>
        <div class="card"><div class="card-head"><h2>Diseases with a known gene</h2><span class="sub">and how many find a phenolog</span></div><div id="funnel-diseases"></div></div>
      </div>
    </section>`;

  // metric cards: the same species, in the same order and colours, in every card
  const metrics = el.querySelector("#metrics");
  const drawMetrics = log => METRICS.forEach(([key]) => speciesBars(metrics.querySelector("#m-" + key), speciesItems(model, key), { log }));
  drawMetrics(false);
  el.querySelectorAll("#metric-scale .chip").forEach(chip => chip.addEventListener("click", () => {
    el.querySelectorAll("#metric-scale .chip").forEach(c => c.classList.toggle("off", c !== chip));
    drawMetrics(chip.dataset.scale === "log");
  }));
  followSpecies(metrics);

  // tree
  const treeEl = el.querySelector("#tree");
  const sel = el.querySelector("#tree-metric");
  const drawTree = () => speciesTree(treeEl, ov.tree, ov.species, { metric: sel.value, metricLabel: TREE_METRICS[sel.value], width: Math.max(620, treeEl.clientWidth - 4) });
  sel.addEventListener("change", drawTree);
  drawTree();

  // fdr table
  const fdrRows = ov.fdr.filter(f => f.fdr === r.fdr).sort((a, b) => b.phenolog_count - a.phenolog_count);
  el.querySelector("#fdr").innerHTML = `<div class="tbl-wrap"><table class="tbl"><thead><tr>
      <th>species</th><th class="num">p cutoff</th><th class="num">phenologs</th><th class="num">expected false</th><th class="num">realised fdr</th></tr></thead><tbody>
    ${fdrRows.map(f => `<tr><td class="nowrap"><a href="${route.species(f.species)}"><span class="dot" style="background:${speciesOf(f.species).color}"></span> <i>${esc(shortBinomial(f.species))}</i></a></td>
      <td class="num mono">${fmt.p(f.pvalue_cutoff)}</td><td class="num">${fmt.int(f.phenolog_count)}</td>
      <td class="num">${Number(f.expected_false_positives).toFixed(1)}</td><td class="num">${fmt.pct(f.realized_fdr)}</td></tr>`).join("")}
    </tbody></table></div>`;

  // from graph to phenologs
  inputsTable(el.querySelector("#inputs"), ov.inputs);
  const k = ov.known_funnel;
  funnel(el.querySelector("#funnel-pairs"), [
    { label: "curated in the graph", n: k.curated },
    { label: "gene has ortholog", n: k.has_family },
    { label: "gene has ≥ 1 ortholog annotated with a phenotype <span class=\"faint\">(in reach)</span>", n: k.in_world },
    { label: "disease has a significant phenolog", n: k.has_phenolog },
    { label: "gene–disease pairs recovered by phenologs method", n: k.recovered },
  ]);
  funnel(el.querySelector("#funnel-diseases"), [
    { label: "diseases with a curated gene", n: k.diseases },
    { label: "with a gene that has ≥ 1 ortholog across species", n: k.diseases_family },
    { label: "with a gene that has ≥ 1 ortholog annotated with a phenotype", n: k.diseases_world },
    { label: "with a significant phenolog", n: k.diseases_phenolog },
  ], { color: "var(--accent-2)" });
  el.querySelectorAll(".tile[data-target]").forEach(tl => tl.addEventListener("click", () => {
    document.getElementById(tl.dataset.target).scrollIntoView({ behavior: "smooth" });
  }));

  // top lists
  const top = ov.top;
  barList(el.querySelector("#top-diseases"), top.diseases.map(d => ({ label: d.label, value: d.value, href: route.disease(d.id), sub: `${d.n_species} sp.`, color: speciesOf(d.top_species).color })));
  barList(el.querySelector("#top-phenotypes"), top.phenotypes.map(d => ({ label: d.label, value: d.value, href: route.phenotype(d.id), sub: `<i>${esc(speciesOf(d.species).binomial)}</i>`, color: speciesOf(d.species).color })));
  barList(el.querySelector("#top-genes"), top.genes.map(d => ({ label: d.label, value: d.value, href: route.gene(d.id), sub: `${d.n_known_diseases} known · ${esc(d.family || "")}` })));
  barList(el.querySelector("#top-families"), top.families.map(d => ({ label: d.label, value: d.value, href: route.family(d.id), sub: esc((d.base_symbols || "").split(", ").slice(0, 3).join(", ")) })));
}

function tile(n, label, cls = "", target = "") {
  return `<div class="tile ${cls} ${target ? "link" : ""}" ${target ? `data-target="${target}"` : ""}><div class="n">${fmt.int(n)}</div><div class="l">${label}</div></div>`;
}

function metricCard(title, sub, id, cls = "") {
  return `<div class="card metric ${cls}"><h2>${title}</h2><div class="sub">${sub}</div><div id="${id}"></div></div>`;
}

function speciesItems(species, key, format = fmt.int) {
  return species.map(s => ({ species: s.species, label: shortBinomial(s.species), value: s[key] || 0, text: format(s[key]),
                             color: s.color, href: route.species(s.species) }));
}

// Hovering a species in any card lights it up in all of them
function followSpecies(root) {
  root.addEventListener("mouseover", e => {
    const sp = e.target.closest("[data-sp]")?.dataset.sp;
    root.classList.toggle("hl", sp != null);
    root.querySelectorAll("[data-sp]").forEach(n => n.classList.toggle("on", n.dataset.sp === sp));
  });
  root.addEventListener("mouseleave", () => root.classList.remove("hl"));
}

function listCard(title, id, browseHref, sub = "") {
  return `<div class="card"><div class="card-head"><h2>${title}</h2><span class="sub">${sub ? sub + " · " : ""}<a href="${browseHref}">browse all</a></span></div><div id="${id}"></div></div>`;
}

// What the graph gave the pipeline: every taxon of the species table, with a hairline bar under each count
// (scaled within its column). Species left out are dimmed with the reason extraction gave
function inputsTable(el, rows) {
  const cols = INPUT_COLUMNS.filter(([key]) => rows.some(r => r[key] != null));    // older species tables have no "with both"
  const max = Object.fromEntries(cols.map(([key]) => [key, Math.max(1, ...rows.map(r => r[key] || 0))]));
  el.innerHTML = `<div class="tbl-wrap"><table class="tbl inputs"><thead><tr><th>Species</th>
      ${cols.map(([, label]) => `<th class="num">${label}</th>`).join("")}</tr></thead><tbody>
    ${rows.map(r => `<tr class="${r.included ? "" : "excluded"}" ${r.notes ? `title="${esc(r.notes)}"` : ""}>
      <td>${r.included ? `<a href="${route.species(r.species)}">` : ""}<span class="dot" style="background:${r.color}"></span> <i>${esc(r.binomial)}</i>${r.included ? "</a>" : ""}
        <div class="faint">${r.included ? esc(r.common) : `left out · ${leftOut(r)}`}</div></td>
      ${cols.map(([key]) => `<td class="num">${fmt.int(r[key])}<span class="hair" style="width:${(100 * (r[key] || 0) / max[key]).toFixed(1)}%; background:${r.color}"></span></td>`).join("")}
    </tr>`).join("")}</tbody></table></div>`;
}

// Why a species sat out, in a few words; the pipeline's own note is the row's tooltip
function leftOut(r) {
  if (r.n_phenotyped_genes === 0) return r.n_ortholog_genes === 0 ? "no phenotype annotations or orthologs" : "no phenotype annotations";
  if (r.n_ortholog_genes === 0) return "no orthologs";
  return esc(r.notes || "too few genes with both a phenotype and an ortholog");
}
