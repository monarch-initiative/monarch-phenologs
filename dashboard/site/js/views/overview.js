// The release page: what this month's graph gave us, species by species
import { state, fmt, esc, speciesOf, shortBinomial, route } from "../format.js";
import { speciesTree } from "../charts/tree.js";
import { barList, speciesBars, funnel } from "../charts/bars.js";
import { FDR_HELP } from "./common.js";

// One card per species statistic: [key, title, what it counts]
const METRICS = [
  ["n_phenologs", "Phenologs", "significant disease–phenotype pairs"],
  ["n_diseases", "Diseases", "human diseases with at least one phenolog"],
  ["n_phenotypes_matched", "Phenotypes with diseases", "model phenotypes in at least one phenolog"],
  ["n_novel_candidates", "Candidate pairs", "novel gene–disease pairs with a supporting phenolog from this species"],
  ["n_phenotyped_genes", "Genes with phenotypes", "genes with at least one phenotype annotation"],
  ["n_common_with_base", "Families shared with human", "PANTHER families both species carry"],
];

// The species table's counts, as columns of the graph inputs table: [key, header, header tooltip]
const INPUT_COLUMNS = [
  ["n_genes", "Genes", "Gene nodes of the taxon in the graph"],
  ["n_phenotyped_genes", "With phenotypes", "Genes with at least one phenotype annotation"],
  ["n_ortholog_genes", "With orthologs", "Genes with an ortholog edge to any species in the graph, which also gives the gene its PANTHER family"],
  ["n_both", "With both", "Genes with a phenotype and an ortholog: what a phenolog comparison is made from"],
  ["n_phenotypes", "Phenotypes on genes with an ortholog", "Distinct phenotypes annotated to at least one gene that has an ortholog"],
  ["n_orthogroups", "Protein families", "PANTHER families of the taxon's genes"],
  ["n_g2p_edges", "Phenotype annotations", "Gene–phenotype edges"],
];

const TREE_METRICS = {
  n_phenologs: "phenologs",
  n_diseases: "diseases with a phenolog",
  n_phenotyped_genes: "genes with phenotypes",
  n_phenotypes_annotated: "annotated phenotypes",
  n_common_with_base: "protein families shared with human",
  n_novel_candidates: "novel candidate pairs supported",
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
          families, each species' p-value cutoff is set from ${fmt.int(trials)} randomised trials at an FDR of ${r.fdr}, and candidates are
          ranked by the ${r.rank_metric === "hg" ? "hypergeometric" : "naive Bayes"} score over each disease's ${fmt.int(r.kneighbs)} most significant phenologs.</p>
        <div class="meta">
          <span>graph <b>${esc(r.graph_version)}</b></span>
          ${r.phenio_version ? `<span>phenio <b>${esc(r.phenio_version)}</b></span>` : ""}
          <span>fdr <b>${r.fdr}</b></span><span>k <b>${r.kneighbs}</b></span><span>metric <b>${esc(r.rank_metric)}</b></span>
          <span>built <b>${esc(r.built)}</b></span>
        </div>
      </div>
      <div class="tiles">
        ${tile(t.diseases, "diseases with phenologs", { cls: "accent", help: "Human diseases with at least one significant phenolog" })}
        ${tile(t.phenologs, "significant phenologs", { cls: "accent", help: "Disease–phenotype pairs at or below their species' FDR cutoff" })}
        ${tile(t.model_phenotypes, "model phenotypes with diseases", { help: "Model species phenotypes in at least one phenolog (a term used by two species counts twice)" })}
        ${tile(t.model_species, "model species", { help: "Species besides human taking part in this release" })}
        ${tile(t.known_recovered, `known gene–disease pairs recovered (${fmt.pct(t.known_recovered / t.known_in_reach)} of ${fmt.int(t.known_in_reach)} in reach)`,
               { target: "pipeline", help: "Known pairs whose gene is in the disease's candidate ranking, out of the pairs in reach (the gene has an ortholog with phenotypes and the disease has a phenolog). Click for each step" })}
        ${tile(t.candidate_genes, "genes with one or more candidate disease associations", { help: "Human genes ranked for at least one disease they are not known for" })}
        ${tile(t.families_implicated, "protein families with a ranked gene", { help: "PANTHER families with a human gene in at least one disease's ranking" })}
        ${tile(t.candidate_pairs, "ranked gene–disease pairs", { help: "Every (disease, gene) pair the ranking scored, known pairs included" })}
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
      ${listCard("Human genes ranked for the most diseases", "top-genes", "#/browse/genes")}
      ${listCard("Protein families ranked for the most diseases", "top-families", "#/browse/families")}
    </section>

    <section class="section" id="pipeline">
      <div class="card-head"><h2>From graph to phenologs</h2><span class="sub">what the knowledge graph gave the pipeline, and how far its curated disease genes get</span></div>
      <div class="card"><div class="card-head"><h2>Species in the graph</h2><span class="sub">the pipeline's species table · a species takes part when enough of its genes have both a phenotype and an ortholog</span></div>
        <div id="inputs"></div></div>
      <div class="grid cols-2" style="margin-top:16px">
        <div class="card"><div class="card-head"><h2>Known gene–disease pairs</h2><span class="sub">curated for ${esc(base.binomial)} · % of the step above</span></div><div id="funnel-pairs"></div></div>
        <div class="card"><div class="card-head"><h2>Diseases with a known gene</h2><span class="sub">and how far they get · % of the step above</span></div><div id="funnel-diseases"></div></div>
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
  el.querySelector("#fdr").innerHTML = `<div class="tbl-wrap"><table class="tbl wrap-head"><thead><tr>
      <th>species</th><th class="num" title="${esc(FDR_HELP.cutoff)}">p cutoff</th><th class="num" title="${esc(FDR_HELP.phenologs)}">phenologs</th>
      <th class="num" title="${esc(FDR_HELP.expected)}">expected false positives</th><th class="num" title="${esc(FDR_HELP.realised)}">realised FDR</th></tr></thead><tbody>
    ${fdrRows.map(f => `<tr><td class="nowrap"><a href="${route.species(f.species)}"><span class="dot" style="background:${speciesOf(f.species).color}"></span> <i>${esc(shortBinomial(f.species))}</i></a></td>
      <td class="num mono">${fmt.p(f.pvalue_cutoff)}</td><td class="num">${fmt.int(f.phenolog_count)}</td>
      <td class="num">${fmt.dec(f.expected_false_positives)}</td><td class="num">${fmt.pct(f.realized_fdr)}</td></tr>`).join("")}
    </tbody></table></div>`;

  // from graph to phenologs
  inputsTable(el.querySelector("#inputs"), ov.inputs);
  const k = ov.known_funnel;
  funnel(el.querySelector("#funnel-pairs"), [
    { label: "curated in the graph", n: k.curated },
    { label: "gene has an ortholog in a model species", n: k.has_ortholog },
    { label: "gene has ≥ 1 ortholog annotated with a phenotype", n: k.gene_in_reach },
    { label: "disease has a significant phenolog <span class=\"faint\">(in reach)</span>", n: k.in_reach },
    { label: "gene–disease pairs recovered by phenologs method", n: k.recovered },
  ]);
  funnel(el.querySelector("#funnel-diseases"), [
    { label: "diseases with a curated gene", n: k.diseases },
    { label: "with a gene that has an ortholog in a model species", n: k.diseases_ortholog },
    { label: "with a gene that has ≥ 1 ortholog annotated with a phenotype", n: k.diseases_gene_in_reach },
    { label: "with a significant phenolog <span class=\"faint\">(in reach)</span>", n: k.diseases_in_reach },
    { label: "with a known gene recovered", n: k.diseases_recovered },
  ], { color: "var(--accent-2)" });
  el.querySelectorAll(".tile[data-target]").forEach(tl => tl.addEventListener("click", () => {
    document.getElementById(tl.dataset.target).scrollIntoView({ behavior: "smooth" });
  }));

  // top lists
  const top = ov.top;
  barList(el.querySelector("#top-diseases"), top.diseases.map(d => ({ label: d.label, value: d.value, href: route.disease(d.id), sub: `${d.n_species} species`, color: speciesOf(d.top_species).color })));
  barList(el.querySelector("#top-phenotypes"), top.phenotypes.map(d => ({ label: d.label, value: d.value, href: route.phenotype(d.id), sub: `<i>${esc(speciesOf(d.species).binomial)}</i>`, color: speciesOf(d.species).color })));
  barList(el.querySelector("#top-genes"), top.genes.map(d => ({ label: d.label, value: d.value, href: route.gene(d.id), sub: `${d.n_known_diseases} known disease${d.n_known_diseases === 1 ? "" : "s"} · ${esc(d.family || "")}` })));
  barList(el.querySelector("#top-families"), top.families.map(d => {
    const genes = (d.base_symbols || "").split(", ").filter(Boolean);
    return { label: d.label, value: d.value, href: route.family(d.id),
             sub: esc(genes.slice(0, 3).join(", ")) + (genes.length > 3 ? ` <span class="faint">+${genes.length - 3}</span>` : "") };
  }));
}

function tile(n, label, { cls = "", target = "", help = "" } = {}) {
  return `<div class="tile ${cls} ${target ? "link" : ""}" ${target ? `data-target="${target}"` : ""} ${help ? `title="${esc(help)}"` : ""}>
    <div class="n">${fmt.int(n)}</div><div class="l">${label}</div></div>`;
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
  el.innerHTML = `<div class="tbl-wrap"><table class="tbl inputs wrap-head"><thead><tr><th>Species</th>
      ${cols.map(([, label, help]) => `<th class="num" title="${esc(help)}">${label}</th>`).join("")}</tr></thead><tbody>
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
