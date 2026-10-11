// The release page: what this month's graph gave us, species by species
import { state, fmt, esc, speciesOf, route } from "../format.js";
import { radar, speciesAxes } from "../charts/radar.js";
import { speciesTree } from "../charts/tree.js";
import { barList } from "../charts/bars.js";
import { showTip, moveTip, hideTip } from "../charts/tooltip.js";

export const RADAR_AXES = [
  ["n_phenologs", "Phenologs"],
  ["n_diseases", "Diseases"],
  ["n_phenotypes_matched", "Phenotypes"],
  ["n_candidates", "Candidates"],
  ["n_common_with_base", "Shared families"],
  ["n_phenotyped_genes", "Genes w/ phenotypes"],
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
        ${tile(t.model_phenotypes, "model phenotypes matched")}
        ${tile(t.model_species, "model species")}
        ${tile(t.known_recovered, `known gene–disease pairs recovered (${fmt.pct(t.known_recovered / t.known_pairs)})`)}
        ${tile(t.candidate_genes, "genes with a candidate score")}
        ${tile(t.families_implicated, "protein families implicated")}
        ${tile(t.candidate_pairs, "ranked gene–disease pairs")}
      </div>
    </section>

    <section class="section">
      <div class="card-head"><h2>Model species</h2><span class="sub">each axis log-scaled against the strongest species · click a card</span></div>
      <div class="grid cards" id="players"></div>
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
        <div class="card-head" style="margin-top:18px"><h2>Phenolog strength</h2><span class="sub">−log10 p of significant phenologs, by species</span></div>
        <div id="pvals"></div>
      </div>
    </section>

    <section class="section">
      <div class="card-head"><h2>Where the hits are</h2></div>
      <div class="grid cols-2">
        ${listCard("Diseases with the most phenologs", "top-diseases", "#/browse/diseases")}
        ${listCard("Model phenotypes matching the most diseases", "top-phenotypes", "#/browse/phenotypes", "phenotypes carried by at least 5 protein families")}
        ${listCard("Human genes predicted for the most diseases", "top-genes", "#/browse/genes")}
        ${listCard("Protein families behind the most predictions", "top-families", "#/browse/families")}
      </div>
    </section>`;

  // player cards
  const players = el.querySelector("#players");
  for (const s of model) {
    const card = document.createElement("div");
    card.className = "card player";
    card.innerHTML = `<div class="radar"></div>
      <div>
        <div class="name" style="color:${s.color}">${esc(s.binomial)}</div>
        <div class="common">${esc(s.common)}${s.ontology ? " · " + esc(s.ontology) : ""}</div>
        <div class="kv">
          <span class="k">phenologs</span><span class="v">${fmt.int(s.n_phenologs)}</span>
          <span class="k">diseases</span><span class="v">${fmt.int(s.n_diseases)}</span>
          <span class="k">phenotypes</span><span class="v">${fmt.int(s.n_phenotypes_matched)}</span>
          <span class="k">candidates</span><span class="v">${fmt.int(s.n_candidates)}</span>
          <span class="k">shared families</span><span class="v">${fmt.int(s.n_common_with_base)}</span>
          <span class="k">p cutoff</span><span class="v mono">${fmt.p(s.fdr_cutoff)}</span>
        </div></div>`;
    card.addEventListener("click", () => { location.hash = route.species(s.species); });
    players.appendChild(card);
    radar(card.querySelector(".radar"), speciesAxes(s, model, RADAR_AXES), { size: 150, color: s.color });
  }

  // tree
  const treeEl = el.querySelector("#tree");
  const sel = el.querySelector("#tree-metric");
  const drawTree = () => speciesTree(treeEl, ov.tree, ov.species, { metric: sel.value, metricLabel: TREE_METRICS[sel.value], width: Math.max(620, treeEl.clientWidth - 4) });
  sel.addEventListener("change", drawTree);
  drawTree();

  // fdr table
  const fdrRows = ov.fdr.filter(f => f.fdr === r.fdr).sort((a, b) => b.phenolog_count - a.phenolog_count);
  el.querySelector("#fdr").innerHTML = `<table class="tbl"><thead><tr>
      <th>species</th><th class="num">p cutoff</th><th class="num">phenologs</th><th class="num">expected false</th><th class="num">realised fdr</th></tr></thead><tbody>
    ${fdrRows.map(f => `<tr><td><a href="${route.species(f.species)}"><span class="dot" style="background:${speciesOf(f.species).color}"></span> <i>${esc(speciesOf(f.species).binomial)}</i></a></td>
      <td class="num mono">${fmt.p(f.pvalue_cutoff)}</td><td class="num">${fmt.int(f.phenolog_count)}</td>
      <td class="num">${Number(f.expected_false_positives).toFixed(1)}</td><td class="num">${fmt.pct(f.realized_fdr)}</td></tr>`).join("")}
    </tbody></table>`;

  pvalStrip(el.querySelector("#pvals"), ov.pval_bins, model);

  // top lists
  const top = ov.top;
  barList(el.querySelector("#top-diseases"), top.diseases.map(d => ({ label: d.label, value: d.value, href: route.disease(d.id), sub: `${d.n_species} sp.`, color: speciesOf(d.top_species).color })));
  barList(el.querySelector("#top-phenotypes"), top.phenotypes.map(d => ({ label: d.label, value: d.value, href: route.phenotype(d.id), sub: `<i>${esc(speciesOf(d.species).binomial)}</i>`, color: speciesOf(d.species).color })));
  barList(el.querySelector("#top-genes"), top.genes.map(d => ({ label: d.label, value: d.value, href: route.gene(d.id), sub: `${d.n_known_diseases} known · ${esc(d.family || "")}` })));
  barList(el.querySelector("#top-families"), top.families.map(d => ({ label: d.label, value: d.value, href: route.family(d.id), sub: esc((d.base_symbols || "").split(", ").slice(0, 3).join(", ")) })));
}

function tile(n, label, cls = "") {
  return `<div class="tile ${cls}"><div class="n">${fmt.int(n)}</div><div class="l">${label}</div></div>`;
}

function listCard(title, id, browseHref, sub = "") {
  return `<div class="card"><div class="card-head"><h2>${title}</h2><span class="sub">${sub ? sub + " · " : ""}<a href="${browseHref}">browse all</a></span></div><div id="${id}"></div></div>`;
}

// One small bar strip per species: how strong its phenologs are
function pvalStrip(el, bins, species) {
  const width = Math.max(300, el.clientWidth || 500), rowH = 30, left = 160, right = 40;
  const maxX = 30;
  bins = d3.rollups(bins, v => d3.sum(v, b => b.n), b => b.species, b => Math.min(b.neglog10, maxX))
           .flatMap(([species, byBin]) => byBin.map(([neglog10, n]) => ({ species, neglog10, n })));
  const x = d3.scaleLinear().domain([2, maxX + 1]).range([0, width - left - right]);
  const svg = d3.select(el).append("svg").attr("width", width).attr("height", species.length * rowH + 28);
  species.forEach((s, i) => {
    const mine = bins.filter(b => b.species === s.species);
    const total = d3.sum(mine, b => b.n) || 1;
    const maxN = d3.max(mine, b => b.n) || 1;
    const y = d3.scaleLinear().domain([0, maxN]).range([0, rowH - 8]);
    const g = svg.append("g").attr("transform", `translate(${left},${i * rowH})`);
    g.append("text").attr("x", -8).attr("y", rowH / 2 + 4).attr("text-anchor", "end").attr("font-size", 11.5).attr("font-style", "italic").attr("fill", "#e6ebf2").text(s.binomial);
    g.selectAll("rect").data(mine).join("rect")
      .attr("x", b => x(b.neglog10)).attr("width", Math.max(1, x(3) - x(2) - 1))
      .attr("y", b => rowH - 4 - y(b.n)).attr("height", b => y(b.n))
      .attr("fill", s.color).attr("fill-opacity", 0.85)
      .on("mousemove", (e, b) => { showTip(e, `<b>${esc(s.binomial)}</b><br>−log10 p in [${b.neglog10}, ${b.neglog10 + 1}): ${fmt.int(b.n)} phenologs (${fmt.pct(b.n / total)})`); moveTip(e); })
      .on("mouseleave", hideTip);
  });
  svg.append("g").attr("class", "axis").attr("transform", `translate(${left},${species.length * rowH + 2})`)
    .call(d3.axisBottom(x).ticks(8).tickFormat(d => d >= maxX ? `${maxX}+` : d)).call(a => a.select(".domain").remove());
  svg.append("text").attr("x", left + x.range()[1] / 2).attr("y", species.length * rowH + 27).attr("text-anchor", "middle").attr("font-size", 10.5).attr("fill", "#5c687a").text("−log10 p (bar height ∝ phenologs in bin, per species)");
}
