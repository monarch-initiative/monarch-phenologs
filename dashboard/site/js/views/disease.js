// One disease: its known genes, its phenologs, the additive decomposition and the candidate ranking
import { ensure, query, lit } from "../db.js";
import { state, fmt, esc, speciesOf, speciesChip, occurrenceChips, route, links } from "../format.js";
import { pageHead, statRow, card, notFound } from "./common.js";
import { table } from "../charts/table.js";
import { additiveChart, ORDERS } from "../charts/stacked.js";

export async function render(el, id) {
  await ensure("diseases", "phenologs", "candidates", "additive", "known_associations");
  const [d] = await query(`SELECT * FROM diseases WHERE disease_id = ${lit(id)}`);
  if (!d) { el.innerHTML = notFound("disease with phenologs", id); return; }

  const [known, phen, cands, add] = await Promise.all([
    query(`SELECT gene_id, symbol, predicate, family FROM known_associations WHERE disease_id = ${lit(id)} ORDER BY symbol`),
    query(`SELECT * FROM phenologs WHERE disease_id = ${lit(id)} ORDER BY rank`),
    query(`SELECT rank, gene_id, symbol, family, score, known, species_occurrence, n_species FROM candidates WHERE disease_id = ${lit(id)} ORDER BY rank`),
    query(`SELECT * FROM additive WHERE disease_id = ${lit(id)} ORDER BY family_rank, phenolog_rank`),
  ]);
  const recovered = new Map(cands.filter(c => c.known).map(c => [c.gene_id, c.rank]));
  const bySpecies = JSON.parse(d.by_species || "{}");
  const dismech = links.dismech(d.disease_id, d.name);

  el.innerHTML = `
    ${pageHead({ kicker: "Disease", title: esc(d.name), id: d.disease_id, idHref: links.monarch(d.disease_id),
                 links: [{ label: "Monarch", href: links.monarch(d.disease_id) }, dismech ? { label: "DisMech", href: dismech } : null],
                 after: `<div class="chips" style="align-self:center">${Object.entries(bySpecies).map(([s, n]) => speciesChip(s, { count: n })).join("")}</div>` })}
    ${statRow([[d.n_known_genes, "known genes"], [d.n_phenologs, "phenologs"], [d.n_species, "species"],
               [fmt.p(d.best_pval), "best phenolog p"], [d.n_candidates, "ranked genes"],
               [`${fmt.int(d.n_known_recovered)} / ${fmt.int(d.n_known_genes)}`, "known genes recovered"]])}

    <div class="grid" style="grid-template-columns: 1fr; gap:16px">
      ${card("Known genes", "known", { sub: "curated gene–disease associations in the graph · amber = recovered by the ranking, with its rank" })}
      ${card("Additive decomposition", "additive", { tools: `<div class="controls" id="add-controls"></div>` })}
      ${card("Candidate genes", "cands", { sub: `ranked by the ${state.overview.release.rank_metric} score over the k=${state.overview.release.kneighbs} nearest phenologs` })}
      ${card("Phenologs", "phens", { sub: "model phenotypes whose protein families overlap this disease's more than chance" })}
    </div>`;

  // known genes
  el.querySelector("#known").innerHTML = known.length
    ? `<div class="chips">${known.map(g => `<a class="chip ${recovered.has(g.gene_id) ? "known" : ""}" href="${route.gene(g.gene_id)}">${esc(g.symbol)}${recovered.has(g.gene_id) ? ` <span class="muted">#${recovered.get(g.gene_id)}</span>` : ""}</a>`).join("")}</div>`
    : `<div class="empty">No curated genes for this disease in the graph</div>`;

  // additive decomposition with its controls
  additiveSection(el.querySelector("#additive"), el.querySelector("#add-controls"), add);

  // candidates
  let novelOnly = false;
  const candCols = [
    { key: "rank", label: "#", num: true, width: "50px" },
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>${r.known ? ' <span class="faint" title="known association">★</span>' : ""}` },
    { key: "family", label: "Family", format: v => `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` },
    { key: "score", label: "Score (p)", num: true, format: v => `<span class="mono">${fmt.p(v)}</span>` },
    { key: "n_species", label: "Species", num: true },
    { key: "species_occurrence", label: "Supporting phenologs by species", format: v => occurrenceChips(v) },
  ];
  const candTable = table(el.querySelector("#cands"), cands, candCols, {
    sortKey: "rank", sortDir: "asc", search: ["symbol", "family"], pageSize: 25,
    tools: `<label class="chip toggle" id="novel-toggle"><input type="checkbox" style="display:none">novel only</label>`,
    onDraw: root => {
      const tog = root.querySelector("#novel-toggle");
      tog.classList.toggle("off", !novelOnly);
      tog.addEventListener("click", () => { novelOnly = !novelOnly; candTable.setRows(novelOnly ? cands.filter(c => !c.known) : cands); });
    },
  });

  // phenologs
  table(el.querySelector("#phens"), phen, [
    { key: "rank", label: "#", num: true, width: "50px" },
    { key: "species", label: "Species", format: v => speciesChip(v) },
    { key: "phenotype_name", label: "Model phenotype", format: (v, r) => `<a href="${route.phenotype(r.phenotype_id)}">${esc(v)}</a> <span class="id faint">${esc(r.phenotype_id)}</span>` },
    { key: "overlap", label: "Overlap", num: true },
    { key: "n_disease_orthogroups", label: "Disease fam.", num: true },
    { key: "n_phenotype_orthogroups", label: "Phenotype fam.", num: true },
    { key: "n_common", label: "Common", num: true },
    { key: "pval", label: "p", num: true, format: v => `<span class="mono">${fmt.p(v)}</span>` },
  ], { sortKey: "rank", sortDir: "asc", search: ["phenotype_name", "species", "phenotype_id"] });
}

// Controls: how many phenologs (k), how many families, which species, and the ordering
function additiveSection(chartEl, controlsEl, rows) {
  if (!rows.length) { chartEl.innerHTML = `<div class="empty">No decomposition available</div>`; return; }
  const maxK = Math.max(...rows.map(r => r.phenolog_rank));
  const speciesPresent = [...new Set(rows.map(r => r.species))];
  const opts = { k: Math.min(20, maxK), topFamilies: 20, order: "total", species: new Set(speciesPresent) };

  controlsEl.innerHTML = `
    <label>phenologs <input type="range" id="k" min="1" max="${maxK}" value="${opts.k}"> <span class="val" id="k-val">${opts.k}</span></label>
    <label>families <select id="fams">${[5, 10, 20, 30, 40].map(n => `<option ${n === 20 ? "selected" : ""}>${n}</option>`).join("")}</select></label>
    <label>order <select id="order">${Object.entries(ORDERS).map(([k, [label]]) => `<option value="${k}">${label}</option>`).join("")}</select></label>
    <span class="chips" id="sp-toggles">${speciesPresent.map(s => `<span class="chip toggle" data-sp="${s}"><span class="dot" style="background:${speciesOf(s).color}"></span><i>${esc(speciesOf(s).binomial)}</i></span>`).join("")}</span>`;

  const draw = () => additiveChart(chartEl, rows, opts);
  controlsEl.querySelector("#k").addEventListener("input", e => { opts.k = +e.target.value; controlsEl.querySelector("#k-val").textContent = opts.k; draw(); });
  controlsEl.querySelector("#fams").addEventListener("change", e => { opts.topFamilies = +e.target.value; draw(); });
  controlsEl.querySelector("#order").addEventListener("change", e => { opts.order = e.target.value; draw(); });
  controlsEl.querySelectorAll("#sp-toggles .chip").forEach(chip => chip.addEventListener("click", () => {
    const s = chip.dataset.sp;
    if (opts.species.has(s)) opts.species.delete(s); else opts.species.add(s);
    chip.classList.toggle("off", !opts.species.has(s));
    draw();
  }));
  draw();
  chartEl.insertAdjacentHTML("beforeend", `<div class="note">Each bar is a PANTHER family; each segment one of the disease's top phenologs
    that the family carries, sized (overlap / disease families) × (1 − p) as in Woods et al. 2013. Human genes of the family are listed at the bar's end,
    <span class="known-mark">known disease genes</span> first.</div>`);
}
