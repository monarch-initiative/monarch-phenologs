// One disease or human phenotype: its known genes, its phenologs, the additive decomposition and the candidate ranking
import { ensure, query, lit } from "./db.js";
import { state, fmt, esc, speciesOf, speciesChip, supportChips, route, links } from "./format.js";
import { pageHead, statRow, card, notFound, legend, nameGrid, pairKey, recoveryKey, supportLabel, HELP, PHENOLOG_HELP } from "./common.js";
import { table } from "./table.js";
import { additiveChart, ORDERS } from "./stacked.js";

const KNOWN_SHOWN = 120;     // known genes in the grid before "show all"

export async function render(el, id) {
  const net = state.net, r = state.overview.release;
  await ensure("terms", "phenologs", "candidates", "additive", "known_associations", "genes");
  const [d] = await query(`SELECT * FROM terms WHERE term_id = ${lit(id)}`);
  if (!d) { el.innerHTML = notFound(`${net.term} with phenologs`, id); return; }

  // a known gene's rank is its place in the whole ranking, also where the release keeps only the top of it
  const [known, phen, cands, add] = await Promise.all([
    query(`SELECT k.gene_id, any_value(k.symbol) AS symbol, min(k.rank) AS rank, coalesce(any_value(g.in_reach), false) AS in_reach
           FROM known_associations k LEFT JOIN genes g ON g.gene_id = k.gene_id
           WHERE k.term_id = ${lit(id)} GROUP BY k.gene_id ORDER BY symbol`),
    query(`SELECT * FROM phenologs WHERE term_id = ${lit(id)} ORDER BY rank`),
    query(`SELECT rank, gene_id, symbol, family, score, known, support, n_phenologs FROM candidates WHERE term_id = ${lit(id)} ORDER BY rank, symbol`),
    query(`SELECT * FROM additive WHERE term_id = ${lit(id)} ORDER BY family_rank, phenolog_rank`),
  ]);
  const bySpecies = JSON.parse(d.by_species || "{}");
  const dismech = net.key === "disease" ? links.dismech(d.term_id, d.name) : null;
  const kn = r.kneighbs;
  const over = d.n_phenologs > kn ? `the ${net.term}'s ${fmt.int(kn)} most significant phenologs` : `all of the ${net.term}'s phenologs`;
  const kept = cands.length < d.n_candidates ? ` · the best ${fmt.int(cands.length)} of the ${fmt.int(d.n_candidates)} genes it ranked` : "";

  el.innerHTML = `
    ${pageHead({ kicker: net.Term, title: esc(d.name), id: d.term_id, idHref: links.monarch(d.term_id),
                 links: [{ label: "Monarch", href: links.monarch(d.term_id) }, dismech ? { label: "DisMech", href: dismech } : null],
                 after: `<div style="align-self:center"><div class="faint" style="font-size:11.5px;margin-bottom:4px">phenologs by species</div>
                   <div class="chips">${Object.entries(bySpecies).map(([s, n]) => speciesChip(s, { count: n })).join("")}</div></div>` })}
    ${statRow([[d.n_known_genes, "known genes"], [d.n_phenologs, "phenologs"], [d.n_species, "species with phenologs"],
               [fmt.p(d.best_pval), "best phenolog p"], [d.n_candidates, "genes ranked"],
               [`${fmt.int(d.n_known_recovered)} / ${fmt.int(d.n_known_genes)}`, "known genes recovered"]])}

    <div class="grid" style="grid-template-columns: 1fr; gap:16px">
      ${card("Known genes", "known", { sub: legend([...recoveryKey(), ["var(--faint)", "gene has no ortholog with phenotypes"]]) })}
      ${card("Additive decomposition", "additive", { tools: `<div class="controls" id="add-controls"></div>` })}
      ${card("Candidate genes", "cands", { sub: `${legend(pairKey())}<span class="faint">ranked by score over ${over}${kept}</span>` })}
      ${card("Phenologs", "phens", { sub: `model phenotypes whose protein families overlap this ${net.term}'s more than chance` })}
    </div>`;

  // known genes; a broad human phenotype has thousands, so the grid starts with the first of them
  const knownEl = el.querySelector("#known");
  const knownGrid = all => {
    const shown = all ? known : known.slice(0, KNOWN_SHOWN);
    knownEl.innerHTML = nameGrid(shown.map(g => ({ label: g.symbol, href: route.gene(g.gene_id), cls: g.rank != null ? "known-mark" : g.in_reach ? "plain" : "out",
                                                   note: g.rank != null ? `#${fmt.int(g.rank)}` : "" }))) +
      (shown.length < known.length ? `<a class="more" href="">show all ${fmt.int(known.length)}</a>` : "");
    knownEl.querySelector("a.more")?.addEventListener("click", e => { e.preventDefault(); knownGrid(true); });
  };
  if (known.length) knownGrid(false);
  else knownEl.innerHTML = `<div class="empty">No curated genes for this ${net.term} in the graph</div>`;

  // additive decomposition with its controls; the known genes in it are marked as in the grid above
  const recovery = new Map(known.map(g => [g.symbol, g.rank != null ? "recovered" : g.in_reach ? "plain" : "out"]));
  additiveSection(el.querySelector("#additive"), el.querySelector("#add-controls"), add, recovery);

  // candidates
  let novelOnly = false;
  const candCols = [
    { key: "rank", label: "Rank", num: true, width: "60px", sortDir: "asc", help: HELP.rank },
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
    { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : "" },
    { key: "score", label: "Score", num: true, sortDir: "asc", help: HELP.score, format: v => `<span class="mono">${fmt.p(v)}</span>` },
    { key: "n_phenologs", label: supportLabel(), sortDir: "desc", help: HELP.support, format: (v, r) => supportChips(r.support, r.symbol) },
  ];
  const candTable = table(el.querySelector("#cands"), cands, candCols, {
    sortKey: "rank", sortDir: "asc", search: ["symbol", "family"], placeholder: "Filter by gene or family…", pageSize: 25,
    download: `${d.term_id}_candidates`,
    tools: `<label class="chip toggle" id="novel-toggle"><input type="checkbox" style="display:none">novel only</label>`,
    onDraw: root => {
      const tog = root.querySelector("#novel-toggle");
      tog.classList.toggle("off", !novelOnly);
      tog.addEventListener("click", () => { novelOnly = !novelOnly; candTable.setRows(novelOnly ? cands.filter(c => !c.known) : cands); });
    },
  });

  // phenologs
  table(el.querySelector("#phens"), phen, [
    { key: "rank", label: "Rank", num: true, width: "74px", sortDir: "asc", help: PHENOLOG_HELP.rank },
    { key: "species", label: "Species", format: v => speciesChip(v) },
    { key: "phenotype_name", label: "Model phenotype", format: (v, r) => `<a href="${route.phenotype(r.phenotype_id)}">${esc(v)}</a> <span class="id faint">${esc(r.phenotype_id)}</span>` },
    { key: "overlap", label: "Shared families", num: true, help: PHENOLOG_HELP.overlap },
    { key: "n_term_orthogroups", label: `${net.Term} families`, num: true, help: PHENOLOG_HELP.n_term_orthogroups },
    { key: "n_phenotype_orthogroups", label: "Phenotype families", num: true, help: PHENOLOG_HELP.n_phenotype_orthogroups },
    { key: "n_common", label: "Families in both species", num: true, help: PHENOLOG_HELP.n_common },
    { key: "pval", label: "p", num: true, sortDir: "asc", help: PHENOLOG_HELP.pval, format: v => `<span class="mono">${fmt.p(v)}</span>` },
  ], { sortKey: "rank", sortDir: "asc", search: ["phenotype_name", "species", "phenotype_id"], placeholder: "Filter by phenotype or species…", wrapHead: true,
       download: `${d.term_id}_phenologs` });
}

// Controls: how many phenologs (k), how many families, which species, and the ordering
function additiveSection(chartEl, controlsEl, rows, recovery) {
  if (!rows.length) { chartEl.innerHTML = `<div class="empty">No decomposition available</div>`; return; }
  const maxK = Math.max(...rows.map(r => r.phenolog_rank));
  const speciesPresent = [...new Set(rows.map(r => r.species))];
  const opts = { k: Math.min(20, maxK), topFamilies: 20, order: "total", species: new Set(speciesPresent), recovery };

  controlsEl.innerHTML = `
    <label>top phenologs <input type="range" id="k" min="1" max="${maxK}" value="${opts.k}"> <span class="val" id="k-val">${opts.k}</span></label>
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
  const t = state.net.term;
  chartEl.insertAdjacentHTML("beforeend", `<div class="note">Each bar is a PANTHER family and each segment one of the ${t}'s top phenologs
    whose phenotype a gene of the family carries in that species, sized (shared families / ${t} families) × (1 − p) (Woods et al. 2013).
    The candidate genes below are ranked by the hypergeometric score instead, which counts only each gene's own orthologs, so the two
    can differ: a family can lead here through genes other than its known one. Human genes of the family follow each bar, known genes
    of this ${t} first, in <span class="known-mark">amber</span> where the ranking recovered them.</div>`);
}
