// Every gene-term pair the release keeps, known and novel, in one list across all diseases or human phenotypes,
// best score first. Up to 1.5M rows stay in DuckDB; each page, sort and filter is a query
import { ensure, query, lit } from "./db.js";
import { state, fmt, esc, supportChips, route } from "./format.js";
import { legend, pairKey, supportLabel, HELP } from "./common.js";
import { pagedTable, EXPORT_MAX } from "./table.js";

const SORTABLE = ["global_rank", "symbol", "term_name", "rank", "score", "family", "n_phenologs"];

export async function render(el) {
  await ensure("candidates", "genes", "terms");
  const ov = state.overview, net = state.net, top = ov.release.top_candidates;
  const model = ov.species.filter(s => !s.is_base);
  const filters = { novelOnly: false, species: null };
  const Pair = net.pair[0].toUpperCase() + net.pair.slice(1);

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Candidates</div><h1>${Pair} candidates</h1>
      <div class="muted">${top ? `Each ${net.term}'s best ${fmt.int(top)} candidate genes` : `Every ${net.pair} pair the candidate ranking scored`},
        across all ${net.terms}, best score first. ${net.Terms} with many phenologs get smaller scores, so they fill the top of the list</div></div></div>
    <div class="card"><div class="card-head"><h2>Global ranking</h2><span class="sub">${legend(pairKey())}</span></div><div id="list"></div></div>`;

  const tools = `<span class="chip toggle off" id="novel-toggle">novel only</span>
    <span class="chips" id="sp-filter"><span class="faint">supported by</span> ${model.map(s => `<span class="chip toggle" data-sp="${s.species}"><span class="dot" style="background:${s.color}"></span><i>${esc(s.binomial)}</i></span>`).join("")}</span>`;

  const columns = [
    { key: "global_rank", label: "Global rank", num: true, sortDir: "asc",
      help: `Place among every ${top ? "kept" : "scored"} ${net.pair} pair of every ${net.term}, by score (ties share a place)` },
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
    { key: "term_name", label: net.human, format: (v, r) => `<a href="${route.term(r.term_id)}">${esc(v)}</a> <span class="id faint">${esc(r.term_id)}</span>` },
    { key: "rank", label: `Rank in ${net.term}`, num: true, sortDir: "asc", help: HELP.rank },
    { key: "score", label: "Score", num: true, sortDir: "asc", help: HELP.score, format: v => `<span class="mono">${fmt.p(v)}</span>` },
    { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : "" },
    { key: "n_phenologs", label: supportLabel(), sortDir: "desc", help: HELP.support, format: (v, r) => supportChips(r.support, r.symbol) },
  ];

  function where(filter) {
    const w = [];
    if (filters.novelOnly) w.push("NOT known");
    if (filters.species) w.push(`support LIKE ${lit("%" + filters.species + ":%")}`);
    const q = filter.trim();
    if (q) {
      // names are matched in the small gene and term tables, then the pairs are filtered by id
      const like = lit("%" + q + "%");
      w.push(`(gene_id IN (SELECT gene_id FROM genes WHERE species = ${lit(ov.release.base_species)} AND (symbol ILIKE ${like} OR gene_id = ${lit(q)}))
               OR term_id IN (SELECT term_id FROM terms WHERE name ILIKE ${like} OR term_id = ${lit(q)}))`);
    }
    return w.length ? "WHERE " + w.join(" AND ") : "";
  }
  const order = (sortKey, sortDir) => `ORDER BY ${SORTABLE.includes(sortKey) ? sortKey : "score"} ${sortDir === "desc" ? "DESC" : "ASC"} NULLS LAST, score, term_id, rank, symbol`;

  async function fetchPage({ sortKey, sortDir, filter, page, pageSize }) {
    const w = where(filter);
    const [rows, [{ n }]] = await Promise.all([
      query(`SELECT * FROM candidates ${w} ${order(sortKey, sortDir)} LIMIT ${pageSize} OFFSET ${page * pageSize}`),
      query(`SELECT count(*)::INTEGER AS n FROM candidates ${w}`),
    ]);
    return { rows, total: n };
  }

  const t = pagedTable(el.querySelector("#list"), columns, fetchPage, {
    sortKey: "score", sortDir: "asc", pageSize: 50, search: true, tools,
    placeholder: `Filter by gene or ${net.term} name, or an exact id…`,
    download: `${net.browse}_candidates`,
    fetchAll: ({ sortKey, sortDir, filter }) => query(`SELECT * FROM candidates ${where(filter)} ${order(sortKey, sortDir)} LIMIT ${EXPORT_MAX + 1}`),
    onDraw: root => {
      const tog = root.querySelector("#novel-toggle");
      tog.classList.toggle("off", !filters.novelOnly);
      tog.addEventListener("click", () => { filters.novelOnly = !filters.novelOnly; t.reload(); });
      root.querySelectorAll("#sp-filter .chip").forEach(chip => {
        chip.classList.toggle("off", filters.species != null && filters.species !== chip.dataset.sp);
        chip.addEventListener("click", () => { filters.species = filters.species === chip.dataset.sp ? null : chip.dataset.sp; t.reload(); });
      });
    },
  });
}
