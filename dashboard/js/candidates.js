// Every gene-disease pair the candidate ranking scored, known and novel, in one list across all diseases,
// best score first. 1.5M rows stay in DuckDB; each page, sort and filter is a query
import { ensure, query, lit } from "./db.js";
import { state, fmt, esc, supportChips, route } from "./format.js";
import { legend, PAIR_KEY, HELP, SUPPORT_LABEL } from "./common.js";
import { pagedTable } from "./table.js";

const SORTABLE = ["global_rank", "symbol", "disease_name", "rank", "score", "family", "n_phenologs"];

export async function render(el) {
  await ensure("candidates", "genes", "diseases");
  const ov = state.overview;
  const model = ov.species.filter(s => !s.is_base);
  const filters = { novelOnly: false, species: null };

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Candidates</div><h1>Gene–disease candidates</h1>
      <div class="muted">Every gene–disease pair the candidate ranking scored, across all diseases, best score first.
        Diseases with many phenologs get smaller scores, so they fill the top of the list</div></div></div>
    <div class="card"><div class="card-head"><h2>Global ranking</h2><span class="sub">${legend(PAIR_KEY)}</span></div><div id="list"></div></div>`;

  const tools = `<span class="chip toggle off" id="novel-toggle">novel only</span>
    <span class="chips" id="sp-filter"><span class="faint">supported by</span> ${model.map(s => `<span class="chip toggle" data-sp="${s.species}"><span class="dot" style="background:${s.color}"></span><i>${esc(s.binomial)}</i></span>`).join("")}</span>`;

  const columns = [
    { key: "global_rank", label: "Global rank", num: true, sortDir: "asc", help: "Place among every scored gene–disease pair of every disease, by score (ties share a place)" },
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
    { key: "disease_name", label: "Human disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}">${esc(v)}</a> <span class="id faint">${esc(r.disease_id)}</span>` },
    { key: "rank", label: "Rank in disease", num: true, sortDir: "asc", help: HELP.rank },
    { key: "score", label: "Score", num: true, sortDir: "asc", help: HELP.score, format: v => `<span class="mono">${fmt.p(v)}</span>` },
    { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : "" },
    { key: "n_phenologs", label: SUPPORT_LABEL, sortDir: "desc", help: HELP.support, format: (v, r) => supportChips(r.support, r.symbol) },
  ];

  async function fetchPage({ sortKey, sortDir, filter, page, pageSize }) {
    const where = [];
    if (filters.novelOnly) where.push("NOT known");
    if (filters.species) where.push(`support LIKE ${lit("%" + filters.species + ":%")}`);
    const q = filter.trim();
    if (q) {
      // names are matched in the small gene and disease tables, then the 1.5M pairs are filtered by id
      const like = lit("%" + q + "%");
      where.push(`(gene_id IN (SELECT gene_id FROM genes WHERE species = ${lit(ov.release.base_species)} AND (symbol ILIKE ${like} OR gene_id = ${lit(q)}))
                   OR disease_id IN (SELECT disease_id FROM diseases WHERE name ILIKE ${like} OR disease_id = ${lit(q)}))`);
    }
    const w = where.length ? "WHERE " + where.join(" AND ") : "";
    const key = SORTABLE.includes(sortKey) ? sortKey : "score";
    const [rows, [{ n }]] = await Promise.all([
      query(`SELECT * FROM candidates ${w} ORDER BY ${key} ${sortDir === "desc" ? "DESC" : "ASC"} NULLS LAST, score, disease_id, rank, symbol
             LIMIT ${pageSize} OFFSET ${page * pageSize}`),
      query(`SELECT count(*)::INTEGER AS n FROM candidates ${w}`),
    ]);
    return { rows, total: n };
  }

  const t = pagedTable(el.querySelector("#list"), columns, fetchPage, {
    sortKey: "score", sortDir: "asc", pageSize: 50, search: true, tools,
    placeholder: "Filter by gene or disease name, or an exact id…",
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
