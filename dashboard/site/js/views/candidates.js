// Every gene-disease candidate pair in one list across all diseases, best score first (the paper's global
// ranking). 1.5M rows stay in DuckDB; each page, sort and filter is a query
import { ensure, query, lit } from "../db.js";
import { state, fmt, esc, speciesOf, occurrenceChips, route } from "../format.js";
import { legend, PAIR_KEY } from "./common.js";
import { pagedTable } from "../charts/table.js";

const SORTABLE = ["global_rank", "symbol", "disease_name", "rank", "score", "family"];

export async function render(el) {
  await ensure("candidates", "genes", "diseases");
  const ov = state.overview;
  const model = ov.species.filter(s => !s.is_base);
  const filters = { known: false, species: null };

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Candidates</div><h1>Gene–disease candidates</h1>
      <div class="muted">Every gene–disease pair the ranking scored, in one list across all diseases, best ${esc(ov.release.rank_metric)} score first</div></div></div>
    <div class="card"><div class="card-head"><h2>Global ranking</h2><span class="sub">${legend(PAIR_KEY)}</span></div><div id="list"></div></div>`;

  const tools = `<span class="chip toggle off" id="known-toggle">include known pairs</span>
    <span class="chips" id="sp-filter">${model.map(s => `<span class="chip toggle" data-sp="${s.species}"><span class="dot" style="background:${s.color}"></span><i>${esc(s.binomial)}</i></span>`).join("")}</span>`;

  const columns = [
    { key: "global_rank", label: "#", num: true, sortDir: "asc", help: "Place among the novel pairs of every disease (ties share a place); known pairs are not ranked",
      format: v => v == null ? `<span class="known-mark">known</span>` : fmt.int(v) },
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
    { key: "disease_name", label: "Human disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}">${esc(v)}</a> <span class="id faint">${esc(r.disease_id)}</span>` },
    { key: "rank", label: "Rank in disease", num: true, sortDir: "asc" },
    { key: "score", label: "Score (p)", num: true, sortDir: "asc", format: v => `<span class="mono">${fmt.p(v)}</span>` },
    { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : "" },
    { key: "species_occurrence", label: "Supporting phenologs by species", sortable: false, format: v => occurrenceChips(v) },
  ];

  async function fetchPage({ sortKey, sortDir, filter, page, pageSize }) {
    const where = [];
    if (!filters.known) where.push("NOT known");
    if (filters.species) where.push(`species_occurrence LIKE ${lit("%" + filters.species + ":%")}`);
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
      query(`SELECT * FROM candidates ${w} ORDER BY ${key} ${sortDir === "desc" ? "DESC" : "ASC"} NULLS LAST, score, disease_id, rank
             LIMIT ${pageSize} OFFSET ${page * pageSize}`),
      query(`SELECT count(*)::INTEGER AS n FROM candidates ${w}`),
    ]);
    return { rows, total: n };
  }

  const t = pagedTable(el.querySelector("#list"), columns, fetchPage, {
    sortKey: "score", sortDir: "asc", pageSize: 50, search: true, tools,
    placeholder: "Filter by gene or disease name, or an exact id…",
    onDraw: root => {
      const tog = root.querySelector("#known-toggle");
      tog.classList.toggle("off", !filters.known);
      tog.addEventListener("click", () => { filters.known = !filters.known; t.reload(); });
      root.querySelectorAll("#sp-filter .chip").forEach(chip => {
        chip.classList.toggle("off", filters.species != null && filters.species !== chip.dataset.sp);
        chip.addEventListener("click", () => { filters.species = filters.species === chip.dataset.sp ? null : chip.dataset.sp; t.reload(); });
      });
    },
  });
}
