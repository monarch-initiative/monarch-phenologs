// Full search results, grouped by kind
import { ensure, query, lit } from "../db.js";
import { esc } from "../format.js";

const ORDER = ["disease", "gene", "phenotype", "family", "species"];

export async function render(el, q) {
  q = (q || "").trim();
  await ensure("search");
  const like = lit("%" + q + "%"), prefix = lit(q + "%");
  const rows = q.length < 2 ? [] : await query(`
    SELECT kind, id, label, detail FROM search
    WHERE label ILIKE ${like} OR id ILIKE ${like}
    ORDER BY (id ILIKE ${prefix}) DESC, (label ILIKE ${prefix}) DESC, weight DESC, label
    LIMIT 400`);
  const groups = d3.group(rows, r => r.kind);

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Search</div><h1>“${esc(q)}”</h1><div class="muted">${rows.length ? `${rows.length}${rows.length === 400 ? "+" : ""} matches` : "no matches"}</div></div></div>
    <div class="grid cols-2">
      ${ORDER.filter(k => groups.has(k)).map(k => `
        <div class="card"><div class="card-head"><h2>${k}${k === "species" ? "" : "s"}</h2><span class="sub">${groups.get(k).length}</span></div>
          <div class="chips">${groups.get(k).slice(0, 80).map(r => `<a class="chip" href="#/${k}/${encodeURIComponent(r.id)}" title="${esc(r.detail || "")}">${esc(r.label)} <span class="faint id">${esc(r.id)}</span></a>`).join("")}
          ${groups.get(k).length > 80 ? `<span class="faint">+${groups.get(k).length - 80} more</span>` : ""}</div></div>`).join("")}
    </div>`;
}
