// Full search results, grouped by kind, in left-aligned columns. Genes and phenotypes show their species,
// since the results span every species
import { ensure, query, lit } from "./db.js";
import { esc, speciesOf, shortBinomial } from "./format.js";

const ORDER = ["disease", "gene", "phenotype", "family", "species"];
const SHOWN = 80;

export async function render(el, q) {
  q = (q || "").trim();
  await ensure("search");
  const like = lit("%" + q + "%"), prefix = lit(q + "%");
  const rows = q.length < 2 ? [] : await query(`
    SELECT kind, id, label, species, detail FROM search
    WHERE label ILIKE ${like} OR id ILIKE ${like}
    ORDER BY (id ILIKE ${prefix}) DESC, (label ILIKE ${prefix}) DESC, weight DESC, label
    LIMIT 400`);
  const groups = d3.group(rows, r => r.kind);

  const item = (k, r) => {
    const sp = r.species && k !== "species" ? speciesOf(r.species) : null;
    return `<a href="#/${k}/${encodeURIComponent(r.id)}" title="${esc(r.detail || "")}">` +
      (sp ? `<span class="dot" style="background:${sp.color}"></span> ${esc(r.label)} <span class="faint"><i>${esc(shortBinomial(r.species))}</i></span>`
          : `${esc(r.label)} <span class="faint id">${esc(r.id)}</span>`) + `</a>`;
  };

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Search</div><h1>“${esc(q)}”</h1><div class="muted">${rows.length ? `${rows.length}${rows.length === 400 ? "+" : ""} matches` : "no matches"}</div></div></div>
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${ORDER.filter(k => groups.has(k)).map(k => `
        <div class="card"><div class="card-head"><h2>${k}${k === "species" ? "" : "s"}</h2><span class="sub">${groups.get(k).length}</span></div>
          <div class="name-grid wide">${groups.get(k).slice(0, SHOWN).map(r => item(k, r)).join("")}</div>
          ${groups.get(k).length > SHOWN ? `<div class="faint" style="margin-top:6px">+${groups.get(k).length - SHOWN} more</div>` : ""}</div>`).join("")}
    </div>`;
}
