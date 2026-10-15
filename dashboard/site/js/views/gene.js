// One gene: what it is known for, what it is predicted for, its phenotypes, its protein family
import { ensure, query, lit } from "../db.js";
import { state, fmt, esc, speciesOf, speciesChip, occurrenceChips, route, links } from "../format.js";
import { pageHead, statRow, card, notFound, legend, PAIR_KEY } from "./common.js";
import { table } from "../charts/table.js";

export async function render(el, id) {
  await ensure("genes", "known_associations", "candidates", "gene_phenotypes", "diseases");
  let rows = await query(`SELECT * FROM genes WHERE gene_id = ${lit(id)}`);
  if (!rows.length) {
    // allow a symbol in the url, base species first
    rows = await query(`SELECT * FROM genes WHERE upper(symbol) = upper(${lit(id)}) ORDER BY (species = ${lit(state.overview.release.base_species)}) DESC LIMIT 1`);
    if (rows.length) { location.replace(route.gene(rows[0].gene_id)); return; }
    el.innerHTML = notFound("gene", id); return;
  }
  const g = rows[0];
  const isBase = g.species === state.overview.release.base_species;

  const [known, predicted, phenos, members] = await Promise.all([
    query(`SELECT k.disease_id, k.disease_name, k.predicate, d.n_phenologs FROM known_associations k LEFT JOIN diseases d USING (disease_id)
           WHERE k.gene_id = ${lit(g.gene_id)} ORDER BY k.disease_name`),
    query(`SELECT c.rank, c.disease_id, d.name, c.score, c.known, c.species_occurrence, d.n_candidates
           FROM candidates c JOIN diseases d USING (disease_id) WHERE c.gene_id = ${lit(g.gene_id)} ORDER BY c.score`),
    query(`SELECT phenotype_id, phenotype_name FROM gene_phenotypes WHERE gene_id = ${lit(g.gene_id)} ORDER BY phenotype_name`),
    g.family ? query(`SELECT species, gene_id, symbol FROM genes WHERE family = ${lit(g.family)} ORDER BY species, symbol`) : [],
  ]);

  el.innerHTML = `
    ${pageHead({ kicker: "Gene", title: `${esc(g.symbol)} <span class="muted" style="font-weight:400;font-size:18px">· <i>${esc(speciesOf(g.species).binomial)}</i></span>`,
                 id: g.gene_id, idHref: links.monarch(g.gene_id),
                 links: [{ label: "Monarch", href: links.monarch(g.gene_id) }, g.family ? { label: "PANTHER " + g.family, href: links.panther(g.family) } : null],
                 after: `<div style="align-self:center">${speciesChip(g.species)} ${g.family ? `<a class="chip mono" href="${route.family(g.family)}">${esc(g.family)}</a>` : ""}</div>` })}
    ${statRow([[g.n_phenotypes, "phenotypes"], [g.n_known_diseases, "known diseases"], [g.n_predicted_diseases, "diseases ranked"],
               [g.n_novel_diseases, "of them novel"], [g.best_rank ?? "–", "best rank"]])}
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${isBase ? card("Known diseases", "known", { sub: "curated associations in the graph" }) : ""}
      ${isBase ? card("Predicted diseases", "pred", { sub: `${legend(PAIR_KEY)}<span class="faint">every disease whose candidate ranking includes this gene</span>` }) : ""}
      ${card("Phenotypes", "phenos", { sub: isBase ? "HPO annotations of this gene" : "annotations that feed the phenolog calculation" })}
      ${g.family ? card("Protein family " + g.family, "members", { sub: "all genes of the family in this release" }) : ""}
    </div>`;

  if (isBase) {
    el.querySelector("#known").innerHTML = known.length
      ? `<div class="chips">${known.map(k => `<a class="chip" href="${route.disease(k.disease_id)}">${esc(k.disease_name)}${k.n_phenologs ? ` <span class="muted">${fmt.int(k.n_phenologs)} ph.</span>` : ""}</a>`).join("")}</div>`
      : `<div class="empty">No curated disease associations</div>`;

    table(el.querySelector("#pred"), predicted, [
      { key: "rank", label: "Rank", num: true, format: (v, r) => `${v} <span class="faint">/ ${fmt.int(r.n_candidates)}</span>` },
      { key: "name", label: "Disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
      { key: "score", label: "Score (p)", num: true, format: v => `<span class="mono">${fmt.p(v)}</span>` },
      { key: "species_occurrence", label: "Supporting phenologs by species", format: v => occurrenceChips(v) },
    ], { sortKey: "score", sortDir: "asc", search: ["name", "disease_id"], empty: "This gene is in no candidate ranking" });
  }

  table(el.querySelector("#phenos"), phenos, [
    { key: "phenotype_name", label: "Phenotype", format: (v, r) => `<a href="${route.phenotype(r.phenotype_id)}">${esc(v)}</a>` },
    { key: "phenotype_id", label: "Id", format: v => `<span class="id muted">${esc(v)}</span>` },
  ], { sortKey: "phenotype_name", sortDir: "asc", search: ["phenotype_name", "phenotype_id"], pageSize: 15, empty: "No phenotype annotations" });

  if (g.family) {
    const bySp = d3.group(members, m => m.species);
    el.querySelector("#members").innerHTML = [...bySp.entries()].map(([sp, ms]) => `
      <div style="margin-bottom:8px">${speciesChip(sp, { count: ms.length })}
        <span class="chips" style="display:inline-flex;margin-left:6px">${ms.slice(0, 30).map(m => `<a class="chip ${m.gene_id === g.gene_id ? "known" : ""}" href="${route.gene(m.gene_id)}">${esc(m.symbol)}</a>`).join("")}
        ${ms.length > 30 ? `<span class="faint">+${ms.length - 30} more</span>` : ""}</span></div>`).join("");
  }
}
