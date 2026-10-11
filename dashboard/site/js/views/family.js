// One PANTHER protein family: its members across species, and the diseases it is implicated in
import { ensure, query, lit } from "../db.js";
import { state, fmt, esc, speciesChip, speciesOf, route, links } from "../format.js";
import { pageHead, statRow, card, notFound } from "./common.js";
import { table } from "../charts/table.js";

export async function render(el, id) {
  await ensure("families", "genes", "candidates", "diseases", "known_associations");
  const [f] = await query(`SELECT * FROM families WHERE family = ${lit(id)}`);
  if (!f) { el.innerHTML = notFound("protein family", id); return; }
  const base = state.overview.release.base_species;

  const [members, predicted, known] = await Promise.all([
    query(`SELECT species, gene_id, symbol, n_phenotypes, n_known_diseases FROM genes WHERE family = ${lit(id)} ORDER BY species, symbol`),
    query(`SELECT c.disease_id, d.name, min(c.score) AS best_score, min(c.rank) AS best_rank, bool_or(c.known) AS known,
                  d.n_candidates, any_value(c.species_occurrence) AS species_occurrence
           FROM candidates c JOIN diseases d USING (disease_id) WHERE c.family = ${lit(id)}
           GROUP BY c.disease_id, d.name, d.n_candidates ORDER BY best_score`),
    query(`SELECT DISTINCT disease_id, disease_name, symbol FROM known_associations WHERE family = ${lit(id)} ORDER BY disease_name`),
  ]);
  const bySp = d3.group(members, m => m.species);
  const order = state.overview.species.map(s => s.species);
  const speciesList = [...bySp.keys()].sort((a, b) => order.indexOf(a) - order.indexOf(b));

  el.innerHTML = `
    ${pageHead({ kicker: "Protein family", title: `<span class="mono" style="font-size:inherit">${esc(f.family)}</span>`,
                 links: [{ label: "PANTHER", href: links.panther(f.family) }],
                 after: `<div class="chips" style="align-self:center;max-width:560px">${(f.base_symbols || "").split(", ").filter(Boolean).slice(0, 12).map(s => `<a class="chip" href="${route.gene(s)}">${esc(s)}</a>`).join("")}</div>` })}
    ${statRow([[f.n_genes, "genes"], [f.n_species, "species"], [f.n_base_genes, "human genes"], [f.n_known_diseases, "known diseases"],
               [f.n_predicted_diseases, "diseases ranked"], [f.n_novel_diseases, "of them novel"], [fmt.p(f.best_score), "best score"]])}
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${card("Members by species", "members", { sub: "every gene of the family in this release · amber = has a curated disease" })}
      ${card("Diseases this family is ranked for", "pred", { sub: "best score of any member gene · ★ = a member is already known for the disease" })}
      ${known.length ? card("Known associations", "known", { sub: "curated gene–disease pairs involving this family" }) : ""}
    </div>`;

  el.querySelector("#members").innerHTML = speciesList.map(sp => {
    const ms = bySp.get(sp);
    return `<div style="margin-bottom:8px">${speciesChip(sp, { count: ms.length })}
      <span class="chips" style="display:inline-flex;margin-left:6px">${ms.slice(0, 40).map(m => `<a class="chip ${m.n_known_diseases ? "known" : ""}" href="${route.gene(m.gene_id)}" title="${m.n_phenotypes} phenotypes">${esc(m.symbol)}</a>`).join("")}
      ${ms.length > 40 ? `<span class="faint">+${ms.length - 40} more</span>` : ""}</span></div>`;
  }).join("");

  table(el.querySelector("#pred"), predicted, [
    { key: "name", label: "Disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>${r.known ? " ★" : ""}` },
    { key: "best_rank", label: "Best rank", num: true, format: (v, r) => `${v} <span class="faint">/ ${fmt.int(r.n_candidates)}</span>` },
    { key: "best_score", label: "Score (p)", num: true, format: v => `<span class="mono">${fmt.p(v)}</span>` },
  ], { sortKey: "best_score", sortDir: "asc", search: ["name", "disease_id"], empty: "Not in any candidate ranking" });

  if (known.length) {
    el.querySelector("#known").innerHTML = `<div class="chips">${known.map(k => `<a class="chip" href="${route.disease(k.disease_id)}"><b>${esc(k.symbol)}</b>&nbsp;· ${esc(k.disease_name)}</a>`).join("")}</div>`;
  }
}
