// One phenotype: the human diseases it is a phenolog of, and the genes annotated with it
import { ensure, query, lit } from "./db.js";
import { state, fmt, esc, speciesChip, route, links } from "./format.js";
import { pageHead, statRow, card, notFound, PHENOLOG_HELP } from "./common.js";
import { table } from "./table.js";

const ORTHOLOGS_SHOWN = 4;

export async function render(el, id) {
  await ensure("phenotypes", "phenologs", "gene_phenotypes", "genes", "human_orthologs");
  const rows = await query(`SELECT * FROM phenotypes WHERE phenotype_id = ${lit(id)} ORDER BY n_phenologs DESC`);
  if (!rows.length) { el.innerHTML = notFound("phenotype", id); return; }
  const p = rows[0];
  const base = state.overview.release.base_species;
  const several = rows.length > 1;                        // an MP term can belong to mouse and rat; species columns only then
  const model = rows.some(r => r.species !== base);

  const [phen, genes] = await Promise.all([
    query(`SELECT rank, species, disease_id, disease_name, overlap, n_disease_orthogroups, n_phenotype_orthogroups, n_common, pval
           FROM phenologs WHERE phenotype_id = ${lit(id)} ORDER BY pval`),
    query(`SELECT gp.species, gp.gene_id, coalesce(g.symbol, gp.gene_id) AS symbol, g.family,
                  string_agg(ho.human_symbol || '|' || ho.human_gene_id, ',' ORDER BY ho.human_symbol) AS human
           FROM gene_phenotypes gp
           LEFT JOIN genes g ON g.gene_id = gp.gene_id AND g.species = gp.species
           LEFT JOIN human_orthologs ho ON ho.gene_id = gp.gene_id
           WHERE gp.phenotype_id = ${lit(id)}
           GROUP BY gp.species, gp.gene_id, g.symbol, g.family ORDER BY gp.species, symbol`),
  ]);

  el.innerHTML = `
    ${pageHead({ kicker: "Phenotype", title: esc(p.name), id: p.phenotype_id, idHref: links.monarch(p.phenotype_id),
                 links: [{ label: "Monarch", href: links.monarch(p.phenotype_id) }],
                 after: `<div class="chips" style="align-self:center">${rows.map(r => speciesChip(r.species)).join("")}</div>` })}
    ${rows.map(r => statRow([[`<i>${esc(r.species.replace(/-/g, " "))}</i>`, "species"], [r.n_genes, "genes annotated"], [r.n_orthogroups, "protein families"],
                             ...(r.species === base ? [[r.n_known_diseases, "known diseases of these genes"]]
                                                    : [[r.n_shared_orthogroups, "families shared with human"], [r.n_diseases, "human disease phenologs"],
                                                       [fmt.p(r.best_pval), "best p"]])])).join("")}
    ${model ? "" : `<p class="muted" style="margin:-6px 0 16px">Phenologs in this release compare human diseases, not human phenotypes, with the model species' phenotypes.</p>`}
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${model ? card("Phenologs", "phens", { sub: "human diseases whose protein families overlap this phenotype's more than chance" }) : ""}
      ${card("Genes annotated", "genes", { sub: model ? "the genes carrying this phenotype, and so the protein families behind every phenolog above" : "the genes carrying this phenotype" })}
    </div>`;

  if (model) {
    table(el.querySelector("#phens"), phen, [
      { key: "disease_name", label: "Human diseases", format: (v, r) => `<a href="${route.disease(r.disease_id)}">${esc(v)}</a> <span class="id faint">${esc(r.disease_id)}</span>` },
      ...(several ? [{ key: "species", label: "Species", format: v => speciesChip(v) }] : []),
      { key: "rank", label: "Phenolog rank in disease", num: true, sortDir: "asc", help: PHENOLOG_HELP.rank },
      { key: "overlap", label: "Shared families", num: true, help: PHENOLOG_HELP.overlap },
      { key: "n_disease_orthogroups", label: "Disease families", num: true, help: PHENOLOG_HELP.n_disease_orthogroups },
      { key: "n_phenotype_orthogroups", label: "Phenotype families", num: true, help: PHENOLOG_HELP.n_phenotype_orthogroups },
      { key: "pval", label: "p", num: true, sortDir: "asc", help: PHENOLOG_HELP.pval, format: v => `<span class="mono">${fmt.p(v)}</span>` },
    ], { sortKey: "pval", sortDir: "asc", search: ["disease_name", "disease_id"], placeholder: "Filter by disease name or Mondo id…",
         empty: "Not a significant phenolog of any disease", wrapHead: true });
  }

  table(el.querySelector("#genes"), genes, [
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}">${esc(v)}</a> <span class="id faint">${esc(r.gene_id)}</span>` },
    ...(several ? [{ key: "species", label: "Species", format: v => speciesChip(v) }] : []),
    { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : `<span class="faint">none</span>` },
    ...(model ? [{ key: "human", label: "Human ortholog(s)", format: orthologLinks, sortValue: r => r.human ? r.human.toLowerCase() : null }] : []),
  ], { sortKey: "symbol", sortDir: "asc", search: ["symbol", "family", "gene_id", "human"],
       placeholder: model ? "Filter by gene, family or human ortholog…" : "Filter by gene or family…" });
}

// "SYMBOL|ID,SYMBOL|ID" --> links to the human genes
function orthologLinks(v) {
  if (!v) return `<span class="faint">none</span>`;
  const items = v.split(",").map(x => x.split("|"));
  const shown = items.slice(0, ORTHOLOGS_SHOWN).map(([symbol, id]) => `<a href="${route.gene(id)}">${esc(symbol)}</a>`).join(", ");
  return items.length > ORTHOLOGS_SHOWN ? `${shown} <span class="faint">+${items.length - ORTHOLOGS_SHOWN}</span>` : shown;
}
