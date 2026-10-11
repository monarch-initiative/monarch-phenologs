// One model phenotype: the diseases it is a phenolog of, and the genes annotated with it
import { ensure, query, lit } from "../db.js";
import { fmt, esc, speciesChip, route, links } from "../format.js";
import { pageHead, statRow, card, notFound } from "./common.js";
import { table } from "../charts/table.js";

export async function render(el, id) {
  await ensure("phenotypes", "phenologs", "gene_phenotypes", "genes");
  const rows = await query(`SELECT * FROM phenotypes WHERE phenotype_id = ${lit(id)} ORDER BY n_phenologs DESC`);
  if (!rows.length) { el.innerHTML = notFound("phenotype", id); return; }
  const p = rows[0];       // an MP term can belong to mouse and rat; stats are shown per species

  const [phen, genes] = await Promise.all([
    query(`SELECT rank, species, disease_id, disease_name, overlap, n_disease_orthogroups, n_phenotype_orthogroups, n_common, pval
           FROM phenologs WHERE phenotype_id = ${lit(id)} ORDER BY pval`),
    query(`SELECT gp.species, gp.gene_id, coalesce(g.symbol, gp.gene_id) AS symbol, g.family, g.n_known_diseases
           FROM gene_phenotypes gp LEFT JOIN genes g ON g.gene_id = gp.gene_id AND g.species = gp.species
           WHERE gp.phenotype_id = ${lit(id)} ORDER BY gp.species, symbol`),
  ]);

  el.innerHTML = `
    ${pageHead({ kicker: "Phenotype", title: esc(p.name), id: p.phenotype_id, idHref: links.monarch(p.phenotype_id),
                 links: [{ label: "Monarch", href: links.monarch(p.phenotype_id) }],
                 after: `<div class="chips" style="align-self:center">${rows.map(r => speciesChip(r.species)).join("")}</div>` })}
    ${rows.map(r => statRow([[`<i>${esc(r.species.replace(/-/g, " "))}</i>`, "species"], [r.n_genes, "genes annotated"], [r.n_orthogroups, "protein families"],
                             [r.n_phenologs, "phenologs"], [r.n_diseases, "diseases"], [fmt.p(r.best_pval), "best p"]])).join("")}
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${card("Phenologs", "phens", { sub: "human diseases whose protein families overlap this phenotype's more than chance" })}
      ${card("Genes annotated", "genes", { sub: "the genes carrying this phenotype, and so the protein families behind every phenolog above" })}
    </div>`;

  table(el.querySelector("#phens"), phen, [
    { key: "disease_name", label: "Disease", format: (v, r) => `<a href="${route.disease(r.disease_id)}">${esc(v)}</a> <span class="id faint">${esc(r.disease_id)}</span>` },
    { key: "species", label: "Species", format: v => speciesChip(v) },
    { key: "rank", label: "Rank in disease", num: true },
    { key: "overlap", label: "Overlap", num: true },
    { key: "n_disease_orthogroups", label: "Disease fam.", num: true },
    { key: "n_phenotype_orthogroups", label: "Phenotype fam.", num: true },
    { key: "pval", label: "p", num: true, format: v => `<span class="mono">${fmt.p(v)}</span>` },
  ], { sortKey: "pval", sortDir: "asc", search: ["disease_name", "disease_id"], empty: "Not a significant phenolog of any disease" });

  table(el.querySelector("#genes"), genes, [
    { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}">${esc(v)}</a> <span class="id faint">${esc(r.gene_id)}</span>` },
    { key: "species", label: "Species", format: v => speciesChip(v) },
    { key: "family", label: "Family", format: v => v ? `<a class="mono" href="${route.family(v)}">${esc(v)}</a>` : `<span class="faint">none</span>` },
  ], { sortKey: "symbol", sortDir: "asc", search: ["symbol", "family", "gene_id"] });
}
