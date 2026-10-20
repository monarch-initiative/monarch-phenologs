// One gene. A base species gene: the diseases or human phenotypes it is known for and ranked for, its phenotypes, its protein family.
// A model species gene: its phenotypes, its human orthologs, its protein family
import { ensure, query, lit } from "./db.js";
import { state, fmt, esc, speciesOf, speciesChip, supportChips, route, links } from "./format.js";
import { pageHead, statRow, card, notFound, legend, nameGrid, pairKey, recoveryKey, supportLabel, HELP } from "./common.js";
import { table } from "./table.js";

const MEMBERS_SHOWN = 60;    // per species, before "+N more"

export async function render(el, id) {
  const net = state.net, top = state.overview.release.top_candidates;
  await ensure("genes", "known_associations", "candidates", "gene_phenotypes", "terms", "human_orthologs");
  let rows = await query(`SELECT * FROM genes WHERE gene_id = ${lit(id)}`);
  if (!rows.length) {
    // allow a symbol in the url, base species first
    rows = await query(`SELECT * FROM genes WHERE upper(symbol) = upper(${lit(id)}) ORDER BY (species = ${lit(state.overview.release.base_species)}) DESC LIMIT 1`);
    if (rows.length) { location.replace(route.gene(rows[0].gene_id)); return; }
    el.innerHTML = notFound("gene", id); return;
  }
  const g = rows[0];
  const isBase = g.species === state.overview.release.base_species;

  const [known, ranked, phenos, members, human] = await Promise.all([
    query(`SELECT k.term_id, any_value(k.term_name) AS term_name, min(k.rank) AS rank,
                  k.term_id IN (SELECT term_id FROM terms) AS has_phenologs
           FROM known_associations k WHERE k.gene_id = ${lit(g.gene_id)} GROUP BY k.term_id ORDER BY term_name`),
    query(`SELECT c.rank, c.term_id, d.name, c.score, c.known, c.support, c.n_phenologs, d.n_candidates
           FROM candidates c JOIN terms d USING (term_id) WHERE c.gene_id = ${lit(g.gene_id)} ORDER BY c.score, d.name`),
    query(`SELECT phenotype_id, phenotype_name FROM gene_phenotypes WHERE gene_id = ${lit(g.gene_id)} ORDER BY phenotype_name`),
    g.family ? query(`SELECT species, gene_id, symbol FROM genes WHERE family = ${lit(g.family)} ORDER BY species, symbol`) : [],
    isBase ? [] : query(`SELECT human_gene_id, human_symbol FROM human_orthologs WHERE gene_id = ${lit(g.gene_id)} ORDER BY human_symbol`),
  ]);

  const stats = isBase
    ? [[g.n_phenotypes, "phenotypes"], [g.n_known_terms, `known ${net.terms}`], [g.n_predicted_terms, `${net.terms} ranked for`],
       [g.n_novel_terms, "of them novel"], [g.best_rank ?? "–", `best rank in a ${net.term}`]]
    : [[g.n_phenotypes, "phenotypes"], [human.length, "human orthologs"], [g.n_ortholog_species, "species with an ortholog"]];

  el.innerHTML = `
    ${pageHead({ kicker: "Gene", title: `${esc(g.symbol)} <span class="muted" style="font-weight:400;font-size:18px">· <i>${esc(speciesOf(g.species).binomial)}</i></span>`,
                 id: g.gene_id, idHref: links.monarch(g.gene_id),
                 links: [{ label: "Monarch", href: links.monarch(g.gene_id) }, g.family ? { label: "PANTHER " + g.family, href: links.panther(g.family) } : null],
                 after: `<div style="align-self:center">${speciesChip(g.species)} ${g.family ? `<a class="chip mono" href="${route.family(g.family)}">${esc(g.family)}</a>` : ""}</div>` })}
    ${statRow(stats)}
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${isBase ? card(`Known ${net.terms}`, "known", { sub: known.length ? legend([...recoveryKey(), ["var(--faint)", "no phenolog associations"]]) : "" }) : ""}
      ${isBase ? card(`${net.Terms} this gene is ranked for`, "ranked", { sub: `${legend(pairKey())}<span class="faint">${top
          ? `every ${net.term} whose best ${fmt.int(top)} candidates include this gene` : `every ${net.term} whose candidate ranking scored this gene`}</span>` }) : ""}
      ${isBase ? "" : card("Human orthologs", "human", { sub: "ortholog edges to human genes in the graph" })}
      ${card("Phenotypes", "phenos", { sub: `${esc(speciesOf(g.species).ontology || "phenotype")} annotations of this gene` })}
      ${g.family ? card("Protein family " + g.family, "members", { sub: "every gene of the family in this release" }) : ""}
    </div>`;

  if (isBase) {
    el.querySelector("#known").innerHTML = known.length
      ? nameGrid(known.map(k => k.has_phenologs
          ? { label: k.term_name, href: route.term(k.term_id), cls: k.rank != null ? "known-mark" : "plain", note: k.rank != null ? `#${fmt.int(k.rank)}` : "" }
          : { label: k.term_name, href: links.monarch(k.term_id), cls: "out", ext: true }), { wide: true })
      : `<div class="empty">No curated ${net.term} associations</div>`;

    table(el.querySelector("#ranked"), ranked, [
      { key: "rank", label: `Rank in ${net.term}`, num: true, sortDir: "asc", help: HELP.rank, format: (v, r) => `${fmt.int(v)} <span class="faint">/ ${fmt.int(r.n_candidates)}</span>` },
      { key: "name", label: net.human, format: (v, r) => `<a href="${route.term(r.term_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
      { key: "score", label: "Score", num: true, sortDir: "asc", help: HELP.score, format: v => `<span class="mono">${fmt.p(v)}</span>` },
      { key: "n_phenologs", label: supportLabel(), sortDir: "desc", help: HELP.support, format: (v, r) => supportChips(r.support, g.symbol) },
    ], { sortKey: "score", sortDir: "asc", search: ["name", "term_id"], placeholder: `Filter by ${net.term} name or ${net.ontology} id…`,
         download: `${g.symbol}_${net.browse}`,
         empty: g.in_reach ? (top ? `Not among the best ${fmt.int(top)} candidates of any ${net.term}` : `No ${net.term}'s phenologs reach this gene`)
              : g.n_ortholog_species ? `No ortholog of ${g.symbol} in the model species carries a phenotype, so no ${net.term}'s ranking can reach it`
              : `${g.symbol} has no ortholog in the model species, so no ${net.term}'s ranking can reach it` });
  } else {
    el.querySelector("#human").innerHTML = human.length
      ? nameGrid(human.map(h => ({ label: h.human_symbol, href: route.gene(h.human_gene_id) })))
      : `<div class="empty">No human ortholog in the graph</div>`;
  }

  table(el.querySelector("#phenos"), phenos, [
    { key: "phenotype_name", label: "Phenotype", format: (v, r) => `<a href="${route.phenotype(r.phenotype_id)}">${esc(v)}</a>` },
    { key: "phenotype_id", label: "Id", format: v => `<span class="id muted">${esc(v)}</span>` },
  ], { sortKey: "phenotype_name", sortDir: "asc", search: ["phenotype_name", "phenotype_id"], pageSize: 15, empty: "No phenotype annotations",
       download: `${g.symbol}_phenotypes` });

  if (g.family) {
    const bySp = d3.group(members, m => m.species);
    const order = state.overview.species.map(s => s.species);
    el.querySelector("#members").innerHTML = [...bySp.keys()].sort((a, b) => order.indexOf(a) - order.indexOf(b)).map(sp => {
      const ms = bySp.get(sp);
      return `<div class="member-row"><div>${speciesChip(sp, { count: ms.length })}</div>
        <div>${nameGrid(ms.slice(0, MEMBERS_SHOWN).map(m => ({ label: m.symbol, href: route.gene(m.gene_id), cls: m.gene_id === g.gene_id ? "current" : "plain" })))}
          ${ms.length > MEMBERS_SHOWN ? `<a class="faint" href="${route.family(g.family)}">+${ms.length - MEMBERS_SHOWN} more</a>` : ""}</div></div>`;
    }).join("");
  }
}
