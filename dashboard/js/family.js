// One PANTHER protein family: its members across species, and the diseases or human phenotypes it is implicated in
import { ensure, query, lit } from "./db.js";
import { state, fmt, esc, speciesChip, route, links } from "./format.js";
import { pageHead, statRow, card, notFound, legend, nameGrid, HELP } from "./common.js";
import { table } from "./table.js";

const MEMBERS_SHOWN = 150;    // per species, before "+N more"

export async function render(el, id) {
  const net = state.net;
  await ensure("families", "genes", "candidates", "terms", "known_associations");
  const [f] = await query(`SELECT * FROM families WHERE family = ${lit(id)}`);
  if (!f) { el.innerHTML = notFound("protein family", id); return; }

  const [members, predicted, known] = await Promise.all([
    query(`SELECT species, gene_id, symbol, n_phenotypes, n_known_terms FROM genes WHERE family = ${lit(id)} ORDER BY species, symbol`),
    // a term is known for the family when one of its human genes is curated for it, scored or not
    query(`SELECT c.term_id, d.name, min(c.score) AS best_score, min(c.rank) AS best_rank, d.n_candidates,
                  c.term_id IN (SELECT term_id FROM known_associations WHERE family = ${lit(id)}) AS known
           FROM candidates c JOIN terms d USING (term_id) WHERE c.family = ${lit(id)}
           GROUP BY c.term_id, d.name, d.n_candidates ORDER BY best_score`),
    query(`SELECT DISTINCT gene_id, symbol, term_id, term_name FROM known_associations WHERE family = ${lit(id)} ORDER BY term_name, symbol`),
  ]);
  const bySp = d3.group(members, m => m.species);
  const order = state.overview.species.map(s => s.species);
  const speciesList = [...bySp.keys()].sort((a, b) => order.indexOf(a) - order.indexOf(b));

  el.innerHTML = `
    ${pageHead({ kicker: "Protein family", title: `<span class="mono" style="font-size:inherit">${esc(f.family)}</span>`,
                 links: [{ label: "PANTHER", href: links.panther(f.family) }] })}
    ${statRow([[f.n_genes, "genes"], [f.n_species, "species"], [f.n_base_genes, "human genes"], [f.n_known_terms, `known ${net.terms}`],
               [f.n_predicted_terms, `${net.terms} ranked for`], [f.n_novel_terms, "of them novel"], [fmt.p(f.best_score), "best score"]])}
    <div class="grid" style="grid-template-columns:1fr; gap:16px">
      ${card("Members by species", "members", { sub: legend([["var(--accent)", `human gene with a curated ${net.term}`], ["var(--text)", `no curated ${net.term}`]]) })}
      ${card(`${net.Terms} this family is ranked for`, "pred", { sub: `${legend([["var(--accent)", `a human member is a known gene of the ${net.term}`], ["var(--accent-2)", "novel for the family"]])}<span class="faint">best rank and score of any member gene</span>` })}
      ${known.length ? card("Known associations", "known", { sub: `curated ${net.pair} pairs involving this family` }) : ""}
    </div>`;

  el.querySelector("#members").innerHTML = speciesList.map(sp => {
    const ms = bySp.get(sp);
    return `<div class="member-row"><div>${speciesChip(sp, { count: ms.length })}</div>
      <div>${nameGrid(ms.slice(0, MEMBERS_SHOWN).map(m => ({ label: m.symbol, href: route.gene(m.gene_id), cls: m.n_known_terms ? "known-mark" : "plain",
                                                             title: `${m.symbol} · ${m.n_phenotypes} phenotypes` })))}
        ${ms.length > MEMBERS_SHOWN ? `<span class="faint">+${ms.length - MEMBERS_SHOWN} more</span>` : ""}</div></div>`;
  }).join("");

  table(el.querySelector("#pred"), predicted, [
    { key: "name", label: net.human, format: (v, r) => `<a href="${route.term(r.term_id)}" class="${r.known ? "known-mark" : ""}">${esc(v)}</a>` },
    { key: "best_rank", label: `Best rank in ${net.term}`, num: true, sortDir: "asc", help: `Best place of any member gene in the ${net.term}'s ranking, out of every gene it scored`,
      format: (v, r) => `${fmt.int(v)} <span class="faint">/ ${fmt.int(r.n_candidates)}</span>` },
    { key: "best_score", label: "Best score", num: true, sortDir: "asc", help: HELP.score, format: v => `<span class="mono">${fmt.p(v)}</span>` },
  ], { sortKey: "best_score", sortDir: "asc", search: ["name", "term_id"], placeholder: `Filter by ${net.term} name or ${net.ontology} id…`,
       empty: "Not in any candidate ranking", download: `${f.family}_${net.browse}` });

  if (known.length) {
    table(el.querySelector("#known"), known, [
      { key: "symbol", label: "Gene", format: (v, r) => `<a href="${route.gene(r.gene_id)}">${esc(v)}</a>` },
      { key: "term_name", label: net.human, format: (v, r) => `<a href="${route.term(r.term_id)}">${esc(v)}</a> <span class="id faint">${esc(r.term_id)}</span>` },
    ], { sortKey: "term_name", sortDir: "asc", search: ["symbol", "term_name", "term_id"], placeholder: `Filter by gene or ${net.term}…`, pageSize: 15,
         download: `${f.family}_known` });
  }
}
