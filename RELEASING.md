# Releasing

There are two kinds of release:

- **Code.** The container image `ghcr.io/monarch-initiative/monarch-phenologs:<date and time>`, which holds the code
  and its locked environment. It is built from the Actions tab whenever the code changes, which is rarely.
- **Data.** One release per Monarch KG version, made by running the image on that KG. It holds both prediction
  networks, a folder each: `disease` (human diseases against the model species' phenotypes) and `phenotype` (human
  phenotypes against them). Published to `gs://data-public-monarchinitiative/phenologs/alpha/<kg version>/<network>/`
  (`https://data.monarchinitiative.org/phenologs/alpha/`).
  - The dashboard (GitHub Pages) lists the data releases from `releases.json` there. A new data release needs no
    site deploy.
  - Each data release records the KG, PHENIO and code versions that made it, the code as the image's tag and
    commit (`2026-10-19-2315 (1a2b3c4)`).

## Making a data release

Run the image's steps in order on a fresh project directory. Scripts are at `/app/python` and run from any working
directory. No versions need to be given:

- Step 01 takes the newest Monarch KG release (`monarch-kg/latest/metadata.yaml`) and the PHENIO it was built with.
- It writes both to `monarch_kg/versions.json`, where the later steps read them.
- `-gv 2026-03-06` asks for another KG. A date Monarch never released is looked for in `monarch-kg-dev/`.

Steps 01 and 02 run once. Steps 03 to 07 run once per network, with `$PRD` set to `disease` and then to `phenotype`.

| step | command | needs (old cluster settings, a starting point) |
|---|---|---|
| 01 | `01_download_and_setup.py -p $PROJECT` | internet, 8 GB, ~6 GB disk |
| 02 | `02_kg_data_extraction.py -p $PROJECT` | 8 GB |
| 03 | `03_compute_phenologs_randomized_trials.py -p $PROJECT -n 1000 -c $CPUS -taxon_id 9606 -prd $PRD -seed 1` | 128 cpus, 256 GB, 6 h |
| 04 | `04_compute_phenologs.py -p $PROJECT -c $CPUS -taxon_id 9606 -prd $PRD` | 8 cpus, 32 GB |
| 05 | `05_compute_phenologs_fdr_data.py -p $PROJECT -taxon_id 9606 -prd $PRD` | 8 GB |
| 06 | `06_phenologs_compute_ortholog_to_phenotype_rankings.py -p $PROJECT -taxon_id 9606 -prd $PRD -fdr 0.05 -kneighbs 1000 -rank_metric hg -rank_type gene` | 8 GB; phenotype 8 GB, 30 min |
| 06, candidates | `phenologs_disease_gene_candidate_rankings.py -p $PROJECT -taxon_id 9606 -prd $PRD -fdr 0.05 -kneighbs 1000 -rank_metric hg` | 8 GB; phenotype 8 GB, 10 min, 6 GB disk |
| 07 | `07_build_release.py -p $PROJECT -prd $PRD -c $CPUS -m <memory, e.g. 8GB>` | 4 cpus, 8 GB, under a minute; phenotype 4 cpus, 12 GB, under a minute |

- **Order.** 03 and 04 both need only 02's output, so they can run side by side. 05 needs both. The two networks
  share only 01 and 02, so they can run side by side too.
- **Output.** `$PROJECT/release/<kg version>/disease/` (about 30 MB) and `.../phenotype/` (about 55 MB). These
  folders are the release; everything else in `$PROJECT` is working data that can be deleted once it is published.
- **Phenotype candidates.** A human phenotype's ranking scores thousands of genes (32 million pairs for 2026-03-06),
  so 07 keeps each human phenotype's best 100 in the release (`-top_candidates`). The full lists stay in `$PROJECT`.
- **Repeats.** The trials are seeded, so the same KG, code and seed give the same release.
- **Measured resource peaks.** The paper runs' Nextflow traces on the cluster (`PHENOLOGS_RESEARCH/<pipeline>/runs/`)
  have the measured peaks. The phenotype numbers for 06 and 07 were measured on a laptop, for 2026-03-06.
- **Checking a release before publishing.** Copy the release folder to a machine with a browser, then:

      python dashboard/serve.py -p $PROJECT      # http://localhost:8765

## Publishing a data release

Needs write access to `gs://data-public-monarchinitiative/phenologs/`. Order matters:

```
python /app/python/08_update_release_list.py -p $PROJECT            # the very first release: add -r none
ALPHA=gs://data-public-monarchinitiative/phenologs/alpha
for PRD in disease phenotype; do
    gcloud storage rsync -r $PROJECT/release/2026-03-06/$PRD $ALPHA/2026-03-06/$PRD
    gcloud storage rsync -r --delete-unmatched-destination-objects $ALPHA/2026-03-06/$PRD $ALPHA/latest/$PRD
done
gcloud storage cp $PROJECT/release/releases.json $ALPHA/releases.json
```

1. **08** fetches the published `releases.json`, puts the project's releases in it (in place of earlier builds of the
   same folder), and writes `$PROJECT/release/releases.json`. Run it right before uploading, so the list it adds to
   is current. It stops if the published list can't be read, so a failed fetch never drops the published releases.
2. **The release folders** go up before the list, so the site never lists a release that isn't there yet.
3. **`latest/<network>/`** is a copy of the newest release of that network, a stable address for downloads. 08 says
   which release is the newest; refresh `latest/` only when it is the one just built (an older KG rebuilt later is not).
4. **`releases.json`** goes up last.

Things to know:

- **Skip.** If the published `releases.json` already lists `<kg version>/disease/` and `<kg version>/phenotype/`,
  that release exists, so nothing needs to run. If it lists one, only the other network needs to run.
- **The alpha folder** also holds older files of other kinds. Never sync it as a whole with deletion.
- **Folder index pages** on data.monarchinitiative.org are regenerated by Monarch's indexer. The site doesn't need them.

## Making an image

- **Building.** Actions tab, Build Image, Run workflow, from any branch. Each build is tagged with its date and time in
  UTC, e.g. `2026-10-19-2315`; the run's summary shows its tags. Nothing is overwritten and nothing needs bumping.
- **latest** moves to every build from main.
- **Trying a branch out.** Give an extra tag (e.g. `dev`) and it moves to the new image too, so a test pipeline can
  follow `dev` while the release pipeline stays pinned. Old images can be deleted in the package settings.
- **Pulling.** The package is private for now (only the org's owners can make it public), so the cluster logs in
  first, with a GitHub classic personal access token that has `read:packages`:

      apptainer registry login --username <github user> docker://ghcr.io      # asks for the token
      apptainer pull docker://ghcr.io/monarch-initiative/monarch-phenologs:2026-10-19-2315

  A job that can't answer a prompt sets `APPTAINER_DOCKER_USERNAME` and `APPTAINER_DOCKER_PASSWORD` (the token)
  instead. Once the package is public, no login is needed.
- **Pinning.** The release pipeline pins a date tag, which never moves. Moving to new code is a one-line change in its
  settings.
