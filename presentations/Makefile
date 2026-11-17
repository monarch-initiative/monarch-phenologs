# Phenologs talk. `make setup` once, then `make deck`.
PY   ?= .venv/bin/python
DECK := out/phenologs_talk.pptx

.PHONY: setup deck text figures crops pdf sheets screenshot clean

setup:            ## slide tooling in its own environment (the pipeline's uv.lock is left alone)
	uv venv --python 3.11 .venv
	uv pip install --python $(PY) -r requirements.txt

deck:             ## figures + deck -> out/phenologs_talk.pptx
	$(PY) -m phenologs_talk.build

text:             ## deck only, reusing the figures already rendered (after editing words in deck.yaml)
	$(PY) -m phenologs_talk.build --no-figures

figures:          ## figures only -> out/figures/
	$(PY) -m phenologs_talk.build figures

crops:            ## the panels cut from the preprint -> out/crops/ (to check the boxes in paper.py)
	$(PY) -m phenologs_talk.build crops

pdf: deck         ## PDF and one PNG per slide, through PowerPoint (macOS)
	tools/export_pdf.sh $(DECK)

sheets: pdf       ## contact sheets of every slide -> out/sheets/
	$(PY) tools/contact_sheet.py 3 9

screenshot:       ## refresh assets/dashboard_bbs.png from the live dashboard (Chrome, Node 22+)
	node tools/screenshot_dashboard.mjs

clean:
	rm -rf out
