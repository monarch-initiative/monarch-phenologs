"""Colours, fonts and slide geometry shared by the figures and the deck builder.

The UNC template (assets/UNC_NEW_Template_2026.pptx) is 13.333 x 7.5 in with Arial throughout.
Its theme colours are Carolina Blue (accent1) and Navy (dk1); everything here is built from those
plus the dashboard's species palette, so slides, figures and the web site agree.
"""

# --- UNC template palette --------------------------------------------------------------------
NAVY = "#13294B"
CAROLINA = "#4B9CD3"
CAROLINA_DARK = "#2F7FB8"     # Carolina Blue darkened for small text on white
TINT = "#EEF5FB"              # panel background, a light tint of Carolina Blue
TINT_2 = "#DBEBF6"            # the template's lt2
WHITE = "#FFFFFF"
INK = "#13294B"               # body text
MUTED = "#5F6B7A"             # secondary text
FAINT = "#9AA5B1"             # tertiary text, hairlines
RULE = "#CBD5DF"
GREY_FILL = "#E9EDF1"

# One accent that is not blue: new candidates, non-obvious phenologs, "look here".
ACCENT = "#F26B43"            # the template's secondary orange
ACCENT_TINT = "#FDE9E2"
ACCENT_DARK = "#D04A1E"         # the accent darkened for text on a tinted panel

FONT = "Arial"

# --- slide geometry, inches ------------------------------------------------------------------
SLIDE_W, SLIDE_H = 13.333, 7.5
MARGIN_X = 0.68               # the template's title and footer start here
TITLE_Y = 0.38
LEDE_Y = 1.08                 # one-sentence takeaway under the title
LEDE_H = 0.50
FIG_X = 0.67
FIG_Y = 1.70                  # figures sit below the takeaway...
FIG_W = 12.0
FIG_H = 5.10                  # ...and end above the footer at 6.99
FIG_H_TALL = 5.62             # a slide with no takeaway line: the figure grows up into that space
FIG_DPI = 220

# --- species ---------------------------------------------------------------------------------
# Filled in from data.yaml by figures.load_data(); listed here so a typo fails loudly.
SPECIES_ORDER = ["mouse", "rat", "frog", "zebrafish", "worm", "pombe", "dicty"]
