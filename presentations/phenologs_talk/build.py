"""Build the talk: data.yaml + deck.yaml + figures -> out/phenologs_talk.pptx on the UNC template.

    python -m phenologs_talk.build            # figures, then the deck
    python -m phenologs_talk.build --no-figures   # reuse out/figures/ (text-only edits)
    python -m phenologs_talk.build crops      # write the preprint crops to out/crops/ and stop

Titles, ledes, footers, cards and notes are real text boxes, so they stay editable in PowerPoint
and survive the import into Google Slides. Diagrams go in as pictures.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import shutil
import string
import sys
from datetime import datetime
from pathlib import Path

import yaml
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

from . import figures as F
from . import paper, render_figures
from . import theme as T
from .draw import sci_plain

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TEMPLATE = ROOT / "assets" / "UNC_NEW_Template_2026.pptx"
OUT = ROOT / "out"

LAYOUT = {
    "cover": "Cover with Full Signature",
    "divider": {"carolina": "Divider with Subtitle Carolina Blue", "navy": "Divider with Subtitle Navy",
                "white": "Divider with Subtitle White"},
    "content": "Content with Blank Area",
    "end": "End Slide with Signature White",
}
# placeholder idx values in the template's layouts
PH_COVER_SUB, PH_COVER_DATE, PH_DIVIDER_LABEL, PH_FOOTER = 1, 10, 16, 13
TITLE_PT, LEDE_PT, LEDE_MAX = 28, 17, 104


def rgb(hex_color: str) -> RGBColor:
    return RGBColor.from_string(hex_color.lstrip("#").upper())


class _Fmt(string.Formatter):
    """str.format with one extra conversion: {x!p} writes a P value as 1.3 × 10⁻¹⁸."""

    def convert_field(self, value, conversion):
        if conversion == "p":
            return sci_plain(float(value))
        return super().convert_field(value, conversion)


def fill(text, D):
    return _Fmt().vformat(text, (), D) if isinstance(text, str) and "{" in text else text


# ---------------------------------------------------------------------------------------------
# template handling
# ---------------------------------------------------------------------------------------------
def open_template(path: Path):
    """Open the template and drop its example slides, keeping the master and layouts."""
    prs = Presentation(str(path))
    ids = prs.slides._sldIdLst
    for sld in list(ids):
        prs.part.drop_rel(sld.rId)
        ids.remove(sld)
    return prs


def add_slide(prs, layout_name):
    layout = next((l for l in prs.slide_layouts if l.name == layout_name), None)
    if layout is None:
        raise KeyError(f"layout {layout_name!r} is not in the template: {[l.name for l in prs.slide_layouts]}")
    return prs.slides.add_slide(layout)


def add_slide_number(slide):
    """python-pptx leaves the layout's slide-number placeholder behind; copy it onto the slide."""
    for ph in slide.slide_layout.placeholders:
        if "SLIDE_NUMBER" in str(ph.placeholder_format.type):
            slide.shapes._spTree.append(copy.deepcopy(ph._element))


def placeholder(slide, idx):
    return next((p for p in slide.placeholders if p.placeholder_format.idx == idx), None)


def drop_empty_placeholders(slide):
    for p in list(slide.placeholders):
        if p.placeholder_format.type is not None and p.has_text_frame and not p.text_frame.text.strip():
            if "SLIDE_NUMBER" in str(p.placeholder_format.type):
                continue
            p._element.getparent().remove(p._element)


def set_text(shape, text, size=None, bold=None, color=None, italic=None):
    tf = shape.text_frame
    tf.text = text
    for para in tf.paragraphs:
        for run in para.runs:
            if size:
                run.font.size = Pt(size)
            if bold is not None:
                run.font.bold = bold
            if italic is not None:
                run.font.italic = italic
            if color:
                run.font.color.rgb = rgb(color)
    return tf


def textbox(slide, x, y, w, h, text, size=14, bold=False, color=T.INK, italic=False, align=PP_ALIGN.LEFT,
            anchor=MSO_ANCHOR.TOP, name=None):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    para = tf.paragraphs[0]
    para.alignment = align
    run = para.add_run()
    run.text = text
    run.font.size, run.font.bold, run.font.italic = Pt(size), bold, italic
    run.font.color.rgb = rgb(color)
    return tb


def add_para(tf, text, size, bold=False, color=T.INK, space_before=0, italic=False):
    para = tf.add_paragraph() if tf.paragraphs[0].runs else tf.paragraphs[0]
    para.space_before = Pt(space_before)
    run = para.add_run()
    run.text = text
    run.font.size, run.font.bold, run.font.italic = Pt(size), bold, italic
    run.font.color.rgb = rgb(color)
    return para


def bullet(para, indent=0.24):
    """Give a paragraph a real bullet with a hanging indent."""
    pPr = para._p.get_or_add_pPr()
    pPr.set("marL", str(int(Inches(indent))))
    pPr.set("indent", str(-int(Inches(indent))))
    for tag, attrs in (("a:buFont", {"typeface": "Arial"}), ("a:buChar", {"char": "•"})):
        pPr.append(pPr.makeelement(qn(tag), attrs))
    return para


def card(slide, x, y, w, h, fill_hex=T.TINT, name=None):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.adjustments[0] = min(0.5, 0.14 / min(w, h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = rgb(fill_hex)
    shp.line.fill.background()
    shp.shadow.inherit = False
    if name:
        shp.name = name
    tf = shp.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.28)
    tf.margin_top = tf.margin_bottom = Inches(0.22)
    tf.vertical_anchor = MSO_ANCHOR.TOP
    tf.paragraphs[0].alignment = PP_ALIGN.LEFT
    return shp


def autofit(shape, font_scale=None, line_reduction=None):
    """Write the shrink PowerPoint itself applies when it opens the deck and the text does not fit its placeholder,
    so the built file looks the same before and after it is saved from PowerPoint. Values are percentages."""
    body = shape.text_frame._txBody.find(qn("a:bodyPr"))
    for child in list(body):
        if child.tag in (qn("a:normAutofit"), qn("a:noAutofit"), qn("a:spAutoFit")):
            body.remove(child)
    fit = body.makeelement(qn("a:normAutofit"), {})
    if font_scale:
        fit.set("fontScale", str(int(font_scale * 1000)))
    if line_reduction:
        fit.set("lnSpcReduction", str(int(line_reduction * 1000)))
    body.append(fit)


def set_alt(picture, text):
    picture._element.nvPicPr.cNvPr.set("descr", text)


# ---------------------------------------------------------------------------------------------
# slide kinds
# ---------------------------------------------------------------------------------------------
def content_frame(prs, spec, meta, D, lede=None):
    """Title, lede and footer on the template's blank content layout."""
    slide = add_slide(prs, LAYOUT["content"])
    title = fill(spec["title"], D)
    set_text(slide.shapes.title, title.replace("\n", "\v"), size=spec.get("title_size", TITLE_PT))
    if "\n" in title:           # a two-line title: grow the box so the first line stays where one-line titles sit
        t = slide.shapes.title                # position is inherited from the layout: set all four, or width becomes 0
        t.left, t.top, t.width, t.height = t.left, t.top, t.width, t.height + Inches(0.47) * title.count("\n")
        if lede is not None or spec.get("lede"):
            print(f"  note: slide {spec['id']!r} has a two-line title and a lede; they will overlap")
    if "title_box" in spec:     # the title placeholder as moved by hand: x, y, width, height in inches
        t = slide.shapes.title
        t.left, t.top, t.width, t.height = (Inches(v) for v in spec["title_box"])
    lede = lede if lede is not None else spec.get("lede")
    if lede:
        text = fill(lede, D)
        by_hand = "lede_box" in spec or spec.get("figure_dy", 0) > 0      # laid out in PowerPoint and carried over
        most = int(LEDE_MAX * LEDE_PT / spec.get("lede_size", LEDE_PT))    # a smaller lede holds more on its one line
        if len(text) > most and not by_hand:
            print(f"  note: lede on {spec['id']!r} is {len(text)} characters and may wrap onto the figure (keep under {most})")
        x, y, w = spec.get("lede_box", (T.MARGIN_X + 0.1, T.LEDE_Y, 11.4))
        textbox(slide, x, y, w, T.LEDE_H, text, size=spec.get("lede_size", LEDE_PT), color="#3E4C5E", name="Lede")
    foot = placeholder(slide, PH_FOOTER)
    if foot is not None:
        section = spec.get("section")
        set_text(foot, (f"{section}   ·   " if section else "") + fill(meta["footer"], D))
        autofit(foot, line_reduction=10)          # what PowerPoint settles on for the one-line footer
    return slide


def slide_cover(prs, spec, meta, D, figs):
    slide = add_slide(prs, LAYOUT["cover"])
    # a newline in cover_title is a line break inside the title, not a new paragraph
    set_text(slide.shapes.title, fill(meta.get("cover_title", meta["title"]), D).replace("\n", "\v"), size=60)
    sub = placeholder(slide, PH_COVER_SUB)
    tf = sub.text_frame
    tf.text = ""
    add_para(tf, fill(meta["subtitle"], D), 20, bold=True, color=T.NAVY)
    add_para(tf, fill(meta["presenter"], D), 15, bold=False, color=T.NAVY, space_before=6)
    # two lines do not fit the template's half-inch box, so PowerPoint shows them at 70% (14 and 10.5 pt)
    autofit(sub, font_scale=70, line_reduction=20)
    set_text(placeholder(slide, PH_COVER_DATE), fill(meta["date"], D))
    return [slide]


def slide_divider(prs, spec, meta, D, figs):
    slide = add_slide(prs, LAYOUT["divider"][spec.get("style", "carolina")])
    set_text(slide.shapes.title, fill(spec["title"], D), size=spec.get("size", 44))
    lab = placeholder(slide, PH_DIVIDER_LABEL)
    if spec.get("label"):
        set_text(lab, fill(spec["label"], D))
    if spec.get("sub"):         # one line under the title; `link` makes its last part clickable
        t = slide.shapes.title                          # the title is centred in its box, so sit just below its text
        size = spec.get("size", 44)
        lines = max(1, -(-len(fill(spec["title"], D)) // 30))
        y = t.top + t.height // 2 + Inches(lines * size * 1.2 / 72 / 2 + 0.3)
        tb = textbox(slide, 0, 0, 1, 0.45, fill(spec["sub"], D) + ("  " if spec.get("link") else ""), size=18, bold=True,
                     color=T.NAVY, name="Divider sub")
        tb.left, tb.top, tb.width = t.left, y, t.width
        if spec.get("link"):
            url = fill(spec["link"], D)
            run = tb.text_frame.paragraphs[0].add_run()
            run.text = url.split("://", 1)[-1]
            run.font.size, run.font.bold = Pt(18), False
            run.hyperlink.address = url
    return [slide]


def slide_figure(prs, spec, meta, D, figs):
    slides = []
    steps = spec.get("steps") or [{}]
    for st in steps:
        key = f"{spec['figure']}-{st['step']}" if "step" in st else spec["figure"]
        if key not in figs:
            raise KeyError(f"slide {spec['id']!r} wants figure {key!r}, which was not rendered")
        slide = content_frame(prs, spec, meta, D, lede=st.get("lede"))
        pic = slide.shapes.add_picture(str(figs[key]), Inches(T.FIG_X), Inches(T.FIG_Y), width=Inches(T.FIG_W))
        # a taller canvas (theme.FIG_H_TALL) grows upward; figure_dy moves the figure down (up when negative)
        pic.top = Inches(T.FIG_Y + T.FIG_H + spec.get("figure_dy", 0)) - pic.height
        if pic.top < Inches(T.FIG_Y - 0.01) and (st.get("lede") or spec.get("lede")) and "lede_box" not in spec:
            print(f"  note: figure {key!r} is taller than the standard frame and slide {spec['id']!r} has a lede under it")
        pic.name = f"Figure: {key}"
        set_alt(pic, F.REGISTRY[spec["figure"]][1])
        slides.append(slide)
    return slides


def slide_cards(prs, spec, meta, D, figs):
    slide = content_frame(prs, spec, meta, D)
    cards = spec["cards"]
    cols = spec.get("columns", 2)
    rows = -(-len(cards) // cols)
    gap = 0.3
    x0, y0, W, H = T.FIG_X, T.FIG_Y + 0.05, T.FIG_W, T.FIG_H - 0.15
    cw, ch = (W - gap * (cols - 1)) / cols, (H - gap * (rows - 1)) / rows
    for i, c in enumerate(cards):
        r, k = divmod(i, cols)
        shp = card(slide, x0 + k * (cw + gap), y0 + r * (ch + gap), cw, ch, name=f"Card {i + 1}")
        tf = shp.text_frame
        if c.get("tag"):
            add_para(tf, fill(c["tag"], D), 11, bold=True, color=T.CAROLINA_DARK)
        if c.get("head"):
            add_para(tf, fill(c["head"], D), 20, bold=True, color=T.NAVY, space_before=4)
        if c.get("body"):
            add_para(tf, fill(c["body"], D), 14.5, color="#3E4C5E", space_before=8)
        for item in c.get("items", []):                 # a card can be a short bulleted list instead
            bullet(add_para(tf, fill(item, D), 16, color=T.NAVY, space_before=8))
    return [slide]


def slide_stats(prs, spec, meta, D, figs):
    slide = content_frame(prs, spec, meta, D)
    stats = spec["stats"]
    gap = 0.3
    x0, y0, W = T.FIG_X, T.FIG_Y - 0.25, T.FIG_W
    cw, ch = (W - gap * (len(stats) - 1)) / len(stats), 2.75
    for i, s in enumerate(stats):
        first = i == 0
        shp = card(slide, x0 + i * (cw + gap), y0, cw, ch, fill_hex=T.NAVY if first else T.TINT, name=f"Stat {i + 1}")
        tf = shp.text_frame
        add_para(tf, fill(s["big"], D), 40, bold=True, color=T.WHITE if first else T.NAVY)
        add_para(tf, fill(s["head"], D), 18, bold=True, color=T.WHITE if first else T.NAVY, space_before=10)
        add_para(tf, fill(s["body"], D), 14, color="#DCE8F5" if first else "#3E4C5E", space_before=8)
    if spec.get("closing"):
        textbox(slide, x0, y0 + ch + 0.4, W, 1.0, fill(spec["closing"], D), size=20, bold=True, color=T.NAVY, name="Closing")
    return [slide]


def slide_columns(prs, spec, meta, D, figs):
    slide = content_frame(prs, spec, meta, D)
    gap = 0.4
    cw = (T.FIG_W - gap) / 2
    for k, side in enumerate(["left", "right"]):
        col = spec[side]
        shp = card(slide, T.FIG_X + k * (cw + gap), T.LEDE_Y + 0.3, cw, 4.4, name=f"Column {side}")
        tf = shp.text_frame
        add_para(tf, fill(col["head"], D), 18, bold=True, color=T.NAVY)
        for item in col["items"]:
            add_para(tf, fill(item, D), 14.5, color="#3E4C5E", space_before=8)
    return [slide]


def slide_end(prs, spec, meta, D, figs):
    return [add_slide(prs, LAYOUT["end"])]


KINDS = {"cover": slide_cover, "divider": slide_divider, "figure": slide_figure, "cards": slide_cards,
         "stats": slide_stats, "columns": slide_columns, "end": slide_end}


# ---------------------------------------------------------------------------------------------
def _sha1(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def keep_hand_edits(out: Path) -> None:
    """A deck that was changed after the last build was edited in PowerPoint: copy it to hand_edited/ before it is
    overwritten. out/<name>.built holds the hash of the file the last build wrote."""
    if not out.exists():
        return
    digest, stamp = _sha1(out), out.with_suffix(".built")
    if stamp.exists() and stamp.read_text().strip() == digest:
        return
    keep = ROOT / "hand_edited"
    keep.mkdir(exist_ok=True)
    if any(_sha1(f) == digest for f in keep.glob("*.pptx")):
        return
    name = f"{out.stem}_{datetime.fromtimestamp(out.stat().st_mtime):%Y-%m-%d_%H%M}_hand_edited.pptx"
    shutil.copy2(out, keep / name)
    print(f"  note: {out.name} was changed since the last build (edited in PowerPoint?). Kept a copy: hand_edited/{name}")


def build(render: bool = True) -> Path:
    D = F.load_data()
    deck = yaml.safe_load((HERE / "deck.yaml").read_text())
    meta = deck["meta"]
    if render:
        figs = render_figures.render(data=D)
    else:
        figs = {p.stem: p for p in (OUT / "figures").glob("*.png")}

    prs = open_template(TEMPLATE)
    count = 0
    for spec in deck["slides"]:
        if spec.get("skip"):
            continue
        if spec["kind"] not in KINDS:
            raise KeyError(f"slide {spec.get('id')!r}: unknown kind {spec['kind']!r}")
        for slide in KINDS[spec["kind"]](prs, spec, meta, D, figs):
            drop_empty_placeholders(slide)
            if spec["kind"] not in ("cover", "end"):
                add_slide_number(slide)
            if spec.get("notes"):
                slide.notes_slide.notes_text_frame.text = fill(spec["notes"], D).strip()
            count += 1

    prs.core_properties.title = fill(meta["title"], D)
    prs.core_properties.author = "Aaron Odell"
    OUT.mkdir(exist_ok=True)
    out = OUT / f"{meta.get('out', 'deck')}.pptx"
    keep_hand_edits(out)
    prs.save(str(out))
    out.with_suffix(".built").write_text(_sha1(out))
    print(f"{count} slides -> {out}")
    return out


def check_inputs() -> None:
    """The UNC template and the preprint PDF are not kept in git: say so when one is missing."""
    missing = [p.name for p in (TEMPLATE, paper.PDF) if not p.exists()]
    if missing:
        sys.exit(f"Missing from assets/: {', '.join(missing)}. These are not kept in git; see \"Two files to add first\" in the README.")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="?", default="deck", choices=["deck", "figures", "crops"])
    ap.add_argument("--no-figures", action="store_true", help="reuse the figures already in out/figures/")
    args = ap.parse_args(argv)
    check_inputs()
    if args.what == "crops":
        for p in paper.crop_all(OUT / "crops"):
            print(p)
    elif args.what == "figures":
        for key, path in render_figures.render().items():
            print(f"{key:24s} {path}")
    else:
        build(render=not args.no_figures)


if __name__ == "__main__":
    sys.exit(main())
