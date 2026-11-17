"""Small drawing layer over matplotlib for the talk's diagrams.

A Canvas is measured in inches with the origin at the top left, the same as a slide, and is
placed on the slide at 1:1. So a 14 pt label here is 14 pt on the slide, and a shape 2 in wide
here is 2 in wide there.
"""
from __future__ import annotations

import re
import math
from dataclasses import dataclass

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Polygon, Rectangle, Wedge  # noqa: E402
from matplotlib.path import Path  # noqa: E402

from . import theme as T  # noqa: E402

plt.rcParams.update({
    "font.family": T.FONT,
    "mathtext.fontset": "custom",
    "mathtext.rm": T.FONT,
    "mathtext.it": f"{T.FONT}:italic",
    "mathtext.bf": f"{T.FONT}:bold",
    "svg.fonttype": "none",
})


def sci(p: float, digits: int = 1, bold: bool = False) -> str:
    """1.3e-18 -> '1.3 × 10⁻¹⁸' (mathtext superscript, in the text font). bold: for use in bold text."""
    if p == 0:
        return "0"
    exp = math.floor(math.log10(abs(p)))
    mant = p / 10 ** exp
    if round(mant, digits) >= 10:
        mant, exp = mant / 10, exp + 1
    mant_s = f"{mant:.{digits}f}"
    face = "mathbf" if bold else "mathregular"
    return f"{mant_s} × 10$^{{\\{face}{{{exp}}}}}$"


def sci_plain(p: float, digits: int = 1) -> str:
    """1.3e-18 -> '1.3 × 10^-18' with unicode superscripts, for slide text (not mathtext)."""
    exp = math.floor(math.log10(abs(p)))
    mant = p / 10 ** exp
    if round(mant, digits) >= 10:
        mant, exp = mant / 10, exp + 1
    sup = str(exp).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{mant:.{digits}f} × 10{sup}"


@dataclass
class Box:
    """Where a drawn shape ended up, for hanging connectors off it."""
    x: float
    y: float
    w: float
    h: float

    @property
    def cx(self): return self.x + self.w / 2
    @property
    def cy(self): return self.y + self.h / 2
    @property
    def l(self): return (self.x, self.cy)
    @property
    def r(self): return (self.x + self.w, self.cy)
    @property
    def t(self): return (self.cx, self.y)
    @property
    def b(self): return (self.cx, self.y + self.h)
    @property
    def right(self): return self.x + self.w
    @property
    def bottom(self): return self.y + self.h


class Canvas:
    def __init__(self, w: float = T.FIG_W, h: float = T.FIG_H, dpi: int = T.FIG_DPI, bg: str | None = T.WHITE):
        self.w, self.h, self.dpi = w, h, dpi
        self.fig = plt.figure(figsize=(w, h), dpi=dpi)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, w)
        self.ax.set_ylim(h, 0)
        self.ax.axis("off")
        if bg:
            self.fig.patch.set_facecolor(bg)
        else:
            self.fig.patch.set_alpha(0)
        self._renderer = self.fig.canvas.get_renderer()

    # --- text ---------------------------------------------------------------------------------
    def text(self, x, y, s, size=13, color=T.INK, weight="normal", ha="left", va="center", style="normal",
             alpha=1.0, rotation=0, zorder=6, linespacing=1.18, **kw):
        return self.ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, fontstyle=style,
                            alpha=alpha, rotation=rotation, zorder=zorder, linespacing=linespacing, **kw)

    def measure(self, s, size=13, weight="normal", style="normal"):
        """Width and height of a text string, in inches."""
        t = self.ax.text(0, 0, s, fontsize=size, fontweight=weight, fontstyle=style)
        bb = t.get_window_extent(self._renderer)
        t.remove()
        return bb.width / self.dpi, bb.height / self.dpi

    def fit(self, s, max_w, size, weight="normal", min_size=9.0, style="normal"):
        """Largest size <= `size` at which `s` is no wider than `max_w`."""
        while size > min_size and self.measure(s, size, weight, style)[0] > max_w:
            size -= 0.5
        return size

    def wrap(self, s, max_w, size, weight="normal"):
        """Greedy word wrap to `max_w` inches. Returns the string with newlines.
        Breaks at spaces only, so a no-break space (\u00a0) keeps two words together."""
        lines, cur = [], ""
        for word in re.split(r"[ \t\n]+", s.strip()):
            trial = (cur + " " + word).strip()
            if cur and self.measure(trial, size, weight)[0] > max_w:
                lines.append(cur)
                cur = word
            else:
                cur = trial
        if cur:
            lines.append(cur)
        return "\n".join(lines)

    def rich(self, x, y, runs, max_w, size=13, color=T.INK, line_h=None, zorder=6):
        """Word-wrapped text with styled runs, e.g. a bold coloured phrase inside a sentence.

        `runs` is [(text, {"color": ..., "weight": ...}), ...]; runs join at spaces. The block is
        centred vertically on `y`. Returns the number of lines.
        """
        line_h = line_h or size / 72 * 1.32
        space = self.measure("a a", size)[0] - self.measure("aa", size)[0]
        words = [(w, st or {}) for text, st in runs for w in text.split()]
        lines, cur, cur_w = [], [], 0.0
        for w, st in words:
            ww = self.measure(w, size, st.get("weight", "normal"))[0]
            if cur and cur_w + space + ww > max_w:
                lines.append(cur)
                cur, cur_w = [], 0.0
            cur_w += (space if cur else 0.0) + ww
            cur.append((w, st, ww))
        if cur:
            lines.append(cur)
        base = y - (len(lines) - 1) * line_h / 2 + size / 72 * 0.36      # baseline that centres the block on y
        for i, line in enumerate(lines):
            cx = x
            for w, st, ww in line:
                self.text(cx, base + i * line_h, w, size=size, color=st.get("color", color),
                          weight=st.get("weight", "normal"), va="baseline", zorder=zorder)
                cx += ww + space
        return len(lines)

    # --- shapes -------------------------------------------------------------------------------
    def box(self, x, y, w, h, fill=T.WHITE, edge=T.NAVY, lw=1.4, r=0.12, zorder=2, alpha=1.0, ls="-"):
        r = min(r, w / 2, h / 2)
        patch = FancyBboxPatch((x + r, y + r), w - 2 * r, h - 2 * r, boxstyle=f"round,pad={r},rounding_size={r}",
                               facecolor=fill if fill else "none", edgecolor=edge if edge else "none",
                               linewidth=lw if edge else 0, zorder=zorder, alpha=alpha, linestyle=ls, mutation_aspect=1)
        self.ax.add_patch(patch)
        return Box(x, y, w, h)

    def panel(self, x, y, w, h, tag=None, title=None, fill=T.TINT, edge=None, tag_color=T.NAVY, title_color=T.MUTED):
        """Soft background panel with an optional 'a)  Title' header."""
        b = self.box(x, y, w, h, fill=fill, edge=edge, lw=1.0, r=0.14, zorder=0)
        if tag or title:
            tx = x + 0.2
            if tag:
                self.text(tx, y + 0.28, tag, size=12.5, weight="bold", color=tag_color)
                tx += self.measure(tag, 12.5, "bold")[0] + 0.12
            if title:
                self.text(tx, y + 0.28, title, size=12.5, weight="bold", color=title_color)
        return b

    def pill(self, cx, cy, label, w=None, h=0.46, size=13, fill=T.WHITE, edge=T.NAVY, lw=1.5, color=T.INK,
             weight="normal", sub=None, sub_size=None, sub_color=None, sub_style="normal", pad=0.24, zorder=3,
             style="normal", r=None):
        """Rounded label centred on (cx, cy). Width fits the text unless given. `sub` adds a second line."""
        sub_size = sub_size or size - 2
        tw = self.measure(label, size, weight, style)[0]
        if sub:
            tw = max(tw, self.measure(sub, sub_size, "normal", sub_style)[0])
            h = max(h, 0.78)
        w = w or tw + 2 * pad
        b = self.box(cx - w / 2, cy - h / 2, w, h, fill=fill, edge=edge, lw=lw, r=h / 2 if r is None else r, zorder=zorder)
        if sub:
            self.text(cx, cy - 0.15, label, size=size, color=color, weight=weight, ha="center", style=style, zorder=zorder + 1)
            self.text(cx, cy + 0.17, sub, size=sub_size, color=sub_color or color, ha="center", style=sub_style, zorder=zorder + 1)
        else:
            self.text(cx, cy, label, size=size, color=color, weight=weight, ha="center", style=style, zorder=zorder + 1)
        return b

    def pill_left(self, x, cy, label, **kw):
        """Pill whose left edge is at x."""
        size, weight, pad = kw.get("size", 13), kw.get("weight", "normal"), kw.get("pad", 0.24)
        w = kw.pop("w", None) or self.measure(label, size, weight, kw.get("style", "normal"))[0] + 2 * pad
        return self.pill(x + w / 2, cy, label, w=w, **kw)

    def circle(self, x, y, r, fill=T.WHITE, edge=T.NAVY, lw=1.2, zorder=3, alpha=1.0, ls="-"):
        self.ax.add_patch(Circle((x, y), r, facecolor=fill if fill else "none", edgecolor=edge if edge else "none",
                                 linewidth=lw if edge else 0, zorder=zorder, alpha=alpha, linestyle=ls))
        return Box(x - r, y - r, 2 * r, 2 * r)

    def square(self, x, y, s, fill=T.WHITE, edge=T.NAVY, lw=1.2, zorder=3, r=0.03, alpha=1.0):
        return self.box(x - s / 2, y - s / 2, s, s, fill=fill, edge=edge, lw=lw, r=r, zorder=zorder, alpha=alpha)

    def rect(self, x, y, w, h, fill=T.NAVY, edge=None, lw=0, zorder=2, alpha=1.0):
        self.ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge if edge else "none", linewidth=lw,
                                    zorder=zorder, alpha=alpha))
        return Box(x, y, w, h)

    # --- connectors ---------------------------------------------------------------------------
    def line(self, x0, y0, x1, y1, color=T.FAINT, lw=1.2, ls="-", zorder=1, alpha=1.0):
        self.ax.plot([x0, x1], [y0, y1], color=color, lw=lw, ls=ls, zorder=zorder, alpha=alpha, solid_capstyle="round")

    def arrow(self, x0, y0, x1, y1, color=T.FAINT, lw=1.5, head=0.11, zorder=4, ls="-"):
        ang = math.atan2(y1 - y0, x1 - x0)
        bx, by = x1 - head * 0.9 * math.cos(ang), y1 - head * 0.9 * math.sin(ang)
        self.ax.plot([x0, bx], [y0, by], color=color, lw=lw, ls=ls, zorder=zorder, solid_capstyle="butt")
        self._head(x1, y1, ang, head, color, zorder)

    def _head(self, x, y, ang, size, color, zorder):
        a, b = ang + math.radians(152), ang - math.radians(152)
        pts = [(x, y), (x + size * math.cos(a), y + size * math.sin(a)), (x + size * math.cos(b), y + size * math.sin(b))]
        self.ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=color, linewidth=0.5, zorder=zorder))

    def scurve(self, x0, y0, x1, y1, color=T.NAVY, lw=1.4, ls="-", zorder=1, alpha=1.0, head=0.0, bend=0.5):
        """Sideways S: leaves (x0, y0) and arrives at (x1, y1) horizontally."""
        dx = (x1 - x0) * bend
        path = Path([(x0, y0), (x0 + dx, y0), (x1 - dx, y1), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        self.ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, ls=ls, zorder=zorder, alpha=alpha,
                                    capstyle="round"))
        if head:
            self._head(x1, y1, 0 if x1 > x0 else math.pi, head, color, zorder)

    def vcurve(self, x0, y0, x1, y1, color=T.NAVY, lw=1.4, ls="-", zorder=1, alpha=1.0, head=0.0, bend=0.5):
        """Upright S: leaves and arrives vertically."""
        dy = (y1 - y0) * bend
        path = Path([(x0, y0), (x0, y0 + dy), (x1, y1 - dy), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        self.ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, ls=ls, zorder=zorder, alpha=alpha,
                                    capstyle="round"))
        if head:
            self._head(x1, y1, math.pi / 2 if y1 > y0 else -math.pi / 2, head, color, zorder)

    def ccurve(self, x0, y0, x1, y1, bulge=0.5, color=T.NAVY, lw=1.4, ls="-", zorder=1, alpha=1.0, head=0.0):
        """C shape: leaves (x0, y0) sideways, swings out and comes back to (x1, y1) from the same side.
        bulge > 0 swings to the right, < 0 to the left; the curve reaches about three quarters of it."""
        bx = (max(x0, x1) if bulge > 0 else min(x0, x1)) + bulge
        path = Path([(x0, y0), (bx, y0), (bx, y1), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        self.ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, ls=ls, zorder=zorder, alpha=alpha,
                                    capstyle="round"))
        if head:
            self._head(x1, y1, math.pi if bulge > 0 else 0, head, color, zorder)

    def arc(self, x0, y0, x1, y1, lift=0.6, color=T.FAINT, lw=1.2, ls="-", zorder=1, alpha=1.0):
        """Arch from one point to another, rising `lift` inches above the higher of the two."""
        top = min(y0, y1) - lift
        path = Path([(x0, y0), (x0, top), (x1, top), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        self.ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, ls=ls, zorder=zorder, alpha=alpha,
                                    capstyle="round"))

    # --- glyphs -------------------------------------------------------------------------------
    def badge(self, x, y, code, color, r=0.19, size=None, text_color=T.WHITE, zorder=5):
        """Species badge: a filled circle with the two-letter code (Mm, Dr, Ce...)."""
        self.circle(x, y, r, fill=color, edge=None, zorder=zorder)
        self.text(x, y + r * 0.04, code, size=size or r * 46, color=text_color, weight="bold", ha="center", zorder=zorder + 1)
        return Box(x - r, y - r, 2 * r, 2 * r)

    def person(self, cx, top, h, color=T.NAVY, alpha=1.0, zorder=3):
        """Plain figure: a head over rounded shoulders. `top` is the top of the head; half-width is 0.18 h."""
        head_r = h * 0.15
        self.ax.add_patch(Circle((cx, top + head_r), head_r, facecolor=color, edgecolor="none", alpha=alpha, zorder=zorder))
        bw = h * 0.36
        by, bottom = top + 2 * head_r + h * 0.05, top + h
        l, r, sh, rr = cx - bw / 2, cx + bw / 2, bw / 2, bw * 0.14
        k = 0.5523 * sh
        verts = [(l, bottom - rr), (l, by + sh), (l, by + sh - k), (cx - k, by), (cx, by), (cx + k, by), (r, by + sh - k),
                 (r, by + sh), (r, bottom - rr), (r, bottom - rr * 0.45), (r - rr * 0.55, bottom), (r - rr, bottom),
                 (l + rr, bottom), (l + rr * 0.55, bottom), (l, bottom - rr * 0.45), (l, bottom - rr), (l, bottom - rr)]
        codes = [Path.MOVETO, Path.LINETO] + [Path.CURVE4] * 6 + [Path.LINETO] + [Path.CURVE4] * 3 + [Path.LINETO] + \
                [Path.CURVE4] * 3 + [Path.CLOSEPOLY]
        self.ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=zorder))
        return Box(cx - bw / 2, top, bw, h)

    def wedge(self, x, y, r, a0, a1, fill, zorder=3):
        self.ax.add_patch(Wedge((x, y), r, a0, a1, facecolor=fill, edgecolor="none", zorder=zorder))

    # --- output -------------------------------------------------------------------------------
    def image(self, path, x, y, w=None, h=None, zorder=2):
        """Place a raster image with its top left at (x, y). Give w or h; the other follows the aspect.
        The image is resampled once, to the canvas's own pixels, so its edges stay clean."""
        from PIL import Image
        img = Image.open(path).convert("RGB")
        iw, ih = img.size
        if w is None:
            w = h * iw / ih
        if h is None:
            h = w * ih / iw
        px = (max(1, round(w * self.dpi)), max(1, round(h * self.dpi)))
        if px != img.size:
            img = img.resize(px, Image.LANCZOS)
        self.ax.imshow(np.asarray(img), extent=(x, x + w, y + h, y), zorder=zorder, interpolation="none", aspect="auto")
        self.ax.set_xlim(0, self.w)
        self.ax.set_ylim(self.h, 0)
        return Box(x, y, w, h)

    def save(self, path):
        # opaque by default: a half-covered transparent edge pixel prints as a grey hairline in PowerPoint's PDF export
        self.fig.savefig(path, dpi=self.dpi, transparent=self.fig.patch.get_alpha() == 0, facecolor=self.fig.get_facecolor())
        plt.close(self.fig)
        return path
