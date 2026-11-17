"""Tile out/slides/*.png into contact sheets for a quick visual check.  python tools/contact_sheet.py [cols] [per_sheet]"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw

root = Path(__file__).resolve().parents[1] / "out"
cols = int(sys.argv[1]) if len(sys.argv) > 1 else 2
per = int(sys.argv[2]) if len(sys.argv) > 2 else 6
files = sorted((root / "slides").glob("slide-*.png"))
(root / "sheets").mkdir(exist_ok=True)
for s in range(0, len(files), per):
    ims = [Image.open(f).convert("RGB") for f in files[s:s + per]]
    w, h = ims[0].size
    rows = -(-len(ims) // cols)
    sheet = Image.new("RGB", (cols * (w + 8), rows * (h + 26)), "#888888")
    d = ImageDraw.Draw(sheet)
    for i, im in enumerate(ims):
        x, y = (i % cols) * (w + 8), (i // cols) * (h + 26)
        sheet.paste(im, (x, y + 22))
        d.text((x + 6, y + 5), f"slide {s + i + 1}", fill="white")
    out = root / "sheets" / f"sheet-{s // per + 1:02d}.png"
    sheet.save(out)
    print(out, sheet.size)
