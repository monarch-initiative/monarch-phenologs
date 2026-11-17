#!/bin/bash
# Export a .pptx to PDF with PowerPoint (macOS), then cut one PNG per slide for a quick look.
#   tools/export_pdf.sh out/phenologs_talk.pptx
# PowerPoint is sandboxed: it is handed the file with `open`, and writes the PDF inside its own
# container, from where the PDF is moved next to the deck. Needs pdftoppm for the PNGs.
set -euo pipefail
deck="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
pdf="${deck%.pptx}.pdf"
name="$(basename "$deck")"
box="$HOME/Library/Group Containers/UBF8T346G9.Office"
tmp="$box/phenologs_talk_export_$$.pdf"

open -g -a "Microsoft PowerPoint" "$deck"
osascript <<OSA
with timeout of 300 seconds
  tell application "Microsoft PowerPoint"
    repeat 60 times
      if (exists presentation "$name") then exit repeat
      delay 0.5
    end repeat
    set p to presentation "$name"
    save p in (POSIX file "$tmp") as save as PDF
    close p saving no
  end tell
end timeout
OSA
mv "$tmp" "$pdf"
echo "$pdf"

if command -v pdftoppm >/dev/null; then
  shots="$(dirname "$deck")/slides"
  rm -rf "$shots" && mkdir -p "$shots"
  pdftoppm -png -r "${SLIDE_DPI:-80}" "$pdf" "$shots/slide"
  echo "$shots/"
fi
