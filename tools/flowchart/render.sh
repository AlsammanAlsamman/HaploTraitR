#!/usr/bin/env bash
# Render the HaploTraitR workflow diagram (flowchart.html) to
# man/figures/haplotraitr_workflow.png with headless Chrome.
# Usage: tools/flowchart/render.sh   (run from the package root)
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
html="file://$here/flowchart.html"
out="$here/../../man/figures/haplotraitr_workflow.png"
chrome="${CHROME:-$(command -v google-chrome || command -v chromium || command -v chromium-browser)}"
flags=(--headless --disable-gpu --no-sandbox --hide-scrollbars --allow-file-access-from-files)

# First pass: let the script lay out the diagram and report its height
height=$("$chrome" "${flags[@]}" --window-size=1400,4000 --virtual-time-budget=3000 --dump-dom "$html" 2>/dev/null |
         grep -o 'data-height="[0-9]*"' | grep -o '[0-9]*')

# Second pass: screenshot at 2x resolution with some spare room (the headless
# viewport is slightly smaller than the window), then crop to the diagram height
"$chrome" "${flags[@]}" --window-size="1400,$((height + 200))" --force-device-scale-factor=2 \
  --virtual-time-budget=3000 --screenshot="$out" "$html" >/dev/null 2>&1
python3 -c "
from PIL import Image
im = Image.open('$out')
im.crop((0, 0, im.width, $height * 2)).save('$out', optimize=True)
"
echo "Wrote $out (1400 x $height, 2x)"
