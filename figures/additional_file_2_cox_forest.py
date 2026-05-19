#!/usr/bin/env python3
"""
Additional file 2 (v15) — Combined Cox proportional hazards forest plot.

Single-panel figure: one row per synergistic miRNA pair, showing
MYCN-adjusted, age-stratified HR for coordinated "both high" expression.

Source: survival/cox_forest/cox_forest_combined_v14.{png,svg}
(produced by survival/cox_forest/cox_forest_combined.py).

This thin wrapper re-saves the source PNG under the standardized
additional-file name. No re-rendering or modification.
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "04_survival"
SRC_SCRIPT = SRC_DIR / "cox_forest_combined.py"
SRC_PNG = SRC_DIR / "cox_forest_combined.png"
SRC_SVG = SRC_DIR / "cox_forest_combined.svg"
SRC_PDF = SRC_DIR / "cox_forest_combined.pdf"
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    if SRC_SCRIPT.exists() and not (SRC_PNG.exists() and SRC_PDF.exists()):
        subprocess.run([sys.executable, str(SRC_SCRIPT)], cwd=SRC_DIR, check=True)
    assert SRC_PNG.exists(), f"Missing source: {SRC_PNG}"
    for src, dst_name in [
        (SRC_PNG, "Additional file 2 v15.png"),
        (SRC_SVG, "Additional file 2 v15.svg"),
        (SRC_PDF, "Additional file 2 v15.pdf"),
    ]:
        if src.exists():
            dst = OUT_DIR / dst_name
            shutil.copy2(src, dst)
            print(f"Saved: {dst}")


if __name__ == "__main__":
    main()
