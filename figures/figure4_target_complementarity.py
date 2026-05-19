#!/usr/bin/env python3
"""
Figure 4 composite (v15) — Synergistic miRNA pairs selectively reinforce
on-target pathways through complementarity.

Three panels (single row):
  (A) Boxplot — incremental coverage by pathway category
      (On-target vs Liability vs Housekeeping)
  (B) Scatter — on-target vs liability incremental coverage per pair
  (C) Scatter — on-target vs housekeeping incremental coverage per pair

This script is a thin wrapper around the existing
scripts/figures/create_final_figure_v14.py — it reuses the same plotting
logic (load_batch_data, categorize_modules, create_figure) and routes the
output to the canonical manuscript figure name.

Output: figures/Figure 4 v15.{png,pdf}
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_FIGURES = ROOT / "figures"

sys.path.insert(0, str(SCRIPTS_FIGURES))
from figure4_panels import (  # noqa: E402
    create_figure,
    load_batch_data,
)

BATCH_DIR = ROOT / "03_target_analysis" / "outputs"
OUT_DIR = ROOT / "figures"


def main() -> None:
    metrics = load_batch_data(BATCH_DIR)
    if metrics.empty:
        raise SystemExit(f"No metrics found in {BATCH_DIR}")

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Figure 4 v15.png"
    out_pdf = OUT_DIR / "Figure 4 v15.pdf"

    # create_figure handles a single output path. Render twice so each
    # format gets its own bbox_inches='tight' pass without re-running the
    # data load.
    create_figure(metrics, str(out_png))
    create_figure(metrics, str(out_pdf))


if __name__ == "__main__":
    main()
