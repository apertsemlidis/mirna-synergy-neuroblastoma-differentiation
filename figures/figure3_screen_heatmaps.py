#!/usr/bin/env python3
"""
Figure 3 composite (v15) — Pairwise screen heatmaps.

Two panels, each a 40x40 miRNA-x-miRNA symmetric square split along the
diagonal:
  (A) Neurite length        — upper-right: HSA-relative; lower-left: absolute
  (B) Cell body cluster area — same triangle convention

Superhits are outlined (thin border) and the six combinations advanced to
follow-on assays are outlined with thick borders.

Inputs (all under data/screen/):
  HSA_dfs/nl_slice(96, 126, None).csv
  ABS_dfs/nl_slice(96, 126, None).csv
  HSA_dfs/cbca_slice(96, 126, None).csv
  ABS_dfs/cbca_slice(96, 126, None).csv
  superhits.csv

Output: figures/Figure 3 v15.{png,pdf}

Port of the heatmap block in scripts/screen/screen_analysis_v14.py
(lines 599-685); refactored into a single composite figure with self-
contained data loads and no implicit notebook globals.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "01_screen"
OUT_DIR = ROOT / "figures"

HITS_TO_ADVANCE = {
    "hsa-miR-124-3p + hsa-miR-363-3p",
    "hsa-miR-124-3p + hsa-miR-34b-5p",
    "hsa-miR-137 + hsa-miR-450b-3p",
    "hsa-miR-137 + hsa-miR-449b-5p",
    "hsa-miR-137 + hsa-miR-17-5p",
    "hsa-miR-19b-3p + hsa-miR-2110",
}




def trim_hsa(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = df.columns.map(lambda x: x[4:] if x.startswith("hsa-") else x)
    df.index = df.index.map(lambda x: x[4:] if x.startswith("hsa-") else x)
    return df


def draw_panel(
    ax: plt.Axes,
    abs_df: pd.DataFrame,
    hsa_df: pd.DataFrame,
    superhits: pd.Index,
    absolute_label: str,
    hsa_label: str,
    hsa_cbar_label: str,
) -> None:
    mask_lower = np.tril(np.ones_like(abs_df, dtype=bool), k=0)
    mask_upper = np.triu(np.ones_like(hsa_df, dtype=bool), k=0)

    sns.heatmap(
        hsa_df,
        ax=ax,
        square=True,
        linewidths=0,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={
            "shrink": 0.45,
            "use_gridspec": False,
            "location": "right",
            "pad": 0.02,
            "label": hsa_cbar_label,
        },
        mask=mask_lower,
        center=0,
        cmap="RdGy_r",
    )
    sns.heatmap(
        abs_df,
        ax=ax,
        square=True,
        linewidths=0,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={
            "shrink": 0.45,
            "use_gridspec": False,
            "location": "bottom",
            "orientation": "horizontal",
            "pad": 0.04,
            "label": absolute_label,
        },
        mask=mask_upper,
        cmap="viridis",
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(
        0.02,
        0.98,
        hsa_label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )

    miR_to_idx = {miR: i for i, miR in enumerate(abs_df.columns)}
    for combo in superhits:
        parts = combo.split(" + ")
        if len(parts) != 2:
            continue
        a, b = parts[0][4:], parts[1][4:]
        if a not in miR_to_idx or b not in miR_to_idx:
            continue
        x = miR_to_idx[a]
        y = miR_to_idx[b]
        thick = combo in HITS_TO_ADVANCE
        lw = 2.5 if thick else 0.9
        if x < y:
            ax.add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor="white", lw=lw))
            ax.add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor="black", lw=lw))
        else:
            ax.add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor="black", lw=lw))
            ax.add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor="white", lw=lw))


def main() -> None:
    superhits_csv = SCREEN / "superhits.csv"
    superhits = pd.read_csv(superhits_csv, index_col=0).index

    panels = [
        (
            "A",
            "nl",
            "Neurite length",
            "neurite length (absolute units)",
            "neurite length / highest single agent",
        ),
        (
            "B",
            "cbca",
            "Cell body cluster area",
            "cell body cluster area (absolute units)",
            "cell body cluster area / lowest single agent",
        ),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(22, 11), constrained_layout=True)

    for ax, (letter, ini, title, abs_lbl, hsa_lbl) in zip(axes, panels):
        abs_df = pd.read_csv(SCREEN / "heatmap_data" / f"{ini}_absolute.csv", index_col=0)
        hsa_df = pd.read_csv(SCREEN / "heatmap_data" / f"{ini}_hsa.csv", index_col=0)
        abs_df = trim_hsa(abs_df)
        hsa_df = trim_hsa(hsa_df)

        draw_panel(
            ax,
            abs_df,
            hsa_df,
            superhits,
            absolute_label=abs_lbl,
            hsa_label=title,
            hsa_cbar_label=hsa_lbl,
        )
        ax.text(
            -0.05,
            1.04,
            f"({letter})",
            transform=ax.transAxes,
            fontsize=18,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Figure 3 v15.png"
    out_pdf = OUT_DIR / "Figure 3 v15.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
