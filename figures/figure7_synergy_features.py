#!/usr/bin/env python3
"""
Figure 7 composite (v15) — Features distinguishing synergistic from
non-synergistic miRNA pairs.

Four-panel boxplot comparison (synergistic dual-positive pairs vs
non-synergistic combinations) across four features:
  (A) Target overlap (Jaccard)
  (B) Combined target set size
  (C) Best single-agent neurite length
  (D) Expression correlation in tumors

Each panel: Mann-Whitney U test (two-sided), significance annotation
above the boxes, jittered individual points overlaid.

Source data: data/target_analysis/synergy_features/all_features.csv
(produced by scripts/target_analysis/synergy_features_v14.py).

Plotting logic follows the in-place regenerate_panel.py reference;
this script adds a PDF output and uses the canonical manuscript figure
filename.

Output: figures/Figure 7 v15.{png,pdf}
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "03_target_analysis" / "outputs" / "all_features.csv"
OUT_DIR = ROOT / "figures"

PANELS = [
    ("jaccard", "(A) Target overlap (Jaccard)"),
    ("size_union", "(B) Combined target set size"),
    ("best_single_nl", "(C) Best single-agent neurite length"),
    ("expr_corr", "(D) Expression correlation in tumors"),
]

COLORS = {"syn": "#FF8C00", "non": "#CCCCCC"}


def fmt_p(pval: float) -> tuple[str, str]:
    if pval < 0.001:
        sig = "***"
    elif pval < 0.01:
        sig = "**"
    elif pval < 0.05:
        sig = "*"
    else:
        sig = "ns"
    text = f"p={pval:.2e} ({sig})" if pval < 0.001 else f"p={pval:.3f} ({sig})"
    return text, sig


def plot_panel(ax, df: pd.DataFrame, column: str, title: str) -> None:
    sub = df[[column, "dual_positive"]].dropna()
    syn = sub.loc[sub["dual_positive"], column]
    non = sub.loc[~sub["dual_positive"], column]

    if len(syn) > 0 and len(non) > 0:
        _, pval = mannwhitneyu(syn, non, alternative="two-sided")
        p_text, _ = fmt_p(pval)
    else:
        p_text = "N/A"

    bp = ax.boxplot(
        [non, syn],
        tick_labels=["Non-syn", "Synergistic"],
        widths=0.5,
        patch_artist=True,
        showfliers=False,
    )
    bp["boxes"][0].set_facecolor(COLORS["non"])
    bp["boxes"][0].set_alpha(0.6)
    bp["boxes"][1].set_facecolor(COLORS["syn"])
    bp["boxes"][1].set_alpha(0.6)
    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.5)

    rng = np.random.default_rng(42)
    for i, (data, color) in enumerate([(non, COLORS["non"]), (syn, COLORS["syn"])]):
        jitter = rng.normal(0, 0.04, size=len(data))
        ax.scatter(
            np.full(len(data), i + 1) + jitter,
            data.values,
            color=color,
            edgecolors="black",
            linewidths=0.3,
            alpha=0.6,
            s=20,
            zorder=3,
        )

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.annotate(
        p_text,
        xy=(0.5, 0.97),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="gray",
            alpha=0.8,
        ),
    )
    ax.set_xlabel(f"n={len(non)}       n={len(syn)}", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    df = pd.read_csv(DATA, index_col=0)
    fig, axes = plt.subplots(1, 4, figsize=(16, 5))
    for ax, (col, title) in zip(axes, PANELS):
        plot_panel(ax, df, col, title)
    plt.tight_layout()

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Figure 7 v15.png"
    out_pdf = OUT_DIR / "Figure 7 v15.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
