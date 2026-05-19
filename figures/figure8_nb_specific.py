#!/usr/bin/env python3
"""
Figure 8 composite (v15) — Incremental coverage of NB-specific
transcriptional modules by synergistic miRNA pairs.

Single-axis boxplot showing incremental coverage (defined as
coverage(A∪B) − max[coverage(A), coverage(B)]) for four NB-specific
gene modules:
  - Adrenergic (ADRN) — neuroblastoma differentiation identity
  - Retinoid response (GO) — differentiation pathway
  - MYCN targets (Wei et al.) — proliferation-driving program
  - Mesenchymal (MES) — undifferentiated cell state

Statistical brackets show ADRN vs MYCN targets and ADRN vs Retinoid
response (one-sided Mann-Whitney U).

Source data:
  data/target_analysis/batch_ts_nb_specific/all_pairs_nb_metrics.csv
(produced by scripts/target_analysis/nb_specific_analysis_v14.py).

Plotting logic follows the in-place regenerate_boxplot.py reference;
this script adds a PDF output and uses the canonical manuscript figure
filename.

Output: figures/Figure 8 v15.{png,pdf}
"""

import warnings
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "03_target_analysis" / "outputs" / "all_pairs_nb_metrics.csv"
OUT_DIR = ROOT / "figures"

CATEGORY_ORDER = ["nb_differentiation", "mycn", "nb_undifferentiated"]
MODULE_DISPLAY_ORDER = [
    "Adrenergic (ADRN)",
    "Retinoid response (GO)",
    "MYCN targets (Wei et al.)",
    "Mesenchymal (MES)",
]

CAT_COLORS = {
    "nb_differentiation": "#2ca02c",
    "mycn": "#d62728",
    "nb_undifferentiated": "#7f7f7f",
}


def fmt_p(p: float) -> str:
    if p < 0.001:
        return "p < 0.001"
    elif p < 0.01:
        return f"p = {p:.3f}"
    else:
        return f"p = {p:.2f}"


def add_bracket(ax, x1: float, x2: float, y: float, h: float, p_text: str) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color="black")
    ax.text(
        (x1 + x2) / 2,
        y + h + 0.002,
        p_text,
        ha="center",
        va="bottom",
        fontsize=8.5,
    )


def main() -> None:
    df = pd.read_csv(DATA)
    df["category"] = pd.Categorical(
        df["category"], categories=CATEGORY_ORDER, ordered=True
    )
    df["module_display"] = pd.Categorical(
        df["module_display"], categories=MODULE_DISPLAY_ORDER, ordered=True
    )
    df = df.sort_values(["category", "module_display"])

    mod_to_cat = df.drop_duplicates("module_display").set_index("module_display")[
        "category"
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    positions = list(range(len(MODULE_DISPLAY_ORDER)))
    box_data = [
        df.loc[df["module_display"] == m, "incremental_coverage"].values
        for m in MODULE_DISPLAY_ORDER
    ]

    bp = ax.boxplot(
        box_data,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        zorder=2,
    )
    for i, (patch, med) in enumerate(zip(bp["boxes"], bp["medians"])):
        cat = mod_to_cat[MODULE_DISPLAY_ORDER[i]]
        color = CAT_COLORS[cat]
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
        patch.set_edgecolor(color)
        med.set_color("black")
        med.set_linewidth(1.5)
    for element in ["whiskers", "caps"]:
        for line in bp[element]:
            line.set_color("gray")
            line.set_linewidth(1)

    rng = np.random.default_rng(42)
    for i, m in enumerate(MODULE_DISPLAY_ORDER):
        vals = df.loc[df["module_display"] == m, "incremental_coverage"].values
        cat = mod_to_cat[m]
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(
            np.full_like(vals, i, dtype=float) + jitter,
            vals,
            color=CAT_COLORS[cat],
            alpha=0.6,
            s=18,
            edgecolors="white",
            linewidths=0.3,
            zorder=3,
        )

    adrn = df.loc[
        df["module_display"] == "Adrenergic (ADRN)", "incremental_coverage"
    ].values
    mycn = df.loc[
        df["module_display"] == "MYCN targets (Wei et al.)", "incremental_coverage"
    ].values
    retinoid = df.loc[
        df["module_display"] == "Retinoid response (GO)", "incremental_coverage"
    ].values

    _, p_adrn_mycn = mannwhitneyu(adrn, mycn, alternative="greater")
    _, p_adrn_ret = mannwhitneyu(adrn, retinoid, alternative="greater")

    ymax = max(np.max(d) for d in box_data if len(d) > 0)
    add_bracket(ax, 0, 2, ymax + 0.015, 0.008, fmt_p(p_adrn_mycn))
    add_bracket(ax, 0, 1, ymax + 0.055, 0.008, fmt_p(p_adrn_ret))

    ax.set_xticks(positions)
    ax.set_xticklabels(MODULE_DISPLAY_ORDER, fontsize=9.5, rotation=15, ha="right")
    ax.set_ylabel("Incremental coverage", fontsize=11)
    ax.set_title(
        "NB-specific modules: incremental coverage across synergistic pairs",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )

    fig.text(
        0.5,
        0.01,
        "Incremental coverage = coverage(A∪B) − max[coverage(A), coverage(B)]",
        ha="center",
        fontsize=8.5,
        fontstyle="italic",
        color="gray",
    )

    legend_handles = [
        mpatches.Patch(
            color="#2ca02c", alpha=0.55, label="Differentiation (ADRN / Retinoid)"
        ),
        mpatches.Patch(color="#d62728", alpha=0.55, label="MYCN targets"),
        mpatches.Patch(color="#7f7f7f", alpha=0.55, label="Undifferentiated (MES)"),
    ]
    ax.legend(handles=legend_handles, fontsize=8.5, loc="upper right", framealpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0, 0.04, 1, 1])

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Figure 8 v15.png"
    out_pdf = OUT_DIR / "Figure 8 v15.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"ADRN vs MYCN targets (one-sided): U p = {p_adrn_mycn:.4f}")
    print(f"ADRN vs Retinoid (one-sided):     U p = {p_adrn_ret:.4f}")


if __name__ == "__main__":
    main()
