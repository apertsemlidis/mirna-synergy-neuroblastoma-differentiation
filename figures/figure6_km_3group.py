#!/usr/bin/env python3
"""
Figure 6 composite (v15) — Three-group Kaplan-Meier survival analysis
for four synergistic miRNA pairs.

Panels (2x2 grid, Option A — parallel content):
  A: miR-124-3p + miR-34b-5p
  B: miR-137-3p + miR-450b-5p
  C: miR-19b-3p + miR-34b-5p
  D: miR-124-3p + miR-363-3p

Re-render path: each panel is drawn natively in matplotlib (not stitched
from per-pair PNGs). Per-panel content mirrors
survival/km_3group/km_3group_all_v14.py — KM curves for 0/1/2 high
groups, Bonferroni-corrected pairwise log-rank p-values, number-at-risk
table — but sized for a 2x2 composite.

Output: figures/Figure 6 v15.{png,pdf}
"""

import warnings
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "04_survival" / "data"
OUT_DIR = ROOT / "figures"

PAIRS = [
    ("A", "124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("B", "137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("C", "19b+34b", ["hsa-miR-19b-3p", "hsa-miR-34b-5p"]),
    ("D", "124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
]

COLORS = {0: "#CC79A7", 1: "#009E73", 2: "#0072B2"}
LABELS = {0: "0 high", 1: "1 high", 2: "2 high"}
TIME_POINTS = [0, 1000, 2000, 3000, 4000, 5000]
N_COMPARISONS = 3

DETECTION_CUTOFF = {
    "hsa-miR-200a-3p": 1.0,
    "hsa-miR-211-5p": 1.0,
    "hsa-miR-429": 1.0,
    "hsa-miR-449a": 1.0,
    "hsa-miR-449b-5p": 1.0,
}


def short_name(mirna: str) -> str:
    return mirna.replace("hsa-", "")


def classify_high(series: pd.Series, mirna: str) -> pd.Series:
    if mirna in DETECTION_CUTOFF:
        return (series > DETECTION_CUTOFF[mirna]).astype(int)
    return (series >= series.median()).astype(int)


def load_data() -> pd.DataFrame:
    expr = pd.read_csv(DATA_DIR / "miRNA_expression_data.csv")
    surv = pd.read_csv(DATA_DIR / "survival_data.csv")
    return surv.merge(expr, on="patient_id")


def plot_km_3group_panel(
    ax, df: pd.DataFrame, panel_label: str, pair_short: str, mirna_pair: list[str]
) -> None:
    """Draw one 3-group KM panel into the given axes."""
    work = df.copy()
    for mirna in mirna_pair:
        work[f"{mirna}_high"] = classify_high(work[mirna], mirna)
    work["n_high"] = work[[f"{m}_high" for m in mirna_pair]].sum(axis=1)

    title = f"{short_name(mirna_pair[0])} + {short_name(mirna_pair[1])}"

    ns = {}
    for g in [0, 1, 2]:
        grp = work[work["n_high"] == g]
        ns[g] = len(grp)
        kmf = KaplanMeierFitter()
        kmf.fit(grp["survival_time"], grp["event"], label=f"{LABELS[g]} (n={ns[g]})")
        kmf.plot_survival_function(
            ax=ax,
            color=COLORS[g],
            linewidth=2.0,
            ci_show=False,
            show_censors=True,
            censor_styles={"ms": 5, "marker": "|"},
        )

    pval_strs = []
    for a, b in combinations([0, 1, 2], 2):
        ga = work[work["n_high"] == a]
        gb = work[work["n_high"] == b]
        lr = logrank_test(
            ga["survival_time"], gb["survival_time"], ga["event"], gb["event"]
        )
        pv_adj = min(lr.p_value * N_COMPARISONS, 1.0)
        pvs = f"{pv_adj:.3f}" if pv_adj >= 0.001 else "< 0.001"
        pval_strs.append(f"{a} vs {b}: p = {pvs}")
    pval_strs.append(f"(Bonferroni × {N_COMPARISONS})")

    ax.text(
        0.98,
        0.98,
        "\n".join(pval_strs),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="gray",
            alpha=0.9,
            linewidth=0.6,
        ),
    )

    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("Time (days)", fontsize=9)
    ax.set_ylabel("Survival Probability", fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(left=0)
    ax.legend(loc="lower left", fontsize=7, frameon=True, framealpha=0.9)
    ax.tick_params(axis="both", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Compact number-at-risk table below the x-axis.
    table_y_start = -0.28
    row_height = 0.06
    ax.text(
        0.5,
        table_y_start + row_height,
        "Number at risk",
        transform=ax.transAxes,
        fontsize=7.5,
        fontweight="bold",
        va="top",
        ha="center",
    )
    xmin, xmax_ax = ax.get_xlim()
    for gi, g in enumerate([0, 1, 2]):
        y_pos = table_y_start - gi * row_height
        ax.text(
            -0.02,
            y_pos,
            LABELS[g],
            transform=ax.transAxes,
            fontsize=6.5,
            fontweight="bold",
            va="top",
            ha="right",
            color=COLORS[g],
        )
        grp = work[work["n_high"] == g]
        for tp in TIME_POINTS:
            nar = (grp["survival_time"] >= tp).sum()
            x_frac = (tp - xmin) / (xmax_ax - xmin)
            ax.text(
                x_frac,
                y_pos,
                str(nar),
                transform=ax.transAxes,
                fontsize=6.5,
                va="top",
                ha="center",
                color=COLORS[g],
            )

    # Bold (A)/(B)/(C)/(D) panel letter, upper-left outside the title.
    ax.text(
        -0.16,
        1.06,
        f"({panel_label})",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def main() -> None:
    df = load_data()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, (letter, pair_short, mirna_pair) in zip(axes.flat, PAIRS):
        plot_km_3group_panel(ax, df, letter, pair_short, mirna_pair)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.55, wspace=0.30, bottom=0.10)

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Figure 6 v15.png"
    out_pdf = OUT_DIR / "Figure 6 v15.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
