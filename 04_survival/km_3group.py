#!/usr/bin/env python3
"""
3-group Kaplan-Meier survival curves for each of the four key miRNA pairs.
Groups: 0 / 1 / 2 miRNAs expressed above threshold. Bonferroni-corrected
pairwise log-rank tests.

Pairs: 124+34b, 137+450b, 19b+34b, 124+363.

History:
  - 2026-04-21 (evening): outputs colocated with scripts; output paths now
    `./km_3group_{pair}_v4.{png,svg}` (relative to this script's directory).
    The centralized `multivariate_results/` sink was removed.
  - 2026-04-21: moved from survival/figure6_km_panels_v4.py to
    survival/km_3group/km_3group_all_v4.py per .state/NAMING_PLAN_v3.md.
    PAIRS now carries a short-pair-ID field; panel letters retained only
    in console output. `os.chdir(Path(__file__).parent)` added so paths
    resolve from this script's own directory.
  - 2026-04-17: created as v4 with median-split unification, Panel C
    pair swap (137+449b → 19b+34b), Bonferroni ×3 on pairwise log-rank,
    and DETECTION_CUTOFF for floor-effect miRNAs. See
    .state/ledger.md 2026-04-17 entry.
"""

import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from itertools import combinations
import warnings

warnings.filterwarnings("ignore")

os.chdir(Path(__file__).parent)

# ── Load and merge data ──────────────────────────────────────────────────────
expr = pd.read_csv("data/miRNA_expression_data.csv")
surv = pd.read_csv("data/survival_data.csv")
df = surv.merge(expr, on="patient_id")

# ── Define pairs (panel_label, short_id, miRNA list) ────────────────────────
PAIRS = [
    ("A", "124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("B", "137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("C", "19b+34b", ["hsa-miR-19b-3p", "hsa-miR-34b-5p"]),
    ("D", "124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
]

# Colorblind-safe palette for 3 groups (0, 1, 2 high)
COLORS = {0: "#CC79A7", 1: "#009E73", 2: "#0072B2"}
LABELS = {0: "0 high miRNAs", 1: "1 high miRNA", 2: "2 high miRNAs"}
TIME_POINTS = [0, 1000, 2000, 3000, 4000, 5000]


def make_short_name(mirna):
    """hsa-miR-124-3p -> miR-124-3p"""
    return mirna.replace("hsa-", "")


# These miRNAs have a floor value of 1.0 in a majority of patients, so
# median-splitting is not informative. Use a detectable-expression cutoff (>1)
# instead. Thresholds were verified against miRNA_expression_data.csv.
DETECTION_CUTOFF = {
    "hsa-miR-200a-3p": 1.0,  # 71/96 tied at floor
    "hsa-miR-211-5p": 1.0,  # 87/96 tied at floor
    "hsa-miR-429": 1.0,  # 81/96 tied at floor
    "hsa-miR-449a": 1.0,  # 61/96 tied at floor
    "hsa-miR-449b-5p": 1.0,  # 88/96 tied at floor
}


def classify_high(series, mirna):
    if mirna in DETECTION_CUTOFF:
        return (series > DETECTION_CUTOFF[mirna]).astype(int)
    return (series >= series.median()).astype(int)


for panel_label, pair_short, mirna_pair in PAIRS:
    # Classify patients
    for mirna in mirna_pair:
        df[f"{mirna}_high"] = classify_high(df[mirna], mirna)

    df["n_high"] = df[[f"{m}_high" for m in mirna_pair]].sum(axis=1)

    short_names = [make_short_name(m) for m in mirna_pair]
    title = f"{short_names[0]} + {short_names[1]}"

    # ── Plot ─────────────────────────────────────────────────────────────────
    # Wider figure to reserve room for legend + p-value box outside the axes.
    fig, ax = plt.subplots(figsize=(9, 6))

    ns = {}
    for g in [0, 1, 2]:
        grp = df[df["n_high"] == g]
        ns[g] = len(grp)
        kmf = KaplanMeierFitter()
        kmf.fit(grp["survival_time"], grp["event"], label=f"{LABELS[g]} (n={ns[g]})")
        kmf.plot_survival_function(
            ax=ax,
            color=COLORS[g],
            linewidth=2.5,
            ci_show=False,
            show_censors=True,
            censor_styles={"ms": 8, "marker": "|"},
        )

    # Pairwise log-rank tests (Bonferroni-corrected for 3 comparisons)
    N_COMPARISONS = 3
    pval_strs = []
    for a, b in combinations([0, 1, 2], 2):
        ga = df[df["n_high"] == a]
        gb = df[df["n_high"] == b]
        lr = logrank_test(
            ga["survival_time"], gb["survival_time"], ga["event"], gb["event"]
        )
        pv_adj = min(lr.p_value * N_COMPARISONS, 1.0)
        pvs = f"{pv_adj:.4f}" if pv_adj >= 0.001 else "< 0.001"
        pval_strs.append(f"{a} vs {b}: p = {pvs}")
    pval_strs.append(f"(Bonferroni × {N_COMPARISONS})")

    # p-value box — outside axes, upper right
    ax.text(
        1.02,
        0.98,
        "\n".join(pval_strs),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.9
        ),
    )

    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.set_ylabel("Survival Probability", fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(left=0)
    # Legend — outside axes, below the p-value box
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 0.62),
        fontsize=10,
        frameon=True,
        borderaxespad=0,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Number-at-risk table ─────────────────────────────────────────────────
    table_y_start = -0.18
    row_height = 0.05

    ax.text(
        0.5,
        table_y_start + 0.05,
        "Number at risk",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="center",
    )

    for gi, g in enumerate([0, 1, 2]):
        y_pos = table_y_start - gi * row_height
        ax.text(
            -0.02,
            y_pos,
            LABELS[g],
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            ha="right",
            color=COLORS[g],
        )
        grp = df[df["n_high"] == g]
        xmin, xmax_ax = ax.get_xlim()
        for tp in TIME_POINTS:
            nar = (grp["survival_time"] >= tp).sum()
            x_frac = (tp - xmin) / (xmax_ax - xmin)
            ax.text(
                x_frac,
                y_pos,
                str(nar),
                transform=ax.transAxes,
                fontsize=8,
                va="top",
                ha="center",
                color=COLORS[g],
            )

    plt.subplots_adjust(bottom=0.22, right=0.78)

    outname = f"km_3group_{pair_short}_v14.png"
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.savefig(outname.replace(".png", ".svg"), bbox_inches="tight")
    plt.close()
    print(
        f"Panel {panel_label} ({pair_short}): {title}  |  "
        f"n = {ns[0]} / {ns[1]} / {ns[2]}  |  Saved: {outname}"
    )

    # Clean up temporary columns
    for mirna in mirna_pair:
        df.drop(columns=[f"{mirna}_high"], inplace=True)
    df.drop(columns=["n_high"], inplace=True)

print("\nAll km_3group panels (v4) generated.")
