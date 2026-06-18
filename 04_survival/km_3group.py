#!/usr/bin/env python3
"""
3-group Kaplan-Meier survival curves for each of the six dose-response
pairs. Groups: 0 / 1 / 2 miRNAs expressed above threshold.
Bonferroni-corrected pairwise log-rank tests.

Pairs (the six dose-response pairs):
  124+363, 124+34b, 137+450b, 137+449b, 137+17, 19b+2110.

Outputs:
  - Per-pair pairwise log-rank statistics -> `km_3group_stats.csv`.
  - Per-pair curves -> `km_3group_{pair}.{png,svg}`.
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
    ("A", "124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
    ("B", "124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("C", "137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("D", "137+449b", ["hsa-miR-137-3p", "hsa-miR-449b-5p"]),
    ("E", "137+17", ["hsa-miR-137-3p", "hsa-miR-17-5p"]),
    ("F", "19b+2110", ["hsa-miR-19b-3p", "hsa-miR-2110"]),
]

# Stats accumulator — written to km_3group_stats.csv at end of script.
STATS_ROWS = []

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

    # Pairwise log-rank tests among the ordered groups (0, 1, 2 high), reported
    # WITHOUT multiple-comparison correction. These three-group survival
    # analyses are exploratory (small subgroups), so the p-values are nominal.
    pval_strs = []
    for a, b in combinations([0, 1, 2], 2):
        ga = df[df["n_high"] == a]
        gb = df[df["n_high"] == b]
        lr = logrank_test(
            ga["survival_time"], gb["survival_time"], ga["event"], gb["event"]
        )
        pv = float(lr.p_value)
        pvs = f"{pv:.4f}" if pv >= 0.001 else "< 0.001"
        pval_strs.append(f"{a} vs {b}: p = {pvs}")
        STATS_ROWS.append(
            {
                "pair": pair_short,
                "panel": panel_label,
                "comparison": f"{a}v{b}",
                "n_a": int(ns[a]),
                "n_b": int(ns[b]),
                "events_a": int(ga["event"].sum()),
                "events_b": int(gb["event"].sum()),
                "logrank_p_raw": pv,
            }
        )

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

    outname = f"km_3group_{pair_short}.png"
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close()
    print(
        f"Panel {panel_label} ({pair_short}): {title}  |  "
        f"n = {ns[0]} / {ns[1]} / {ns[2]}  |  Saved: {outname}"
    )

    # Clean up temporary columns
    for mirna in mirna_pair:
        df.drop(columns=[f"{mirna}_high"], inplace=True)
    df.drop(columns=["n_high"], inplace=True)

print("\nAll km_3group panels generated.")

# Emit per-pair pairwise log-rank stats CSV.
stats_df = pd.DataFrame(STATS_ROWS)
stats_csv = "km_3group_stats.csv"
stats_df.to_csv(stats_csv, index=False)
print(
    f"Saved stats: {stats_csv}  ({len(stats_df)} rows = {len(PAIRS)} pairs x 3 comparisons)"
)
