#!/usr/bin/env python3
"""
Combined forest plot — one row per miRNA pair showing MYCN-adjusted HR
for 'both high' co-expression, stratified on age. Uses penalized Cox
for pairs with complete separation.

Pairs: 124+34b, 137+450b, 19b+34b, 124+363.

History:
  - 2026-04-21 (evening): outputs colocated with scripts; paths now
    `./cox_forest_combined_v4.{png,svg}`. Centralized
    `multivariate_results/` sink removed.
  - 2026-04-21: moved from survival/figure6E_combined_forest_v4.py to
    survival/cox_forest/cox_forest_combined_v4.py per
    .state/NAMING_PLAN_v3.md. `os.chdir(Path(__file__).parent)` added.
  - 2026-04-17: created as v4 with median-split unification, Panel C
    pair swap (137+449b → 19b+34b), penalized Cox on complete
    separation, and DETECTION_CUTOFF. See .state/ledger.md 2026-04-17.
"""

import os
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from lifelines import CoxPHFitter
import warnings

warnings.filterwarnings("ignore")

os.chdir(Path(__file__).parent)

expr = pd.read_csv("data/miRNA_expression_data.csv")
surv = pd.read_csv("data/survival_data.csv")
df = surv.merge(expr, on="patient_id")
df["MYCN_amp"] = df["mycn_amplified_4.0"].astype(int)
df["age_over_18mo"] = df["over_18_months_age_of_diagnosis"].astype(int)

PAIRS = [
    ("miR-124-3p + miR-34b-5p", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("miR-137-3p + miR-450b-5p", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("miR-19b-3p + miR-34b-5p", ["hsa-miR-19b-3p", "hsa-miR-34b-5p"]),
    ("miR-124-3p + miR-363-3p", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
]

results = []

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


for label, mirna_pair in PAIRS:
    for m in mirna_pair:
        df[f"{m}_high"] = classify_high(df[m], m)
    df["both_high"] = (
        (df[f"{mirna_pair[0]}_high"] == 1) & (df[f"{mirna_pair[1]}_high"] == 1)
    ).astype(int)

    n_events_bh = df[df["both_high"] == 1]["event"].sum()
    penalizer = 0.1 if n_events_bh == 0 else 0.0

    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=0.0)
    cox_df = df[
        ["survival_time", "event", "both_high", "MYCN_amp", "age_over_18mo"]
    ].dropna()
    cph.fit(
        cox_df,
        duration_col="survival_time",
        event_col="event",
        strata=["age_over_18mo"],
        formula="both_high + MYCN_amp",
    )

    s = cph.summary
    results.append(
        {
            "label": label,
            "hr": s.loc["both_high", "exp(coef)"],
            "ci_lo": s.loc["both_high", "exp(coef) lower 95%"],
            "ci_hi": s.loc["both_high", "exp(coef) upper 95%"],
            "p": s.loc["both_high", "p"],
            "penalized": penalizer > 0,
        }
    )

    print(
        f"{label}: HR={results[-1]['hr']:.2f} "
        f"({results[-1]['ci_lo']:.2f}-{results[-1]['ci_hi']:.2f}), "
        f"p={results[-1]['p']:.4f}{' [penalized]' if penalizer > 0 else ''}"
    )

    for m in mirna_pair:
        df.drop(columns=[f"{m}_high"], inplace=True)
    df.drop(columns=["both_high"], inplace=True)

# ── Combined forest plot ─────────────────────────────────────────────────────
n_rows = len(results)
fig, ax = plt.subplots(figsize=(6, 3.0))
fig.subplots_adjust(right=0.45)

y_positions = np.arange(n_rows)[::-1]

for i, r in enumerate(results):
    y = y_positions[i]
    color = "#0072B2" if r["hr"] < 1 else "#D55E00"
    ax.plot(
        [r["ci_lo"], r["ci_hi"]],
        [y, y],
        color=color,
        linewidth=2.5,
        solid_capstyle="round",
    )
    marker = "D" if r["penalized"] else "s"
    ax.plot(
        r["hr"],
        y,
        marker,
        color=color,
        markersize=8,
        zorder=5,
        markeredgecolor="white",
        markeredgewidth=0.5,
    )

ax.axvline(x=1, color="#888888", linestyle="--", linewidth=0.8)

ax.set_xscale("log")
ax.set_xticks([0.01, 0.1, 1, 10])
ax.get_xaxis().set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:g}"))
ax.xaxis.set_minor_formatter(mticker.NullFormatter())
ax.set_xlim(0.02, 5)

ax.set_yticks(y_positions)
ax.set_yticklabels([r["label"] for r in results], fontsize=9)
ax.set_xlabel("Hazard Ratio (log scale)", fontsize=9)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", labelsize=8)

for i, r in enumerate(results):
    y = y_positions[i]
    pstr = f"{r['p']:.3f}" if r["p"] >= 0.001 else "< 0.001"
    hr_text = f"{r['hr']:.2f} ({r['ci_lo']:.2f}\u2013{r['ci_hi']:.2f})"
    if r["penalized"]:
        hr_text += "\u2020"
    ax.annotate(
        hr_text,
        xy=(1.05, y),
        xycoords=("axes fraction", "data"),
        fontsize=8.5,
        va="center",
        ha="left",
    )
    ax.annotate(
        pstr,
        xy=(1.72, y),
        xycoords=("axes fraction", "data"),
        fontsize=8.5,
        va="center",
        ha="left",
        fontweight="bold" if r["p"] < 0.05 else "normal",
    )

ax.annotate(
    "HR (95% CI)",
    xy=(1.05, 1.0),
    xycoords="axes fraction",
    fontsize=7.5,
    fontweight="bold",
    va="bottom",
    ha="left",
    color="#555555",
)
ax.annotate(
    "P",
    xy=(1.72, 1.0),
    xycoords="axes fraction",
    fontsize=7.5,
    fontweight="bold",
    va="bottom",
    ha="left",
    color="#555555",
)

if any(r["penalized"] for r in results):
    ax.text(
        0.0,
        -0.25,
        "\u2020 Penalized estimate (0 events in both-high group)",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#666666",
        va="top",
    )

plt.savefig("cox_forest_combined.png", dpi=300, bbox_inches="tight")
plt.savefig("cox_forest_combined.svg", bbox_inches="tight")
plt.savefig("cox_forest_combined.pdf", bbox_inches="tight")
plt.close()
print("\nSaved: cox_forest_combined_v14.{png,svg,pdf}")
