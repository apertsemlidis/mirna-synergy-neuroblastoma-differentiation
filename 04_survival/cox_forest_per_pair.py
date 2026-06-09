#!/usr/bin/env python3
"""
Per-pair Cox forest plots for the four key miRNA pairs. Model:
miRNA_both_high + MYCN_amp, stratified on age > 18 months. Uses
penalized Cox (Firth-like, penalizer=0.1) for pairs with complete
separation. Log-scale HR axis; hand-rolled HR-space rendering with
HR/p sidecar columns.

Pairs: 124+34b, 137+450b, 19b+34b, 124+363.

History:
  - 2026-04-21 (evening): outputs colocated with scripts; paths now
    `./cox_forest_{pair}.png`. Centralized `multivariate_results/`
    sink removed.
  - 2026-04-21: moved from survival/figure6E_forest_all_v4.py to
    survival/cox_forest/cox_forest_all_v4.py per
    an internal directory reorganization. `os.chdir(Path(__file__).parent)` added.
  - 2026-04-17: created as v4 with median-split unification, Panel C
    pair swap (137+449b → 19b+34b), penalized Cox on complete
    separation, and DETECTION_CUTOFF. See the project ledger 2026-04-17.
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
    ("124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("19b+34b", ["hsa-miR-19b-3p", "hsa-miR-34b-5p"]),
    ("124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
]


def short_name(m):
    return m.replace("hsa-", "")


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


for pair_label, mirna_pair in PAIRS:
    for m in mirna_pair:
        df[f"{m}_high"] = classify_high(df[m], m)
    df["both_high"] = (
        (df[f"{mirna_pair[0]}_high"] == 1) & (df[f"{mirna_pair[1]}_high"] == 1)
    ).astype(int)

    title = f"{short_name(mirna_pair[0])} + {short_name(mirna_pair[1])}"
    n_events_bh = df[df["both_high"] == 1]["event"].sum()

    penalizer = 0.1 if n_events_bh == 0 else 0.0
    if penalizer > 0:
        print(f"\n{title}: 0 events in both_high — penalized Cox (Firth-like)")

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
    covariates = s.index.tolist()
    hrs = s["exp(coef)"].values
    ci_lo = s["exp(coef) lower 95%"].values
    ci_hi = s["exp(coef) upper 95%"].values
    pvals = s["p"].values

    label_map = {
        "both_high": f"{short_name(mirna_pair[0])} +\n{short_name(mirna_pair[1])} both high",
        "MYCN_amp": "MYCN amplification",
    }
    labels = [label_map.get(c, c) for c in covariates]
    colors = ["#0072B2" if hr < 1 else "#D55E00" for hr in hrs]

    print(
        f"{title}: both_high HR={hrs[0]:.2f} ({ci_lo[0]:.2f}-{ci_hi[0]:.2f}), "
        f"p={pvals[0]:.4f}{' [penalized]' if penalizer > 0 else ''}"
    )

    fig, ax = plt.subplots(figsize=(5, 2.2))
    fig.subplots_adjust(right=0.42)

    y_pos = np.arange(len(covariates))[::-1]

    for i in range(len(covariates)):
        y = y_pos[i]
        ax.plot(
            [ci_lo[i], ci_hi[i]],
            [y, y],
            color=colors[i],
            linewidth=2.5,
            solid_capstyle="round",
        )
        ax.plot(
            hrs[i],
            y,
            "s",
            color=colors[i],
            markersize=8,
            zorder=5,
            markeredgecolor="white",
            markeredgewidth=0.5,
        )

    ax.axvline(x=1, color="#888888", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.01, 0.1, 1, 10, 100])
    ax.get_xaxis().set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Hazard Ratio (log scale)", fontsize=9)
    ax.set_xlim(0.008, 200)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=8)

    for i in range(len(covariates)):
        y = y_pos[i]
        pv = pvals[i]
        pstr = f"{pv:.3f}" if pv >= 0.001 else "< 0.001"
        hr_text = f"{hrs[i]:.2f} ({ci_lo[i]:.2f}\u2013{ci_hi[i]:.2f})"
        if penalizer > 0 and covariates[i] == "both_high":
            hr_text += "*"
        ax.annotate(
            hr_text,
            xy=(1.08, y),
            xycoords=("axes fraction", "data"),
            fontsize=8.5,
            va="center",
            ha="left",
        )
        ax.annotate(
            pstr,
            xy=(1.88, y),
            xycoords=("axes fraction", "data"),
            fontsize=8.5,
            va="center",
            ha="left",
            fontweight="bold" if pv < 0.05 else "normal",
        )

    ax.annotate(
        "HR (95% CI)",
        xy=(1.08, 1.0),
        xycoords="axes fraction",
        fontsize=7.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#555555",
    )
    ax.annotate(
        "P",
        xy=(1.88, 1.0),
        xycoords="axes fraction",
        fontsize=7.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#555555",
    )

    outname = f"cox_forest_{pair_label}.png"
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outname}")

    for m in mirna_pair:
        df.drop(columns=[f"{m}_high"], inplace=True)
    df.drop(columns=["both_high"], inplace=True)

print("\nAll 6E forest plots (v4) generated.")
