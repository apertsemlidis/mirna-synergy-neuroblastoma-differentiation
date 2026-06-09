#!/usr/bin/env python3
"""
MYCN-stratified Kaplan-Meier survival for each of the six dose-response
pairs. Each pair produces a two-panel figure (MYCN non-amplified |
amplified). Also fits the final age-stratified Cox model per pair.

Pairs (v16 master set, all six dose-response plates):
  124+363, 124+34b, 137+450b, 137+449b, 137+17, 19b+2110.

v16 changes vs v14:
  - PAIRS swap: drop 19b+34b (no dose-response data); add 137+449b,
    137+17, 19b+2110.
  - Per-pair per-stratum log-rank + Cox HR statistics dumped to
    `km_mycn_stratified_stats.csv` (Path 4 refactor).
  - Per-pair output: `km_mycn_stratified_{pair}.{png,svg}` (unsuffixed).
"""

import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import warnings

warnings.filterwarnings("ignore")

os.chdir(Path(__file__).parent)

# ── Load and merge data ──────────────────────────────────────────────────────
expr = pd.read_csv("data/miRNA_expression_data.csv")
surv = pd.read_csv("data/survival_data.csv")
df = surv.merge(expr, on="patient_id")

df["MYCN_amp"] = df["mycn_amplified_4.0"].astype(int)
df["age_over_18mo"] = df["over_18_months_age_of_diagnosis"].astype(int)

# ── Define pairs ─────────────────────────────────────────────────────────────
PAIRS = [
    ("124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
    ("124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("137+449b", ["hsa-miR-137-3p", "hsa-miR-449b-5p"]),
    ("137+17", ["hsa-miR-137-3p", "hsa-miR-17-5p"]),
    ("19b+2110", ["hsa-miR-19b-3p", "hsa-miR-2110"]),
]

# Per-pair per-stratum log-rank + Cox HR stats — written at script end.
STATS_ROWS = []

TIME_POINTS = [0, 1000, 2000, 3000, 4000, 5000]
COLORS = {"not_both_high": "#0072B2", "both_high": "#D55E00"}


def short_name(mirna):
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


for pair_label, mirna_pair in PAIRS:
    # Classify patients
    for m in mirna_pair:
        df[f"{m}_high"] = classify_high(df[m], m)

    df["both_high"] = (
        (df[f"{mirna_pair[0]}_high"] == 1) & (df[f"{mirna_pair[1]}_high"] == 1)
    ).astype(int)

    title = f"{short_name(mirna_pair[0])} + {short_name(mirna_pair[1])}"

    # ── Cox PH summary ───────────────────────────────────────────────────────
    n_both_high = df["both_high"].sum()
    n_events_both_high = df[df["both_high"] == 1]["event"].sum()
    print(f"\n{'=' * 60}")
    print(
        f"{title}: n={len(df)}, both_high={n_both_high}, "
        f"events in both_high={n_events_both_high}"
    )

    # Penalized Cox for complete separation (Firth-like), standard otherwise
    penalizer = 0.1 if n_events_both_high == 0 else 0.0
    if penalizer > 0:
        print("  ** 0 events in both_high — penalized Cox (Firth-like)")
    try:
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
        hr = s.loc["both_high", "exp(coef)"]
        lo = s.loc["both_high", "exp(coef) lower 95%"]
        hi = s.loc["both_high", "exp(coef) upper 95%"]
        pv = s.loc["both_high", "p"]
        pstr = f"{pv:.3f}" if pv >= 0.001 else "< 0.001"
        pen_mark = "*" if penalizer > 0 else ""
        cox_hr_text = f"Cox HR {hr:.2f} ({lo:.2f}\u2013{hi:.2f}){pen_mark}"
        cox_p_text = f"p = {pstr}"
        pen_log = " [penalized]" if penalizer > 0 else ""
        print(
            f"  Cox HR (both_high): {hr:.2f} ({lo:.2f}-{hi:.2f}), p={pv:.4f}{pen_log}"
        )
        cph.print_summary()
        STATS_ROWS.append(
            {
                "pair": pair_label,
                "stratum": "cox_overall",
                "n_not_both": int((df["both_high"] == 0).sum()),
                "n_both": int(n_both_high),
                "events_not_both": int(df[df["both_high"] == 0]["event"].sum()),
                "events_both": int(n_events_both_high),
                "cox_HR": float(hr),
                "cox_CI_lo": float(lo),
                "cox_CI_hi": float(hi),
                "cox_p": float(pv),
                "penalized": bool(penalizer > 0),
                "logrank_p_raw": None,
                "logrank_p_bonferroni": None,
                "status": "OK",
            }
        )
    except Exception as e:
        print(f"  Cox model failed: {e}")
        cox_hr_text = "Cox: model failed"
        cox_p_text = ""
        STATS_ROWS.append(
            {
                "pair": pair_label,
                "stratum": "cox_overall",
                "n_not_both": int((df["both_high"] == 0).sum()),
                "n_both": int(n_both_high),
                "events_not_both": int(df[df["both_high"] == 0]["event"].sum()),
                "events_both": int(n_events_both_high),
                "cox_HR": None,
                "cox_CI_lo": None,
                "cox_CI_hi": None,
                "cox_p": None,
                "penalized": None,
                "logrank_p_raw": None,
                "logrank_p_bonferroni": None,
                "status": f"FAILED:{type(e).__name__}",
            }
        )

    # ── MYCN-stratified KM (compact layout) ─────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.02)

    for idx, (mycn_val, mycn_label) in enumerate(
        [(0, "MYCN Non-Amp"), (1, "MYCN Amp")]
    ):
        ax = axes[idx]
        subset = df[df["MYCN_amp"] == mycn_val]

        grp0 = subset[subset["both_high"] == 0]
        grp1 = subset[subset["both_high"] == 1]

        n0, n1 = len(grp0), len(grp1)

        kmf0 = KaplanMeierFitter()
        kmf1 = KaplanMeierFitter()

        # Per-stratum stats row — populated below depending on group sizes.
        stratum_row = {
            "pair": pair_label,
            "stratum": "non_amp" if mycn_val == 0 else "amp",
            "n_not_both": int(n0),
            "n_both": int(n1),
            "events_not_both": int(grp0["event"].sum()),
            "events_both": int(grp1["event"].sum()),
            "cox_HR": None,
            "cox_CI_lo": None,
            "cox_CI_hi": None,
            "cox_p": None,
            "penalized": None,
            "logrank_p_raw": None,
            "logrank_p_bonferroni": None,
            "status": "OK",
        }

        # Handle edge case: empty group
        if n1 == 0:
            stratum_row["status"] = "EMPTY:both_high"
            kmf0.fit(
                grp0["survival_time"], grp0["event"], label=f"Not both high (n={n0})"
            )
            kmf0.plot_survival_function(
                ax=ax, color=COLORS["not_both_high"], linewidth=2, ci_show=False
            )
            ax.text(
                0.95,
                0.95,
                "Both high: n=0",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="gray",
                    alpha=0.8,
                ),
            )
        elif n0 == 0:
            stratum_row["status"] = "EMPTY:not_both_high"
            kmf1.fit(grp1["survival_time"], grp1["event"], label=f"Both high (n={n1})")
            kmf1.plot_survival_function(
                ax=ax, color=COLORS["both_high"], linewidth=2, ci_show=False
            )
        else:
            kmf0.fit(
                grp0["survival_time"], grp0["event"], label=f"Not both high (n={n0})"
            )
            kmf1.fit(grp1["survival_time"], grp1["event"], label=f"Both high (n={n1})")

            kmf0.plot_survival_function(
                ax=ax, color=COLORS["not_both_high"], linewidth=2, ci_show=False
            )
            kmf1.plot_survival_function(
                ax=ax, color=COLORS["both_high"], linewidth=2, ci_show=False
            )

            lr = logrank_test(
                grp0["survival_time"],
                grp1["survival_time"],
                grp0["event"],
                grp1["event"],
            )
            # Bonferroni for 2 strata (non-amp + amp) per pair
            N_STRATA = 2
            pval_raw = lr.p_value
            pval_adj = min(pval_raw * N_STRATA, 1.0)
            pstr = f"p = {pval_adj:.4f}" if pval_adj >= 0.001 else "p < 0.001"
            pstr += f"\n(Bonferroni × {N_STRATA})"
            ax.text(
                0.95,
                0.95,
                pstr,
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="gray",
                    alpha=0.8,
                ),
            )
            print(
                f"  {mycn_label}: n={n0} vs {n1}, "
                f"log-rank raw p={pval_raw:.4f}, "
                f"Bonferroni-adj p={pval_adj:.4f}"
            )
            stratum_row["logrank_p_raw"] = float(pval_raw)
            stratum_row["logrank_p_bonferroni"] = float(pval_adj)

        STATS_ROWS.append(stratum_row)

        ax.set_title(mycn_label, fontsize=11, fontweight="bold")
        ax.set_xlabel("Time (days)", fontsize=10)
        if idx == 0:
            ax.set_ylabel("Survival Probability", fontsize=10)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlim(left=0)
        # Legend in lower-right — curves tend to flatten in the right-tail
        # region, while lower-left typically holds the crashing "not both
        # high" curve in the MYCN-amp panel. Per-panel (not shared) so
        # stratum-specific n values appear in each legend.
        ax.legend(loc="lower right", fontsize=8, frameon=True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.subplots_adjust(wspace=0.08)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outname = f"km_mycn_stratified_{pair_label}.png"
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.savefig(outname.replace(".png", ".svg"), bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outname}")

    # Clean up
    for m in mirna_pair:
        df.drop(columns=[f"{m}_high"], inplace=True)
    df.drop(columns=["both_high"], inplace=True)

print("\nAll MYCN-stratified panels generated.")

stats_df = pd.DataFrame(STATS_ROWS)
stats_csv = "km_mycn_stratified_stats.csv"
stats_df.to_csv(stats_csv, index=False)
print(
    f"Saved stats: {stats_csv}  ({len(stats_df)} rows = {len(PAIRS)} pairs × 3 records each)"
)
