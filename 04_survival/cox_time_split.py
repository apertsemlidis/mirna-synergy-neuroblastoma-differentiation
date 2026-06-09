#!/usr/bin/env python3
"""
Generate the time-split Cox HR sensitivity table (Supplementary Table Sx).

For each of the four key miRNA pairs, fit the final age-stratified Cox
model (miRNA_both_high + MYCN_amp, strata=age) on:
  - full follow-up
  - early window  (event times <= SPLIT_DAYS, rest censored at SPLIT_DAYS)
  - late window   (patients surviving past SPLIT_DAYS, with entry time shifted)

Writes a tidy CSV and prints a markdown table to stdout.

Purpose: demonstrate that for pairs flagged by the Schoenfeld residual test
(Panels A and B in the current layout), the protective direction of the
miRNA effect is consistent in both halves of follow-up; Panel A shows a
delayed-effect pattern (stronger late), while Panel B's miRNA effect is
stable across time (its PH violation is in the MYCN nuisance covariate).

Pairs (the six dose-response pairs):
  124+363, 124+34b, 137+450b, 137+449b, 137+17, 19b+2110.

Output:
  ./cox_time_split_table.csv (colocated; already-tidy stats output)
"""

import os
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter
import warnings

warnings.filterwarnings("ignore")

os.chdir(Path(__file__).parent)

SPLIT_DAYS = 702  # median event time across all 21 events

# Final-model Figure 6 pairs (Panel label → two miRNAs)
PAIRS = [
    ("A", "124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
    ("B", "124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("C", "137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("D", "137+449b", ["hsa-miR-137-3p", "hsa-miR-449b-5p"]),
    ("E", "137+17", ["hsa-miR-137-3p", "hsa-miR-17-5p"]),
    ("F", "19b+2110", ["hsa-miR-19b-3p", "hsa-miR-2110"]),
]

# Kept consistent with the survival KM panels.
DETECTION_CUTOFF = {
    "hsa-miR-200a-3p": 1.0,
    "hsa-miR-211-5p": 1.0,
    "hsa-miR-429": 1.0,
    "hsa-miR-449a": 1.0,
    "hsa-miR-449b-5p": 1.0,
}


def classify_high(series, mirna):
    if mirna in DETECTION_CUTOFF:
        return (series > DETECTION_CUTOFF[mirna]).astype(int)
    return (series >= series.median()).astype(int)


def fit_cox(cox_df, penalizer):
    """Fit the final age-stratified Cox model; retry with penalizer=0.1 on
    convergence failure. Returns (summary_row_for_miRNA, penalized_flag) or
    (None, None)."""

    def _fit(pen):
        cph = CoxPHFitter(penalizer=pen, l1_ratio=0.0)
        cph.fit(
            cox_df,
            duration_col="T",
            event_col="E",
            strata=["age_over_18mo"],
            formula="miRNA_both_high + MYCN_amp",
        )
        return cph

    try:
        cph = _fit(penalizer)
        return cph.summary.loc["miRNA_both_high"], penalizer > 0
    except Exception:
        try:
            cph = _fit(0.1)
            return cph.summary.loc["miRNA_both_high"], True
        except Exception:
            return None, None


def run_window(cox_df, window):
    """Subset/transform cox_df for the given window, then fit."""
    if window == "full":
        sub = cox_df.copy()
    elif window == "early":
        sub = cox_df.copy()
        sub["E"] = ((sub["T"] <= SPLIT_DAYS) & (sub["E"] == 1)).astype(int)
        sub["T"] = sub["T"].clip(upper=SPLIT_DAYS)
    elif window == "late":
        sub = cox_df[cox_df["T"] > SPLIT_DAYS].copy()
        sub["T"] = sub["T"] - SPLIT_DAYS
    else:
        raise ValueError(f"unknown window {window!r}")

    n_total = len(sub)
    n_events = int(sub["E"].sum())
    n_bh = int((sub["miRNA_both_high"] == 1).sum())
    n_events_bh = int(sub[sub["miRNA_both_high"] == 1]["E"].sum())

    if n_bh == 0 or n_total == 0:
        return {
            "window": window,
            "n_total": n_total,
            "n_events": n_events,
            "n_bh": n_bh,
            "n_events_bh": n_events_bh,
            "HR": None,
            "CI_lo": None,
            "CI_hi": None,
            "p": None,
            "penalized": None,
            "note": "skipped (empty group)",
        }

    penalizer = 0.1 if n_events_bh == 0 else 0.0
    row, penalized = fit_cox(sub, penalizer)
    if row is None:
        return {
            "window": window,
            "n_total": n_total,
            "n_events": n_events,
            "n_bh": n_bh,
            "n_events_bh": n_events_bh,
            "HR": None,
            "CI_lo": None,
            "CI_hi": None,
            "p": None,
            "penalized": None,
            "note": "fit failed",
        }
    return {
        "window": window,
        "n_total": n_total,
        "n_events": n_events,
        "n_bh": n_bh,
        "n_events_bh": n_events_bh,
        "HR": float(row["exp(coef)"]),
        "CI_lo": float(row["exp(coef) lower 95%"]),
        "CI_hi": float(row["exp(coef) upper 95%"]),
        "p": float(row["p"]),
        "penalized": bool(penalized),
        "note": "",
    }


def main():
    expr = pd.read_csv("data/miRNA_expression_data.csv")
    surv = pd.read_csv("data/survival_data.csv")
    df = surv.merge(expr, on="patient_id")
    df["MYCN_amp"] = df["mycn_amplified_4.0"].astype(int)
    df["age_over_18mo"] = df["over_18_months_age_of_diagnosis"].astype(int)

    out_rows = []
    for panel, short, pair in PAIRS:
        df["a"] = classify_high(df[pair[0]], pair[0])
        df["b"] = classify_high(df[pair[1]], pair[1])
        df["miRNA_both_high"] = ((df["a"] == 1) & (df["b"] == 1)).astype(int)

        cox_df = (
            df[
                [
                    "survival_time",
                    "event",
                    "miRNA_both_high",
                    "MYCN_amp",
                    "age_over_18mo",
                ]
            ]
            .dropna()
            .copy()
        )
        cox_df.columns = ["T", "E", "miRNA_both_high", "MYCN_amp", "age_over_18mo"]

        for window in ("full", "early", "late"):
            result = run_window(cox_df, window)
            result["panel"] = panel
            result["pair"] = short
            result["pair_mirnas"] = " + ".join(pair)
            out_rows.append(result)

    out = pd.DataFrame(out_rows)[
        [
            "panel",
            "pair",
            "pair_mirnas",
            "window",
            "n_total",
            "n_events",
            "n_bh",
            "n_events_bh",
            "HR",
            "CI_lo",
            "CI_hi",
            "p",
            "penalized",
            "note",
        ]
    ]

    out_path = "cox_time_split_table.csv"
    out.to_csv(out_path, index=False)

    # Markdown rendering for supplement draft / stdout
    print("\n# Supplementary Table Sx — Time-split Cox HR sensitivity analysis")
    print(
        f"\nSplit at median event time = {SPLIT_DAYS} days. Final model: "
        "miRNA_both_high + MYCN_amp, stratified on age (>18mo at diagnosis)."
    )
    print(
        "\n| Panel | Pair | Window | n | events | n(both high) | events(bh) | "
        "HR (95% CI) | p | Penalized |"
    )
    print("|---|---|---|---:|---:|---:|---:|---|---:|:---:|")
    for _, r in out.iterrows():
        if r["HR"] is None:
            hr_ci = "—"
            pstr = "—"
        else:
            hr_ci = f"{r['HR']:.2f} ({r['CI_lo']:.2f}\u2013{r['CI_hi']:.2f})"
            pstr = f"{r['p']:.4f}" if r["p"] >= 0.001 else "< 0.001"
        pen = "*" if r["penalized"] else ""
        print(
            f"| {r['panel']} | {r['pair']} | {r['window']} | {r['n_total']} | "
            f"{r['n_events']} | {r['n_bh']} | {r['n_events_bh']} | {hr_ci} | {pstr} | {pen} |"
        )
    print(
        "\n`*` penalized Cox (ridge penalty 0.1) applied when the both-high "
        "group has zero events in the window, to avoid complete-separation "
        "divergence."
    )
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
