#!/usr/bin/env python3
"""
Screen 3-group KM / both-high Cox effect sizes across all 31 synergistic
miRNA pairs in `mirna_pairs.csv`.

Per-pair output:
  - n per group (n_0 / n_1 / n_2), events in the "2 high" group
  - HR, 95% CI for `both_high` from an MYCN-adjusted, age-stratified
    Cox model (matches Model 4 of cox_multivariate/multivariate_survival_v14.py)
  - log-rank p (unadjusted, 2-high vs <2-high) — for reference only
  - `penalized` flag if Firth-like (penalizer=0.1) fallback was used

Framing: per-pair HR + CI reported for screening, not pair-level
significance testing. Formal hypothesis tests are reserved for the four
main-figure pairs (km_3group_all_v14.py).

Floor-effect miRNAs use the `DETECTION_CUTOFF` dict from
`km_3group_all_v14.py`; pairs containing miRNAs absent from the
expression dataset are reported with `status="MISSING:..."` and no HR.

Output: `km_3group_screen_v14.csv` (sorted by HR ascending).
"""

import os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")
os.chdir(Path(__file__).parent)

# ── Load data ────────────────────────────────────────────────────────────────
pairs_df = pd.read_csv("../03_target_analysis/synergistic_pairs.csv")
expr = pd.read_csv("data/miRNA_expression_data.csv")
surv = pd.read_csv("data/survival_data.csv")
df = surv.merge(expr, on="patient_id")

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


def short(m):
    return m.replace("hsa-miR-", "").replace("-5p", "").replace("-3p", "")


def fit_stratified_cox(data):
    """Fit MYCN-adjusted, age-stratified Cox on `both_high` + `MYCN_amp`
    (matches Model 4 of multivariate_survival_v14.py). Retry with penalizer=0.1
    on convergence failure.
    """
    for pen in (0.0, 0.1):
        try:
            cph = CoxPHFitter(penalizer=pen, l1_ratio=0.0)
            cph.fit(
                data,
                duration_col="T",
                event_col="E",
                strata=["age_over_18mo"],
            )
            return cph, pen > 0.0
        except Exception:
            continue
    return None, None


results = []
for _, row in pairs_df.iterrows():
    m1, m2 = row["mirA"], row["mirB"]
    pair_id = f"{short(m1)}+{short(m2)}"
    missing = [m for m in (m1, m2) if m not in expr.columns]
    if missing:
        results.append(
            {
                "pair": pair_id,
                "mirA": m1,
                "mirB": m2,
                "status": f"MISSING:{','.join(m.replace('hsa-', '') for m in missing)}",
                "n_0": None,
                "n_1": None,
                "n_2": None,
                "events_2": None,
                "HR": None,
                "CI_lo": None,
                "CI_hi": None,
                "p_logrank": None,
                "penalized": None,
            }
        )
        continue

    d = df.copy()
    d[f"{m1}_high"] = classify_high(d[m1], m1)
    d[f"{m2}_high"] = classify_high(d[m2], m2)
    d["n_high"] = d[f"{m1}_high"] + d[f"{m2}_high"]
    d["both_high"] = (d["n_high"] == 2).astype(int)

    n0 = int((d["n_high"] == 0).sum())
    n1 = int((d["n_high"] == 1).sum())
    n2 = int((d["n_high"] == 2).sum())
    ev2 = int(d.loc[d["both_high"] == 1, "event"].sum())

    # Log-rank (both_high vs not both_high)
    g1 = d[d["both_high"] == 1]
    g0 = d[d["both_high"] == 0]
    if n2 < 2 or n2 == 96:
        lr_p = None
    else:
        lr_p = float(
            logrank_test(
                g1["survival_time"], g0["survival_time"], g1["event"], g0["event"]
            ).p_value
        )

    # MYCN-adjusted, age-stratified Cox on both_high (Model 4 analog)
    cox_df = d[
        [
            "survival_time",
            "event",
            "both_high",
            "mycn_amplified_4.0",
            "over_18_months_age_of_diagnosis",
        ]
    ].copy()
    cox_df.columns = ["T", "E", "both_high", "MYCN_amp", "age_over_18mo"]
    cox_df = cox_df.dropna()

    hr = ci_lo = ci_hi = None
    penalized = None
    status = "OK"
    if n2 < 1 or ev2 < 1:
        status = "NO_EVENTS" if ev2 < 1 else "NO_BOTHHIGH"
    else:
        cph, penalized = fit_stratified_cox(cox_df)
        if cph is not None:
            hr = float(cph.hazard_ratios_["both_high"])
            ci = cph.confidence_intervals_
            ci_lo = float(np.exp(ci.loc["both_high"].iloc[0]))
            ci_hi = float(np.exp(ci.loc["both_high"].iloc[1]))
        else:
            status = "COX_FAILED"

    results.append(
        {
            "pair": pair_id,
            "mirA": m1,
            "mirB": m2,
            "status": status,
            "n_0": n0,
            "n_1": n1,
            "n_2": n2,
            "events_2": ev2,
            "HR": hr,
            "CI_lo": ci_lo,
            "CI_hi": ci_hi,
            "p_logrank": lr_p,
            "penalized": penalized,
        }
    )

out = pd.DataFrame(results)
# Sort: OK rows by HR ascending, then other statuses
out["_ok"] = (out["status"] == "OK").astype(int)
out["_hr_sort"] = out["HR"].fillna(np.inf)
out = (
    out.sort_values(["_ok", "_hr_sort"], ascending=[False, True])
    .drop(columns=["_ok", "_hr_sort"])
    .reset_index(drop=True)
)

out.to_csv("km_3group_screen_v14.csv", index=False)

# ── Console summary ──────────────────────────────────────────────────────────
print(
    f"\n{'pair':<16} {'n0/n1/n2':<12} {'ev2':>4}  {'HR':>6} {'95% CI':>16}  {'p_LR':>7}  status"
)
print("-" * 80)
for _, r in out.iterrows():
    n_str = f"{r['n_0']}/{r['n_1']}/{r['n_2']}" if pd.notna(r["n_0"]) else "—"
    ev_str = f"{int(r['events_2'])}" if pd.notna(r["events_2"]) else "—"
    if pd.notna(r["HR"]):
        hr_str = f"{r['HR']:.2f}"
        ci_str = f"({r['CI_lo']:.2f}–{r['CI_hi']:.2f})"
    else:
        hr_str, ci_str = "—", "—"
    p_str = f"{r['p_logrank']:.3f}" if pd.notna(r["p_logrank"]) else "—"
    pen = "*" if r["penalized"] else ""
    print(
        f"{r['pair']:<16} {n_str:<12} {ev_str:>4}  {hr_str:>6}{pen:<1} {ci_str:>15}  {p_str:>7}  {r['status']}"
    )
print("\n* = penalized Cox (Firth-like, penalizer=0.1) used")
print(f"\nSaved: km_3group_screen_v14.csv  ({len(out)} pairs)")
