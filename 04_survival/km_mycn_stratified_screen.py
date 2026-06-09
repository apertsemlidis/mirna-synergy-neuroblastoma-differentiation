#!/usr/bin/env python3
"""
Screen MYCN-stratified Cox effect sizes across all 31 synergistic miRNA pairs.

Per-pair output (direct analog of the two-panel main figure):
  - n per stratum × group (n_nonamp_both_high / n_amp_both_high)
  - HR + 95% CI for `both_high` within MYCN-nonamp stratum (age-stratified Cox)
  - HR + 95% CI for `both_high` within MYCN-amp stratum (age-stratified Cox)
  - log-rank p per stratum (unadjusted reference)
  - interaction_p: LR test for `both_high × MYCN_amp` (pooled model) —
    formally tests whether the miRNA effect differs between MYCN strata
  - penalized flags per stratum

Caveat: the MYCN-amp stratum has n=17 total and ~11 events; with a further
split by `both_high`, many pairs will hit NO_EVENTS in the amp stratum and
the interaction test will be underpowered. Report these as descriptive
effect sizes, not pair-level significance tests.

Framing matches km_3group_screen.py: screen, not hypothesis test.
Output: `km_mycn_stratified_screen.csv` (sorted by HR_amp ascending).
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
df["MYCN_amp"] = df["mycn_amplified_4.0"].astype(int)
df["age_over_18mo"] = df["over_18_months_age_of_diagnosis"].astype(int)

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


def fit_stratum_cox(data):
    """Age-stratified Cox on `both_high` within one MYCN stratum.
    Retries with penalizer=0.1 on convergence failure; returns (HR, CI_lo,
    CI_hi, penalized_flag, status)."""
    n_bh = int(data["both_high"].sum())
    events_bh = int(data.loc[data["both_high"] == 1, "E"].sum())
    if n_bh < 1 or events_bh < 1 or n_bh == len(data):
        return (
            None,
            None,
            None,
            None,
            (
                "NO_EVENTS"
                if events_bh < 1
                else "NO_BOTHHIGH"
                if n_bh < 1
                else "ALL_BOTHHIGH"
            ),
        )
    # Check age variation within stratum (can't stratify on a constant)
    strata_kw = {}
    if data["age_over_18mo"].nunique() > 1:
        strata_kw["strata"] = ["age_over_18mo"]
    for pen in (0.0, 0.1):
        try:
            cph = CoxPHFitter(penalizer=pen, l1_ratio=0.0)
            cph.fit(data, duration_col="T", event_col="E", **strata_kw)
            hr = float(cph.hazard_ratios_["both_high"])
            ci = cph.confidence_intervals_
            ci_lo = float(np.exp(ci.loc["both_high"].iloc[0]))
            ci_hi = float(np.exp(ci.loc["both_high"].iloc[1]))
            return hr, ci_lo, ci_hi, pen > 0.0, "OK"
        except Exception:
            continue
    return None, None, None, None, "COX_FAILED"


def interaction_p(data):
    """LR test for both_high × MYCN_amp interaction in a pooled, age-stratified
    Cox model. Returns None if the interaction can't be fit (e.g., no events
    in a cell)."""
    # Fit reduced (no interaction) and full (with interaction) models
    base = data[["T", "E", "both_high", "MYCN_amp", "age_over_18mo"]].copy()
    base["bh_x_mycn"] = base["both_high"] * base["MYCN_amp"]
    for pen in (0.0, 0.1):
        try:
            m_red = CoxPHFitter(penalizer=pen).fit(
                base.drop(columns=["bh_x_mycn"]),
                duration_col="T",
                event_col="E",
                strata=["age_over_18mo"],
            )
            m_full = CoxPHFitter(penalizer=pen).fit(
                base, duration_col="T", event_col="E", strata=["age_over_18mo"]
            )
            # Likelihood-ratio test
            lr_stat = 2 * (m_full.log_likelihood_ - m_red.log_likelihood_)
            from scipy.stats import chi2

            return float(chi2.sf(lr_stat, df=1))
        except Exception:
            continue
    return None


# ── Screen ───────────────────────────────────────────────────────────────────
rows = []
for _, pair_row in pairs_df.iterrows():
    m1, m2 = pair_row["mirA"], pair_row["mirB"]
    pair_id = f"{short(m1)}+{short(m2)}"
    missing = [m for m in (m1, m2) if m not in expr.columns]
    if missing:
        rows.append(
            {
                "pair": pair_id,
                "mirA": m1,
                "mirB": m2,
                "status": f"MISSING:{','.join(m.replace('hsa-', '') for m in missing)}",
            }
        )
        continue

    d = df.copy()
    d[f"{m1}_high"] = classify_high(d[m1], m1)
    d[f"{m2}_high"] = classify_high(d[m2], m2)
    d["both_high"] = ((d[f"{m1}_high"] == 1) & (d[f"{m2}_high"] == 1)).astype(int)

    result = {"pair": pair_id, "mirA": m1, "mirB": m2, "status": "OK"}

    for stratum_val, stratum_name in [(0, "nonamp"), (1, "amp")]:
        sub = d[d["MYCN_amp"] == stratum_val][
            ["survival_time", "event", "both_high", "age_over_18mo"]
        ].copy()
        sub.columns = ["T", "E", "both_high", "age_over_18mo"]
        sub = sub.dropna()

        n_bh = int(sub["both_high"].sum())
        ev_bh = int(sub.loc[sub["both_high"] == 1, "E"].sum())
        result[f"n_{stratum_name}_not_bh"] = len(sub) - n_bh
        result[f"n_{stratum_name}_bh"] = n_bh
        result[f"events_{stratum_name}_bh"] = ev_bh

        hr, ci_lo, ci_hi, pen, status = fit_stratum_cox(sub)
        result[f"HR_{stratum_name}"] = hr
        result[f"CI_lo_{stratum_name}"] = ci_lo
        result[f"CI_hi_{stratum_name}"] = ci_hi
        result[f"penalized_{stratum_name}"] = pen
        result[f"status_{stratum_name}"] = status

        # Log-rank within stratum (unadjusted reference)
        if n_bh > 0 and n_bh < len(sub) and ev_bh > 0:
            g1 = sub[sub["both_high"] == 1]
            g0 = sub[sub["both_high"] == 0]
            result[f"p_logrank_{stratum_name}"] = float(
                logrank_test(g1["T"], g0["T"], g1["E"], g0["E"]).p_value
            )
        else:
            result[f"p_logrank_{stratum_name}"] = None

    # Interaction test (pooled model)
    pool = d[
        ["survival_time", "event", "both_high", "MYCN_amp", "age_over_18mo"]
    ].copy()
    pool.columns = ["T", "E", "both_high", "MYCN_amp", "age_over_18mo"]
    pool = pool.dropna()
    result["interaction_p"] = interaction_p(pool)

    rows.append(result)

out = pd.DataFrame(rows)
# Sort: OK rows by HR_amp ascending, others to the end
out["_ok"] = out["status"].fillna("").eq("OK").astype(int)
out["_hr_sort"] = out.get("HR_amp", pd.Series(dtype=float)).fillna(np.inf)
out = (
    out.sort_values(["_ok", "_hr_sort"], ascending=[False, True])
    .drop(columns=["_ok", "_hr_sort"])
    .reset_index(drop=True)
)
out.to_csv("km_mycn_stratified_screen.csv", index=False)

# ── Console summary ──────────────────────────────────────────────────────────
print(
    f"\n{'pair':<16} {'n_na(bh)':<10} {'HR_na (CI)':<20} "
    f"{'n_a(bh)':<9} {'HR_a (CI)':<20} {'int_p':>7}  status"
)
print("-" * 100)


def fmt_hr(hr, lo, hi):
    if pd.isna(hr):
        return "—"
    return f"{hr:.2f} ({lo:.2f}–{hi:.2f})"


for _, r in out.iterrows():
    if r["status"] != "OK":
        print(f"{r['pair']:<16}  {r['status']}")
        continue
    n_na = f"{r['n_nonamp_not_bh']}({r['n_nonamp_bh']})"
    n_a = f"{r['n_amp_not_bh']}({r['n_amp_bh']})"
    hr_na = fmt_hr(r["HR_nonamp"], r["CI_lo_nonamp"], r["CI_hi_nonamp"])
    hr_a = fmt_hr(r["HR_amp"], r["CI_lo_amp"], r["CI_hi_amp"])
    ip = f"{r['interaction_p']:.3f}" if pd.notna(r["interaction_p"]) else "—"
    print(f"{r['pair']:<16} {n_na:<10} {hr_na:<20} {n_a:<9} {hr_a:<20} {ip:>7}  OK")

print(
    "\nColumns: n_na(bh) = nonamp total (both_high); "
    "n_a(bh) = amp total (both_high). CI are 95%."
)
print(f"\nSaved: km_mycn_stratified_screen.csv  ({len(out)} pairs)")
