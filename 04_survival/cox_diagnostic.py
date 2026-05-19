#!/usr/bin/env python3
"""
Single-pair Cox PH diagnostic for miR-137-3p + miR-450b-5p. Same
methodology as cox_forest/cox_forest_all_v4.py (the pair's manuscript
Cox forest). Console-only output; no figures written.

History:
  - 2026-04-21: moved from survival/cox_mir137_mir450b_v4.py to
    survival/cox_diagnostic/cox_diagnostic_137+450b_v4.py per
    .state/NAMING_PLAN_v3.md. `os.chdir(Path(__file__).parent)` added;
    input paths updated to `../miRNA_expression_data.csv` and
    `../survival_data.csv`. The `_diagnostic_` infix now identifies
    this as the single-pair deep-dive script category (distinct from
    the `_preview_` tag on lifelines-native renders in cox_multivariate/).
  - 2026-04-17: created as v4 with penalized Cox on complete
    separation, captured PH assumption test, unified `>=` median-split
    convention, and Bonferroni x2 on MYCN-stratum log-rank tests.
    See .state/ledger.md 2026-04-17.
"""

import os
from pathlib import Path

import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
import warnings

warnings.filterwarnings("ignore")

os.chdir(Path(__file__).parent)

# ── Load and merge data ─────────────────────────────────────────────────────
expr = pd.read_csv("data/miRNA_expression_data.csv")
surv = pd.read_csv("data/survival_data.csv")
df = surv.merge(
    expr[["patient_id", "hsa-miR-137-3p", "hsa-miR-450b-5p"]], on="patient_id"
)

# ── Classify patients (>= median; neither miRNA is in DETECTION_CUTOFF) ────
med_137 = df["hsa-miR-137-3p"].median()
med_450b = df["hsa-miR-450b-5p"].median()
print(f"Median miR-137-3p: {med_137:.4f}")
print(f"Median miR-450b-5p: {med_450b:.4f}")

df["miR137_high"] = (df["hsa-miR-137-3p"] >= med_137).astype(int)
df["miR450b_high"] = (df["hsa-miR-450b-5p"] >= med_450b).astype(int)
df["miRNA_both_high"] = ((df["miR137_high"] == 1) & (df["miR450b_high"] == 1)).astype(
    int
)
df["n_high"] = df["miR137_high"] + df["miR450b_high"]

df["MYCN_amp"] = df["mycn_amplified_4.0"].astype(int)
df["age_over_18mo"] = df["over_18_months_age_of_diagnosis"].astype(int)

n_total = len(df)
n_events = int(df["event"].sum())
print(f"\nTotal patients: {n_total}, Events: {n_events}")
print(
    f"Both high: {df['miRNA_both_high'].sum()}, "
    f"Not both high: {(~df['miRNA_both_high'].astype(bool)).sum()}"
)
print(
    f"MYCN amplified: {df['MYCN_amp'].sum()}, "
    f"Non-amplified: {(~df['MYCN_amp'].astype(bool)).sum()}"
)
print(
    f"Age > 18mo: {df['age_over_18mo'].sum()}, "
    f"Age <= 18mo: {(~df['age_over_18mo'].astype(bool)).sum()}"
)

print("\n=== Group sizes ===")
print(
    f"  Both high & event: {((df['miRNA_both_high'] == 1) & (df['event'] == 1)).sum()}"
)
print(
    f"  Both high & no event: {((df['miRNA_both_high'] == 1) & (df['event'] == 0)).sum()}"
)
print(
    f"  Not both high & event: {((df['miRNA_both_high'] == 0) & (df['event'] == 1)).sum()}"
)
print(
    f"  Not both high & no event: {((df['miRNA_both_high'] == 0) & (df['event'] == 0)).sum()}"
)

# ── Cox PH model, stratified on age_over_18mo ──────────────────────────────
print("\n" + "=" * 70)
print("Cox PH Model (stratified on age_over_18mo)")
print("Covariates: miRNA_both_high + MYCN_amp")
print("=" * 70)

cox_df = df[
    ["survival_time", "event", "miRNA_both_high", "MYCN_amp", "age_over_18mo"]
].dropna()
n_events_bh = int(cox_df[cox_df["miRNA_both_high"] == 1]["event"].sum())
penalizer = 0.1 if n_events_bh == 0 else 0.0
if penalizer > 0:
    print(f"** 0 events in both_high group — penalized Cox (penalizer={penalizer})")

cph = CoxPHFitter(penalizer=penalizer, l1_ratio=0.0)
cph.fit(
    cox_df,
    duration_col="survival_time",
    event_col="event",
    strata=["age_over_18mo"],
    formula="miRNA_both_high + MYCN_amp",
)

cph.print_summary()

summary = cph.summary
for cov in summary.index:
    hr = summary.loc[cov, "exp(coef)"]
    ci_lo = summary.loc[cov, "exp(coef) lower 95%"]
    ci_hi = summary.loc[cov, "exp(coef) upper 95%"]
    p = summary.loc[cov, "p"]
    pen_mark = "*" if penalizer > 0 else ""
    print(
        f"\n  {cov}: HR={hr:.3f} (95% CI: {ci_lo:.3f}-{ci_hi:.3f}){pen_mark}, p={p:.4f}"
    )

n_events_cox = int(cox_df["event"].sum())
n_params = 2
epv = n_events_cox / n_params
print(f"\nEvents-per-variable (EPV): {n_events_cox}/{n_params} = {epv:.1f}")

# ── PH assumption: Schoenfeld residuals test, captured and evaluated ───────
print("\n" + "=" * 70)
print("Proportional Hazards Assumption Test (Schoenfeld residuals)")
print("=" * 70)

ph_result = proportional_hazard_test(cph, cox_df, time_transform="rank")
ph_summary = ph_result.summary
print(ph_summary.to_string())

PH_THRESHOLD = 0.05
violations = ph_summary[ph_summary["p"] < PH_THRESHOLD]
if len(violations) == 0:
    print(
        f"\n[PH OK] All covariate p-values > {PH_THRESHOLD}; "
        "proportional hazards assumption SATISFIED."
    )
    ph_status = "satisfied"
else:
    print(
        f"\n[PH VIOLATION] {len(violations)} covariate(s) with Schoenfeld "
        f"p < {PH_THRESHOLD}: {list(violations.index)}"
    )
    ph_status = "violated"

# ── MYCN-stratified KM with log-rank (Bonferroni × 2 for 2 strata) ─────────
print("\n" + "=" * 70)
print("MYCN-Stratified Kaplan-Meier + Log-Rank Tests (Bonferroni \u00d7 2)")
print("=" * 70)

for mycn_status, label in [(0, "MYCN non-amplified"), (1, "MYCN amplified")]:
    subset = df[df["MYCN_amp"] == mycn_status]
    n_sub = len(subset)
    n_ev = int(subset["event"].sum())

    grp_high = subset[subset["miRNA_both_high"] == 1]
    grp_low = subset[subset["miRNA_both_high"] == 0]

    print(f"\n--- {label} (n={n_sub}, events={n_ev}) ---")
    print(f"  Both high: n={len(grp_high)}, events={int(grp_high['event'].sum())}")
    print(f"  Not both high: n={len(grp_low)}, events={int(grp_low['event'].sum())}")

    if len(grp_high) > 0 and len(grp_low) > 0:
        result = logrank_test(
            grp_high["survival_time"],
            grp_low["survival_time"],
            event_observed_A=grp_high["event"],
            event_observed_B=grp_low["event"],
        )
        p_raw = result.p_value
        p_adj = min(p_raw * 2, 1.0)
        print(
            f"  Log-rank test: chi2={result.test_statistic:.3f}, "
            f"raw p={p_raw:.4f}, Bonferroni-adj p={p_adj:.4f}"
        )

        kmf = KaplanMeierFitter()
        kmf.fit(grp_high["survival_time"], grp_high["event"])
        med_high = kmf.median_survival_time_
        kmf.fit(grp_low["survival_time"], grp_low["event"])
        med_low = kmf.median_survival_time_
        print(
            f"  Median survival — both high: {med_high:.0f}, "
            f"not both high: {med_low:.0f}"
        )
    else:
        print("  Skipped (empty group)")

# ── Overall (unstratified) KM ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("Overall KM: both_high vs not_both_high (unstratified)")
print("=" * 70)

grp_high = df[df["miRNA_both_high"] == 1]
grp_low = df[df["miRNA_both_high"] == 0]

result = logrank_test(
    grp_high["survival_time"],
    grp_low["survival_time"],
    event_observed_A=grp_high["event"],
    event_observed_B=grp_low["event"],
)
print(f"Log-rank test: chi2={result.test_statistic:.3f}, p={result.p_value:.4f}")

kmf = KaplanMeierFitter()
kmf.fit(grp_high["survival_time"], grp_high["event"])
med_high = kmf.median_survival_time_
kmf.fit(grp_low["survival_time"], grp_low["event"])
med_low = kmf.median_survival_time_
print(f"Median survival — both high: {med_high:.0f}, not both high: {med_low:.0f}")

# ── 3-group KM (0, 1, 2 high) ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("3-Group KM: 0 vs 1 vs 2 miRNAs high")
print("=" * 70)

for n in [0, 1, 2]:
    grp = df[df["n_high"] == n]
    print(f"  n_high={n}: n={len(grp)}, events={int(grp['event'].sum())}")
    kmf = KaplanMeierFitter()
    kmf.fit(grp["survival_time"], grp["event"])
    print(f"    Median survival: {kmf.median_survival_time_:.0f}")

print(f"\nDone. PH assumption: {ph_status}.")
