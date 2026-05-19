#!/usr/bin/env python3
"""
Multivariate survival analysis for miRNA co-expression in neuroblastoma.

Tests whether the prognostic association between synergistic miRNA pair
co-expression and overall survival holds after adjusting for known clinical
risk factors (MYCN amplification, age at diagnosis). Fits four Cox models
per pair (univariate, +MYCN, +MYCN+age as covariate, +MYCN strata=age)
plus PH check on the final (strata=age) model.

Pairs: 124+363, 124+34b, 137+450b, 19b+34b.

Data: GSE155945 (Misiak et al., 2021) — 96 neuroblastoma tumors.

History:
  - 2026-04-21 (evening): outputs colocated with scripts; CLI default
    --output-dir changed to "." (this script's directory after chdir).
    Centralized `multivariate_results/` sink removed.
  - 2026-04-21: moved from survival/multivariate_survival_v4.py to
    survival/cox_multivariate/multivariate_survival_v4.py per
    .state/NAMING_PLAN_v3.md. `os.chdir(Path(__file__).parent)` added;
    CLI default --data-dir changed from "." to ".." (data CSVs still
    live at survival/ root). Output filenames changed: long miRNA
    names replaced with short pair IDs; non-unique renders
    (km_univariate_, cox_forest_, km_mycn_stratified_) tagged with
    `_preview_` to distinguish from the manuscript-grade km_3group /
    cox_forest / km_mycn_stratified outputs produced by their
    dedicated scripts.
  - 2026-04-17: created as v4 with median-split unification, Panel C
    pair swap (137+449b → 19b+34b), penalized Cox on complete
    separation, PH assumption capture via proportional_hazard_test,
    Model 4 (age-stratified) added, and _safe_fit helper for
    convergence robustness. See .state/ledger.md 2026-04-17.
"""

import argparse
import os
import warnings
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test

warnings.filterwarnings("ignore")

os.chdir(Path(__file__).parent)

# Short-pair IDs used in output filenames (consistent across all v4 scripts).
SHORT_PAIR_ID = {
    ("hsa-miR-124-3p", "hsa-miR-363-3p"): "124+363",
    ("hsa-miR-124-3p", "hsa-miR-34b-5p"): "124+34b",
    ("hsa-miR-137-3p", "hsa-miR-450b-5p"): "137+450b",
    ("hsa-miR-19b-3p", "hsa-miR-34b-5p"): "19b+34b",
}


def short_pair(mirna_pair):
    """Short ID for a miRNA pair, e.g. '124+363'. Preserves order."""
    key = tuple(mirna_pair)
    return SHORT_PAIR_ID.get(
        key, "_".join(m.replace("hsa-miR-", "") for m in mirna_pair)
    )


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


# --- Key miRNA pairs matching Figure 6 panels A–D ---
KEY_PAIRS = [
    ["hsa-miR-124-3p", "hsa-miR-363-3p"],  # Main exemplar pair (Panel D)
    ["hsa-miR-124-3p", "hsa-miR-34b-5p"],  # Panel A
    ["hsa-miR-137-3p", "hsa-miR-450b-5p"],  # Panel B
    ["hsa-miR-19b-3p", "hsa-miR-34b-5p"],  # Panel C (replaces 137+449b)
]


def load_data(data_dir):
    """Load and merge expression + survival data."""
    expr = pd.read_csv(os.path.join(data_dir, "miRNA_expression_data.csv"))
    surv = pd.read_csv(os.path.join(data_dir, "survival_data.csv"))
    data = pd.merge(expr, surv, on="patient_id")
    return data, expr


def compute_coexpression(data, expr, mirna_pair):
    """Classify patients by number of high-expression miRNAs (0, 1, or 2)."""
    for mirna in mirna_pair:
        data[f"{mirna}_high"] = classify_high(data[mirna], mirna)
    data["high_count"] = data[[f"{m}_high" for m in mirna_pair]].sum(axis=1)
    data["both_high"] = (data["high_count"] == len(mirna_pair)).astype(int)
    return data


def run_univariate_km(data, mirna_pair, output_dir):
    """KM curves stratified by co-expression."""
    pair_label = " + ".join([m.replace("hsa-", "") for m in mirna_pair])
    groups = sorted(data["high_count"].unique())

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = sns.color_palette("husl", len(groups))

    for i, group in enumerate(groups):
        gd = data[data["high_count"] == group]
        kmf = KaplanMeierFitter()
        kmf.fit(
            gd["survival_time"],
            gd["event"],
            label=f"{int(group)} high miRNAs (n={len(gd)})",
        )
        kmf.plot_survival_function(ax=ax, ci_show=False, linewidth=2, color=colors[i])

    # Pairwise log-rank with Bonferroni correction
    pairs_to_test = [(g1, g2) for g1 in groups for g2 in groups if g1 < g2]
    n_tests = len(pairs_to_test)
    pvals = {}
    for g1, g2 in pairs_to_test:
        d1 = data[data["high_count"] == g1]
        d2 = data[data["high_count"] == g2]
        result = logrank_test(
            d1["survival_time"], d2["survival_time"], d1["event"], d2["event"]
        )
        pvals[f"{int(g1)} vs {int(g2)}"] = min(result.p_value * n_tests, 1.0)

    pstring = "\n".join([f"{k}: p = {v:.2e}" for k, v in pvals.items()])
    pstring += f"\n(Bonferroni \u00d7 {n_tests})"
    ax.text(
        0.65,
        0.05,
        pstring,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        fontfamily="monospace",
    )

    ax.set_title(f"Univariate KM: {pair_label}", fontsize=13)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Survival probability")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, f"km_3group_preview_{short_pair(mirna_pair)}_v14.png"),
        dpi=300,
    )
    plt.close()

    print(f"\n--- Univariate KM: {pair_label} ---")
    for k, v in pvals.items():
        print(f"  {k}: p = {v:.4f}")


def run_multivariate_cox(data, mirna_pair, output_dir):
    """Cox proportional hazards model adjusting for clinical covariates.

    Applies penalized Cox (penalizer=0.1) when the both_high group has zero
    events, to avoid complete-separation divergence.
    """
    pair_label = " + ".join([m.replace("hsa-", "") for m in mirna_pair])

    cox_df = data[
        [
            "survival_time",
            "event",
            "both_high",
            "mycn_amplified_4.0",
            "over_18_months_age_of_diagnosis",
        ]
    ].copy()
    cox_df.columns = ["T", "E", "miRNA_both_high", "MYCN_amp", "age_over_18mo"]
    cox_df = cox_df.dropna()

    n_events_bh = cox_df[cox_df["miRNA_both_high"] == 1]["E"].sum()
    penalizer = 0.1 if n_events_bh == 0 else 0.0
    if penalizer > 0:
        print(
            f"  ** {pair_label}: 0 events in both_high — "
            "penalized Cox (Firth-like) applied to all models"
        )

    def _safe_fit(label, fit_fn):
        """Fit a Cox model; if it diverges, retry with penalizer=0.1."""
        try:
            return fit_fn(penalizer), False
        except Exception as e:
            print(
                f"  ** {label}: convergence issue ({type(e).__name__}); "
                "retrying with penalizer=0.1"
            )
            try:
                return fit_fn(0.1), True
            except Exception as e2:
                print(f"  ** {label}: penalized fit also failed ({e2})")
                return None, None

    # --- Model 1: miRNA co-expression only ---
    print(f"\n--- Cox Model 1 (univariate): {pair_label} ---")
    cph1, _ = _safe_fit(
        "Model 1",
        lambda pen: CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(
            cox_df[["T", "E", "miRNA_both_high"]], duration_col="T", event_col="E"
        ),
    )
    if cph1 is not None:
        cph1.print_summary()

    # --- Model 2: miRNA + MYCN ---
    print(f"\n--- Cox Model 2 (+ MYCN): {pair_label} ---")
    cph2, _ = _safe_fit(
        "Model 2",
        lambda pen: CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(
            cox_df[["T", "E", "miRNA_both_high", "MYCN_amp"]],
            duration_col="T",
            event_col="E",
        ),
    )
    if cph2 is not None:
        cph2.print_summary()

    # --- Model 3: miRNA + MYCN + age (age as covariate) ---
    print(f"\n--- Cox Model 3 (+ MYCN + age): {pair_label} ---")
    cph3, _ = _safe_fit(
        "Model 3",
        lambda pen: CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(
            cox_df[["T", "E", "miRNA_both_high", "MYCN_amp", "age_over_18mo"]],
            duration_col="T",
            event_col="E",
        ),
    )
    if cph3 is not None:
        cph3.print_summary()

    # --- Model 4: miRNA + MYCN, stratified by age (authoritative final model,
    #     matches figure6E_forest_all_v4.py and the manuscript forest plot) ---
    print(f"\n--- Cox Model 4 (+ MYCN, strata=age): {pair_label} ---")
    cph4, _ = _safe_fit(
        "Model 4",
        lambda pen: CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(
            cox_df,
            duration_col="T",
            event_col="E",
            strata=["age_over_18mo"],
            formula="miRNA_both_high + MYCN_amp",
        ),
    )
    if cph4 is not None:
        cph4.print_summary()

    # --- PH assumption test on the final model (Schoenfeld residuals) ---
    print(f"\n--- PH assumption test (Model 4, Schoenfeld residuals): {pair_label} ---")
    ph_status = "not tested"
    if cph4 is not None:
        try:
            ph_result = proportional_hazard_test(cph4, cox_df, time_transform="rank")
            ph_summary = ph_result.summary
            print(ph_summary.to_string())
            PH_THRESHOLD = 0.05
            violations = ph_summary[ph_summary["p"] < PH_THRESHOLD]
            if len(violations) == 0:
                ph_status = "satisfied"
                print(f"[PH OK] all Schoenfeld p > {PH_THRESHOLD}")
            else:
                ph_status = f"violated: {list(violations.index)}"
                print(
                    f"[PH VIOLATION] covariates with p < {PH_THRESHOLD}: {list(violations.index)}"
                )
        except Exception as e:
            ph_status = f"test failed: {e}"
            print(f"[PH TEST FAILED] {e}")
    else:
        print("[PH SKIPPED] Model 4 did not fit")

    if cph4 is not None:
        fig, ax = plt.subplots(figsize=(8, 4))
        cph4.plot(ax=ax)
        ax.set_title(
            f"Cox PH Model: {pair_label}\nadjusted for MYCN, stratified by age",
            fontsize=12,
        )
        ax.axvline(x=0, color="black", linestyle="--", linewidth=0.5)
        fig.tight_layout()
        fig.savefig(
            os.path.join(
                output_dir, f"cox_forest_preview_{short_pair(mirna_pair)}_v14.png"
            ),
            dpi=300,
        )
        plt.close()

    summary = []
    model_list = [
        ("Univariate", cph1),
        ("+ MYCN", cph2),
        ("+ MYCN + Age", cph3),
        ("+ MYCN, strata=age", cph4),
    ]
    for name, model in model_list:
        if model is None:
            continue
        s = model.summary
        if "miRNA_both_high" in s.index:
            row = s.loc["miRNA_both_high"]
            summary.append(
                {
                    "Model": name,
                    "HR": row["exp(coef)"],
                    "HR_lower": row["exp(coef) lower 95%"],
                    "HR_upper": row["exp(coef) upper 95%"],
                    "p": row["p"],
                    "penalized": penalizer > 0,
                    "ph": ph_status if name == "+ MYCN, strata=age" else "",
                }
            )
    summary_df = pd.DataFrame(summary)
    print(f"\n--- Summary: {pair_label} ---")
    print(summary_df.to_string(index=False))
    summary_df.to_csv(
        os.path.join(output_dir, f"cox_summary_{short_pair(mirna_pair)}_v14.csv"),
        index=False,
    )

    return summary_df


def run_stratified_km(data, mirna_pair, output_dir):
    """KM curves stratified by MYCN status."""
    pair_label = " + ".join([m.replace("hsa-", "") for m in mirna_pair])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax, (mycn_val, mycn_label) in zip(
        axes, [(0, "MYCN non-amplified"), (1, "MYCN amplified")]
    ):
        subset = data[data["mycn_amplified_4.0"] == mycn_val]
        n_both = (subset["both_high"] == 1).sum()
        n_not = (subset["both_high"] == 0).sum()

        if n_both < 3 or n_not < 3:
            ax.set_title(f"{mycn_label}\n(too few patients to stratify)")
            ax.text(
                0.5,
                0.5,
                f"n={len(subset)}\nboth_high={n_both}",
                transform=ax.transAxes,
                ha="center",
            )
            continue

        for val, label, color in [
            (1, f"Both high (n={n_both})", "tab:red"),
            (0, f"Not both high (n={n_not})", "tab:blue"),
        ]:
            gd = subset[subset["both_high"] == val]
            kmf = KaplanMeierFitter()
            kmf.fit(gd["survival_time"], gd["event"], label=label)
            kmf.plot_survival_function(ax=ax, ci_show=True, linewidth=2, color=color)

        d1 = subset[subset["both_high"] == 1]
        d0 = subset[subset["both_high"] == 0]
        lr = logrank_test(
            d1["survival_time"], d0["survival_time"], d1["event"], d0["event"]
        )
        # Bonferroni for 2 strata (non-amp + amp)
        pval_adj = min(lr.p_value * 2, 1.0)
        ax.text(
            0.65,
            0.05,
            f"log-rank p = {pval_adj:.3f}\n(Bonferroni \u00d7 2)",
            transform=ax.transAxes,
            fontsize=10,
        )

        ax.set_title(f"{mycn_label} (n={len(subset)})")
        ax.set_xlabel("Time (days)")
        ax.legend(frameon=False, fontsize=9)

    axes[0].set_ylabel("Survival probability")
    fig.suptitle(f"MYCN-stratified KM: {pair_label}", fontsize=14)
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            output_dir,
            f"km_mycn_stratified_preview_{short_pair(mirna_pair)}_v14.png",
        ),
        dpi=300,
    )
    plt.close()

    print(f"\n--- MYCN-stratified KM saved: {pair_label} ---")


def main():
    parser = argparse.ArgumentParser(
        description="Multivariate survival analysis for miRNA pairs (v4)"
    )
    parser.add_argument(
        "--data-dir", default="..", help="Directory with expression/survival CSVs"
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Output directory",
    )
    parser.add_argument(
        "--pairs",
        default="key",
        choices=["key", "all"],
        help="Analyze key pairs only or all pairs from the study",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    data, expr = load_data(args.data_dir)
    print(f"Loaded {len(data)} patients, {data.event.sum()} events")
    print(f"MYCN amplified: {data['mycn_amplified_4.0'].sum()}")
    print(f"Age >18mo: {data['over_18_months_age_of_diagnosis'].sum()}")

    pairs = KEY_PAIRS

    all_summaries = []
    for pair in pairs:
        pair_label = " + ".join(pair)
        missing = [m for m in pair if m not in expr.columns]
        if missing:
            print(f"\nSkipping {pair_label}: missing {missing}")
            continue

        data = compute_coexpression(data, expr, pair)
        run_univariate_km(data, pair, args.output_dir)
        summary = run_multivariate_cox(data, pair, args.output_dir)
        run_stratified_km(data, pair, args.output_dir)
        all_summaries.append((pair_label, summary))

    print("\n" + "=" * 70)
    print("OVERALL SUMMARY (v4)")
    print("=" * 70)
    for pair_label, summary in all_summaries:
        print(f"\n{pair_label}:")
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
