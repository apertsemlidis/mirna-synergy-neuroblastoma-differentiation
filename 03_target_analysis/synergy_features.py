#!/usr/bin/env python3
"""
Systematic analysis of features distinguishing synergistic from
non-synergistic miRNA pairs (Suggestion #3).

Compares 34 dual-positive synergistic pairs against ~912 non-synergistic
combinations across multiple features: target overlap, individual potency,
target set size, and expression correlation in patient tumors.

v16 changes vs v14:
  - Per-feature Mann-Whitney U stats (medians, p-values, sample sizes)
    written to `synergy_features_stats.csv` (Path 4 refactor — figure
    composite consumes this rather than recomputing inline).
  - all_features.csv continues to be the per-pair raw feature table.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def load_screen_data(screen_dir):
    """Load NL and CBCA screen results."""
    nl = pd.read_csv(os.path.join(screen_dir, "nl_hsa_scores.csv"), index_col=0)
    cbca = pd.read_csv(os.path.join(screen_dir, "cbca_scores.csv"), index_col=0)
    return nl, cbca


def classify_pairs(nl_df, cbca_df):
    """Classify each pair as dual-positive synergistic or not."""
    # ATRA CBCA reference
    atra_cols = [str(i) for i in range(96)]
    valid_cols = [c for c in atra_cols if c in cbca_df.columns]
    atra_cbca = cbca_df.iloc[0][valid_cols].astype(float).values.mean()

    results = {}
    for combo in nl_df.index:
        # NL synergy: combo > HSA (two-sided p < 0.05)
        combo_nl = nl_df.loc[combo, ["combo_rep_1", "combo_rep_2", "combo_rep_3"]].astype(float).values
        hsa_nl = nl_df.loc[combo, ["HSA_rep_1", "HSA_rep_2", "HSA_rep_3"]].astype(float).values
        _, p_nl = stats.ttest_ind(combo_nl, hsa_nl, equal_var=False)
        ci = hsa_nl.mean() / combo_nl.mean() if combo_nl.mean() > 0 else np.nan
        nl_synergy = (ci < 1.0) and (p_nl < 0.05)

        # CBCA: combo < ATRA (one-sided p < 0.05)
        try:
            combo_cbca = cbca_df.loc[combo, ["combo_rep_1", "combo_rep_2", "combo_rep_3"]].astype(float).values
            _, p_cbca_2s = stats.ttest_1samp(combo_cbca, atra_cbca)
            p_cbca = p_cbca_2s / 2 if combo_cbca.mean() < atra_cbca else 1 - (p_cbca_2s / 2)
            cbca_pass = (combo_cbca.mean() < atra_cbca) and (p_cbca < 0.05)
        except:
            cbca_pass = False

        dual_pos = nl_synergy and cbca_pass

        results[combo] = {
            "dual_positive": dual_pos,
            "nl_synergy": nl_synergy,
            "combo_nl_mean": combo_nl.mean(),
            "hsa_nl_mean": hsa_nl.mean(),
            "nl_ci": ci,
            "nl_p": p_nl,
        }

    return pd.DataFrame(results).T


def compute_target_features(ts_path, combos):
    """Compute target-space features for each pair."""
    ts = pd.read_csv(ts_path, sep="\t")
    mir_col, gene_col = "miRNA", "GeneSymbol"

    # Build target sets for all miRNAs
    target_cache = {}
    for mirna in ts[mir_col].unique():
        targets = set(ts[ts[mir_col] == mirna][gene_col].dropna().unique())
        target_cache[mirna] = targets
        # Also store without hsa- prefix for matching
        short = mirna.replace("hsa-", "")
        target_cache[short] = targets

    features = {}
    for combo in combos:
        parts = combo.split(" + ")
        if len(parts) != 2:
            continue
        mirA, mirB = parts

        # Try multiple name formats
        tA = _get_targets(target_cache, mirA)
        tB = _get_targets(target_cache, mirB)

        if tA is None or tB is None:
            features[combo] = {
                "jaccard": np.nan, "overlap_count": np.nan,
                "size_A": len(tA) if tA else 0, "size_B": len(tB) if tB else 0,
                "size_union": np.nan, "size_ratio": np.nan,
            }
            continue

        union = tA | tB
        overlap = tA & tB
        jaccard = len(overlap) / len(union) if union else np.nan
        size_ratio = min(len(tA), len(tB)) / max(len(tA), len(tB)) if max(len(tA), len(tB)) > 0 else np.nan

        features[combo] = {
            "jaccard": jaccard,
            "overlap_count": len(overlap),
            "size_A": len(tA),
            "size_B": len(tB),
            "size_union": len(union),
            "size_ratio": size_ratio,
        }

    return pd.DataFrame(features).T


def _get_targets(cache, mirna):
    """Try multiple name formats to find targets."""
    for name in [mirna, mirna.replace("hsa-", ""), mirna.rsplit("-", 1)[0] if mirna.endswith(("3p", "5p")) else mirna]:
        if name in cache:
            return cache[name]
    return None


def compute_expression_correlation(expr_path, combos):
    """Compute expression correlation between miRNA pairs in patient tumors."""
    expr = pd.read_csv(expr_path)
    expr = expr.set_index("patient_id")

    corrs = {}
    for combo in combos:
        parts = combo.split(" + ")
        if len(parts) != 2:
            continue
        mirA, mirB = parts

        # Try exact match, then without -3p/-5p
        colA = _find_column(expr, mirA)
        colB = _find_column(expr, mirB)

        if colA is not None and colB is not None:
            r, p = stats.pearsonr(expr[colA], expr[colB])
            corrs[combo] = {"expr_corr": r, "expr_corr_p": p}
        else:
            corrs[combo] = {"expr_corr": np.nan, "expr_corr_p": np.nan}

    return pd.DataFrame(corrs).T


def _find_column(df, mirna):
    """Find miRNA column in expression data."""
    if mirna in df.columns:
        return mirna
    # Try without hsa- prefix variations
    for col in df.columns:
        if mirna.replace("hsa-", "") in col.replace("hsa-", ""):
            return col
    return None


def compute_potency_features(nl_df, combos):
    """Compute individual miRNA potency features."""
    features = {}
    for combo in combos:
        parts = combo.split(" + ")
        if len(parts) != 2:
            continue
        mirA, mirB = parts

        hsa_reps = nl_df.loc[combo, ["HSA_rep_1", "HSA_rep_2", "HSA_rep_3"]].astype(float).values
        best_single = hsa_reps.mean()

        features[combo] = {
            "best_single_nl": best_single,
        }

    return pd.DataFrame(features).T


def plot_feature_comparison(df, feature, ylabel, title, outpath):
    """Box plot comparing a feature between synergistic and non-synergistic pairs."""
    fig, ax = plt.subplots(figsize=(6, 5))

    syn = df[df["dual_positive"] == True][feature].dropna()
    nonsyn = df[df["dual_positive"] == False][feature].dropna()

    bp = ax.boxplot([nonsyn, syn], labels=["Non-synergistic", "Dual-positive\nsynergistic"],
                     patch_artist=True, widths=0.5)
    bp["boxes"][0].set_facecolor("#CCCCCC")
    bp["boxes"][1].set_facecolor("#FF8C00")
    bp["boxes"][0].set_alpha(0.6)
    bp["boxes"][1].set_alpha(0.6)

    # Overlay points
    for i, (data, color) in enumerate([(nonsyn, "gray"), (syn, "darkorange")]):
        jitter = np.random.normal(0, 0.04, len(data))
        ax.scatter(np.repeat(i + 1, len(data)) + jitter, data, color=color, alpha=0.3, s=15, zorder=3)

    # Stats
    if len(syn) > 3 and len(nonsyn) > 3:
        u, p = stats.mannwhitneyu(syn, nonsyn, alternative="two-sided")
        ax.set_title(f"{title}\n(p = {p:.3e}, Mann-Whitney U)", fontsize=11)
    else:
        ax.set_title(title)

    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

    return syn, nonsyn


def main():
    parser = argparse.ArgumentParser(description="Feature analysis of synergistic vs non-synergistic pairs")
    parser.add_argument("--screen-dir", required=True, help="Directory with screen CSVs")
    parser.add_argument("--targetscan-tsv", required=True, help="TargetScan predictions")
    parser.add_argument("--expr-csv", default=None, help="Patient miRNA expression CSV (optional)")
    parser.add_argument("--outdir", default="synergy_features", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- Load and classify ---
    print("Loading screen data...")
    nl_df, cbca_df = load_screen_data(args.screen_dir)
    classification = classify_pairs(nl_df, cbca_df)
    n_syn = classification["dual_positive"].sum()
    n_nonsyn = len(classification) - n_syn
    print(f"  {n_syn} dual-positive synergistic, {n_nonsyn} non-synergistic")

    # --- Target features ---
    print("Computing target-space features...")
    target_features = compute_target_features(args.targetscan_tsv, classification.index)
    df = classification.join(target_features)

    # --- Potency features ---
    print("Computing potency features...")
    potency = compute_potency_features(nl_df, classification.index)
    df = df.join(potency)

    # --- Expression correlation ---
    if args.expr_csv:
        print("Computing expression correlations...")
        expr_corr = compute_expression_correlation(args.expr_csv, classification.index)
        df = df.join(expr_corr)

    # --- Save combined data ---
    df.to_csv(os.path.join(args.outdir, "all_features.csv"))

    # --- Generate plots ---
    print("Generating comparison plots...")

    results = {}

    # 1. Jaccard similarity
    syn, nonsyn = plot_feature_comparison(
        df, "jaccard", "Jaccard similarity", "Target overlap (Jaccard)",
        os.path.join(args.outdir, "jaccard_comparison.png"))
    results["Jaccard similarity"] = (syn.median(), nonsyn.median(),
                                      stats.mannwhitneyu(syn.dropna(), nonsyn.dropna())[1])

    # 2. Target set size (union)
    syn, nonsyn = plot_feature_comparison(
        df, "size_union", "Union target count", "Combined target set size",
        os.path.join(args.outdir, "target_size_comparison.png"))
    if len(syn.dropna()) > 0 and len(nonsyn.dropna()) > 0:
        results["Union target size"] = (syn.median(), nonsyn.median(),
                                         stats.mannwhitneyu(syn.dropna(), nonsyn.dropna())[1])

    # 3. Size ratio (balance between partners)
    syn, nonsyn = plot_feature_comparison(
        df, "size_ratio", "Size ratio (min/max)", "Target set balance",
        os.path.join(args.outdir, "size_ratio_comparison.png"))
    if len(syn.dropna()) > 0 and len(nonsyn.dropna()) > 0:
        results["Size ratio"] = (syn.median(), nonsyn.median(),
                                  stats.mannwhitneyu(syn.dropna(), nonsyn.dropna())[1])

    # 4. Best single agent NL
    syn, nonsyn = plot_feature_comparison(
        df, "best_single_nl", "Best single-agent NL", "Individual miRNA potency",
        os.path.join(args.outdir, "potency_comparison.png"))
    results["Best single NL"] = (syn.median(), nonsyn.median(),
                                  stats.mannwhitneyu(syn.dropna(), nonsyn.dropna())[1])

    # 5. Expression correlation (if available)
    if "expr_corr" in df.columns:
        syn, nonsyn = plot_feature_comparison(
            df, "expr_corr", "Pearson r (tumor expression)", "miRNA co-expression in tumors",
            os.path.join(args.outdir, "expression_corr_comparison.png"))
        valid_syn = syn.dropna()
        valid_nonsyn = nonsyn.dropna()
        if len(valid_syn) > 0 and len(valid_nonsyn) > 0:
            results["Expression correlation"] = (valid_syn.median(), valid_nonsyn.median(),
                                                  stats.mannwhitneyu(valid_syn, valid_nonsyn)[1])

    # --- Multi-panel summary figure ---
    print("Generating summary panel...")
    features_to_plot = ["jaccard", "size_union", "best_single_nl"]
    labels = ["Jaccard\nsimilarity", "Union target\nset size", "Best single\nagent NL"]
    if "expr_corr" in df.columns:
        features_to_plot.append("expr_corr")
        labels.append("Expression\ncorrelation")

    n_panels = len(features_to_plot)
    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    for ax, feat, label in zip(axes, features_to_plot, labels):
        syn_vals = df[df["dual_positive"] == True][feat].dropna()
        nonsyn_vals = df[df["dual_positive"] == False][feat].dropna()

        bp = ax.boxplot([nonsyn_vals, syn_vals],
                         labels=["Non-syn\n(n={})".format(len(nonsyn_vals)),
                                  "Synergistic\n(n={})".format(len(syn_vals))],
                         patch_artist=True, widths=0.5)
        bp["boxes"][0].set_facecolor("#CCCCCC")
        bp["boxes"][1].set_facecolor("#FF8C00")
        bp["boxes"][0].set_alpha(0.6)
        bp["boxes"][1].set_alpha(0.6)

        if len(syn_vals) > 0 and len(nonsyn_vals) > 0:
            _, p = stats.mannwhitneyu(syn_vals, nonsyn_vals)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            ax.set_title(f"{label}\n{sig} (p={p:.2e})", fontsize=10)
        else:
            ax.set_title(label, fontsize=10)

    fig.suptitle("Features of synergistic vs non-synergistic miRNA pairs", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "features_summary_panel.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- Print summary ---
    print("\n" + "=" * 70)
    print("SUMMARY: Synergistic vs Non-Synergistic Pair Features")
    print("=" * 70)
    for feat, (syn_med, nonsyn_med, p) in results.items():
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {feat:30s}  syn={syn_med:.4f}  nonsyn={nonsyn_med:.4f}  p={p:.4e}  {sig}")

    # --- v16 stats CSV (Path 4 refactor) ---
    stats_rows = []
    feature_key_map = {
        "Jaccard similarity":      "jaccard",
        "Union target size":       "size_union",
        "Size ratio":              "size_ratio",
        "Best single NL":          "best_single_nl",
        "Expression correlation":  "expr_corr",
    }
    for feat_label, (syn_med, nonsyn_med, p) in results.items():
        feat_key = feature_key_map.get(feat_label, feat_label)
        n_syn = int(df[df["dual_positive"] == True][feat_key].dropna().shape[0])
        n_nonsyn = int(df[df["dual_positive"] == False][feat_key].dropna().shape[0])
        stats_rows.append({
            "feature": feat_label,
            "feature_key": feat_key,
            "n_synergistic": n_syn,
            "n_nonsynergistic": n_nonsyn,
            "median_synergistic": syn_med,
            "median_nonsynergistic": nonsyn_med,
            "mannwhitney_p": p,
        })
    stats_df = pd.DataFrame(stats_rows)
    stats_csv_path = os.path.join(args.outdir, "synergy_features_stats.csv")
    stats_df.to_csv(stats_csv_path, index=False)
    print(f"\nSaved stats: {stats_csv_path}  ({len(stats_df)} features)")

    print(f"\nResults saved to {args.outdir}/")


if __name__ == "__main__":
    main()
