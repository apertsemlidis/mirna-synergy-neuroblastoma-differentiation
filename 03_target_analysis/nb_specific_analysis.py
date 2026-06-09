#!/usr/bin/env python3
"""
NB-specific target-space analysis for synergistic miRNA pairs.

Extends the generic pathway analysis (v8) with neuroblastoma-specific gene
modules: MYCN transcriptional targets, adrenergic/mesenchymal cell identity
signatures, and retinoid response genes.

Reuses the same coverage/complementarity metrics as the generic analysis
but loads gene sets from local files rather than Enrichr.

v16 changes vs v14:
  - Per-category Mann-Whitney U stats (medians, p-values, sample sizes)
    written to `nb_specific_stats.csv` (Path 4 refactor — figure
    composite consumes this rather than recomputing inline).
  - all_pairs_nb_metrics.csv continues to be the per-pair per-module
    incremental-coverage table.

Gene signatures:
  - MYCN targets: WEI_MYCN_TARGETS_WITH_E_BOX (MSigDB, Wei et al. 2008)
  - ADRN identity: van Groningen et al. 2017 Nat Genet, Table S2
  - MES identity: van Groningen et al. 2017 Nat Genet, Table S2
  - Retinoid response: GO_RESPONSE_TO_RETINOIC_ACID (MSigDB/GO)

Dependencies: pandas, numpy, matplotlib, scipy
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# --- Target database specs ---

@dataclass(frozen=True)
class TargetSource:
    name: str
    mir_col: str
    gene_col: str
    score_col: str
    score_threshold: float
    score_direction: str  # "lt" or "gt"


TARGETSCAN_SPEC = TargetSource(
    name="TargetScan",
    mir_col="miRNA",
    gene_col="GeneSymbol",
    score_col=None,
    score_threshold=0,
    score_direction="lt",
)


# --- Module definitions ---

NB_MODULES = [
    {
        "module": "MYCN transcriptional targets",
        "file": "mycn_targets_wei.txt",
        "display": "MYCN targets (Wei et al.)",
        "category": "mycn",
        "source": "MSigDB: WEI_MYCN_TARGETS_WITH_E_BOX",
    },
    {
        "module": "Adrenergic identity (ADRN)",
        "file": "adrn_van_groningen.txt",
        "display": "Adrenergic (ADRN)",
        "category": "nb_differentiation",
        "source": "van Groningen et al. 2017 Table S2",
    },
    {
        "module": "Mesenchymal identity (MES)",
        "file": "mes_van_groningen.txt",
        "display": "Mesenchymal (MES)",
        "category": "nb_undifferentiated",
        "source": "van Groningen et al. 2017 Table S2",
    },
    {
        "module": "Retinoid response",
        "file": "retinoid_response_gobp.txt",
        "display": "Retinoid response (GO)",
        "category": "nb_differentiation",
        "source": "GO_RESPONSE_TO_RETINOIC_ACID",
    },
]

CATEGORY_ORDER = ["nb_differentiation", "mycn", "nb_undifferentiated"]
CATEGORY_LABELS = {
    "nb_differentiation": "NB differentiation",
    "mycn": "MYCN program",
    "nb_undifferentiated": "NB undifferentiated",
}


# --- Data loading ---

def load_targets(tsv_path: str, spec: TargetSource) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    for col in [spec.mir_col, spec.gene_col]:
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found. Available: {list(df.columns)[:10]}")
    if spec.score_col and spec.score_col in df.columns:
        df[spec.score_col] = pd.to_numeric(df[spec.score_col], errors="coerce")
    return df


def get_targets(df: pd.DataFrame, spec: TargetSource, mirna: str) -> Set[str]:
    # Try exact match first, then without -3p/-5p suffix
    query = mirna.replace("hsa-", "")
    sub = df[df[spec.mir_col].str.contains(query, case=False, na=False)]
    if len(sub) == 0:
        # Strip -3p/-5p suffix for TargetScan family-level matching
        base = query.rsplit("-", 1)[0] if query.endswith(("3p", "5p")) else query
        sub = df[df[spec.mir_col].str.contains(base, case=False, na=False)]
    if spec.score_col and spec.score_col in df.columns:
        if spec.score_direction == "lt":
            sub = sub[sub[spec.score_col] < spec.score_threshold]
        else:
            sub = sub[sub[spec.score_col] > spec.score_threshold]
    return set(sub[spec.gene_col].dropna().unique())


def get_universe(df: pd.DataFrame, spec: TargetSource) -> Set[str]:
    return set(df[spec.gene_col].dropna().unique())


def load_gene_list(filepath: str) -> Set[str]:
    with open(filepath) as f:
        return {line.strip() for line in f if line.strip()}


def load_modules(sig_dir: str) -> Tuple[List[dict], Dict[str, Set[str]]]:
    module_meta = []
    module_gene_sets = {}

    for md in NB_MODULES:
        filepath = os.path.join(sig_dir, md["file"])
        if not os.path.exists(filepath):
            print(f"WARNING: {filepath} not found, skipping {md['module']}", file=sys.stderr)
            continue
        genes = load_gene_list(filepath)
        module_gene_sets[md["module"]] = genes
        module_meta.append({
            "module": md["module"],
            "display": md["display"],
            "category": md["category"],
            "source": md["source"],
            "file": md["file"],
            "n_genes": len(genes),
        })

    return module_meta, module_gene_sets


# --- Metrics ---

def safe_div(num: float, den: float) -> float:
    return float(num) / float(den) if den else float("nan")


def compute_metrics(
    targets_a: Set[str],
    targets_b: Set[str],
    module_gene_sets: Dict[str, Set[str]],
    module_meta: List[dict],
    universe: Set[str],
) -> pd.DataFrame:
    rows = []
    for meta in module_meta:
        name = meta["module"]
        gset = module_gene_sets[name] & universe
        N = len(gset)
        mA = targets_a & gset
        mB = targets_b & gset
        union = mA | mB
        overlap = mA & mB

        cov_a = safe_div(len(mA), N)
        cov_b = safe_div(len(mB), N)
        cov_union = safe_div(len(union), N)
        incremental = cov_union - max(cov_a, cov_b)
        complementarity = 1.0 - safe_div(len(overlap), len(union)) if union else float("nan")

        rows.append({
            "module": name,
            "module_display": meta["display"],
            "category": meta["category"],
            "n_module_genes": N,
            "n_targets_A": len(mA),
            "n_targets_B": len(mB),
            "n_union": len(union),
            "n_overlap": len(overlap),
            "coverage_A": cov_a,
            "coverage_B": cov_b,
            "coverage_union": cov_union,
            "incremental_coverage": incremental,
            "complementarity": complementarity,
        })

    return pd.DataFrame(rows)


# --- Plotting ---

def plot_incremental_coverage(metrics: pd.DataFrame, title: str, outpath: str) -> None:
    df = metrics.copy()
    df["category"] = pd.Categorical(df["category"], categories=CATEGORY_ORDER, ordered=True)
    df = df.sort_values(["category", "incremental_coverage"], ascending=[True, False])

    colors = {
        "nb_differentiation": "#2ca02c",
        "mycn": "#d62728",
        "nb_undifferentiated": "#7f7f7f",
    }
    bar_colors = [colors.get(c, "#333333") for c in df["category"]]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(range(len(df)), df["incremental_coverage"], color=bar_colors)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["module_display"], rotation=25, ha="right")
    ax.set_ylabel("Incremental coverage")
    ax.set_title(title)

    # Legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=colors[k], label=CATEGORY_LABELS[k]) for k in CATEGORY_ORDER if k in colors]
    ax.legend(handles=handles, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_coverage_bars(metrics: pd.DataFrame, title: str, outpath: str) -> None:
    df = metrics.copy()
    df["category"] = pd.Categorical(df["category"], categories=CATEGORY_ORDER, ordered=True)
    df = df.sort_values(["category", "module_display"])

    x = np.arange(len(df))
    w = 0.25
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - w, df["coverage_A"], width=w, label="miRNA A", alpha=0.8)
    ax.bar(x,     df["coverage_B"], width=w, label="miRNA B", alpha=0.8)
    ax.bar(x + w, df["coverage_union"], width=w, label="Union (A∪B)", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(df["module_display"], rotation=25, ha="right")
    ax.set_ylabel("Fraction of module genes targeted")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


# --- Batch analysis ---

def run_single_pair(
    mirA: str, mirB: str,
    df_ts: pd.DataFrame, spec: TargetSource,
    module_meta: List[dict], module_gene_sets: Dict[str, Set[str]],
    universe: Set[str], outdir: str,
) -> pd.DataFrame:
    targets_a = get_targets(df_ts, spec, mirA)
    targets_b = get_targets(df_ts, spec, mirB)

    if not targets_a:
        print(f"  WARNING: no targets found for {mirA}", file=sys.stderr)
    if not targets_b:
        print(f"  WARNING: no targets found for {mirB}", file=sys.stderr)

    metrics = compute_metrics(targets_a, targets_b, module_gene_sets, module_meta, universe)

    pair_label = f"{mirA.replace('hsa-', '')} + {mirB.replace('hsa-', '')}"
    os.makedirs(outdir, exist_ok=True)
    metrics.to_csv(os.path.join(outdir, "metrics.csv"), index=False)
    plot_incremental_coverage(metrics, f"NB-specific incremental coverage\n{pair_label}",
                               os.path.join(outdir, "incremental_coverage_nb.png"))
    plot_coverage_bars(metrics, f"NB-specific module coverage\n{pair_label}",
                        os.path.join(outdir, "coverage_bars_nb.png"))

    # Save audit trail
    with open(os.path.join(outdir, "module_terms_nb.json"), "w") as f:
        json.dump(module_meta, f, indent=2)

    return metrics


def run_batch(pairs_csv: str, df_ts: pd.DataFrame, spec: TargetSource,
              module_meta: List[dict], module_gene_sets: Dict[str, Set[str]],
              universe: Set[str], outdir: str) -> pd.DataFrame:
    pairs = pd.read_csv(pairs_csv)
    all_metrics = []

    for _, row in pairs.iterrows():
        if "mirA" in pairs.columns:
            mirA = row["mirA"].strip()
            mirB = row["mirB"].strip()
        else:
            mirA = row.iloc[0].strip()
            mirB = row.iloc[1].strip()
        pair_label = f"{mirA}__{mirB}".replace("hsa-", "")
        pair_dir = os.path.join(outdir, "per_pair", pair_label)

        print(f"  {mirA} + {mirB}")
        m = run_single_pair(mirA, mirB, df_ts, spec, module_meta, module_gene_sets, universe, pair_dir)
        m["mirA"] = mirA
        m["mirB"] = mirB
        all_metrics.append(m)

    combined = pd.concat(all_metrics, ignore_index=True)
    combined.to_csv(os.path.join(outdir, "all_pairs_nb_metrics.csv"), index=False)
    return combined


def generate_summary_figure(combined: pd.DataFrame, outdir: str) -> None:
    """Box plot of incremental coverage by module across all pairs."""
    df = combined.copy()
    df["category"] = pd.Categorical(df["category"], categories=CATEGORY_ORDER, ordered=True)

    colors = {
        "nb_differentiation": "#2ca02c",
        "mycn": "#d62728",
        "nb_undifferentiated": "#7f7f7f",
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    modules = df.sort_values("category")["module_display"].unique()

    positions = range(len(modules))
    for i, mod in enumerate(modules):
        subset = df[df["module_display"] == mod]
        cat = subset["category"].iloc[0]
        bp = ax.boxplot(subset["incremental_coverage"].dropna(), positions=[i],
                         widths=0.6, patch_artist=True,
                         boxprops=dict(facecolor=colors.get(cat, "#333"), alpha=0.6))
        # Overlay points
        ax.scatter(np.repeat(i, len(subset)), subset["incremental_coverage"],
                    color=colors.get(cat, "#333"), alpha=0.7, s=20, zorder=3)

    ax.set_xticks(list(positions))
    ax.set_xticklabels(modules, rotation=25, ha="right")
    ax.set_ylabel("Incremental coverage")
    ax.set_title("NB-specific modules: incremental coverage across synergistic pairs")

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=colors[k], alpha=0.6, label=CATEGORY_LABELS[k]) for k in CATEGORY_ORDER]
    ax.legend(handles=handles, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "nb_summary_boxplot.png"), dpi=300)
    plt.close(fig)

    # Statistical tests
    diff_vals = df[df["category"] == "nb_differentiation"]["incremental_coverage"].dropna()
    mycn_vals = df[df["category"] == "mycn"]["incremental_coverage"].dropna()
    undiff_vals = df[df["category"] == "nb_undifferentiated"]["incremental_coverage"].dropna()

    print("\n--- Statistical Summary ---")
    print(f"NB differentiation modules:  mean = {diff_vals.mean():.4f}, median = {diff_vals.median():.4f}, n = {len(diff_vals)}")
    print(f"MYCN targets:                mean = {mycn_vals.mean():.4f}, median = {mycn_vals.median():.4f}, n = {len(mycn_vals)}")
    print(f"NB undifferentiated modules: mean = {undiff_vals.mean():.4f}, median = {undiff_vals.median():.4f}, n = {len(undiff_vals)}")

    stats_rows = []
    for label, vals in [("nb_differentiation", diff_vals), ("mycn", mycn_vals), ("nb_undifferentiated", undiff_vals)]:
        stats_rows.append({
            "category": label,
            "n": int(len(vals)),
            "mean_inc_coverage": float(vals.mean()) if len(vals) else None,
            "median_inc_coverage": float(vals.median()) if len(vals) else None,
            "comparison": None,
            "mannwhitney_U": None,
            "mannwhitney_p_onesided_greater": None,
        })
    if len(diff_vals) > 0 and len(mycn_vals) > 0:
        u, p = stats.mannwhitneyu(diff_vals, mycn_vals, alternative="greater")
        print(f"\nDifferentiation vs MYCN: U = {u:.1f}, p = {p:.4e} (one-sided)")
        stats_rows.append({
            "category": None, "n": None,
            "mean_inc_coverage": None, "median_inc_coverage": None,
            "comparison": "nb_differentiation_vs_mycn",
            "mannwhitney_U": float(u), "mannwhitney_p_onesided_greater": float(p),
        })

    if len(diff_vals) > 0 and len(undiff_vals) > 0:
        u, p = stats.mannwhitneyu(diff_vals, undiff_vals, alternative="greater")
        print(f"Differentiation vs Undifferentiated: U = {u:.1f}, p = {p:.4e} (one-sided)")
        stats_rows.append({
            "category": None, "n": None,
            "mean_inc_coverage": None, "median_inc_coverage": None,
            "comparison": "nb_differentiation_vs_nb_undifferentiated",
            "mannwhitney_U": float(u), "mannwhitney_p_onesided_greater": float(p),
        })

    # --- v16 addition: per-module Mann-Whitney tests ---
    # The figure-rendering composite displays brackets between specific
    # modules (ADRN vs MYCN, ADRN vs Retinoid), not aggregated categories.
    # Emit per-module summaries and the specific tests used on the figure.
    print("\n--- Per-module statistics ---")
    module_vals = {}
    for module_name in df["module_display"].unique():
        vals = df.loc[df["module_display"] == module_name, "incremental_coverage"].dropna()
        module_vals[module_name] = vals
        print(f"  {module_name:30s}  mean = {vals.mean():.4f}  median = {vals.median():.4f}  n = {len(vals)}")
        stats_rows.append({
            "category": None, "n": int(len(vals)),
            "mean_inc_coverage": float(vals.mean()) if len(vals) else None,
            "median_inc_coverage": float(vals.median()) if len(vals) else None,
            "comparison": f"module:{module_name}",
            "mannwhitney_U": None,
            "mannwhitney_p_onesided_greater": None,
        })

    # The two bracket comparisons that appear on the figure:
    bracket_comparisons = [
        ("Adrenergic (ADRN)", "MYCN targets (Wei et al.)", "ADRN_vs_MYCN_targets"),
        ("Adrenergic (ADRN)", "Retinoid response (GO)",   "ADRN_vs_Retinoid"),
    ]
    for mod_a, mod_b, key in bracket_comparisons:
        if mod_a in module_vals and mod_b in module_vals:
            a = module_vals[mod_a]
            b = module_vals[mod_b]
            if len(a) > 0 and len(b) > 0:
                u, p = stats.mannwhitneyu(a, b, alternative="greater")
                print(f"  {key}: U = {u:.1f}, p = {p:.4e} (one-sided ADRN > comparator)")
                stats_rows.append({
                    "category": None, "n": None,
                    "mean_inc_coverage": None, "median_inc_coverage": None,
                    "comparison": key,
                    "mannwhitney_U": float(u),
                    "mannwhitney_p_onesided_greater": float(p),
                })

    stats_df = pd.DataFrame(stats_rows)
    stats_csv_path = os.path.join(outdir, "nb_specific_stats.csv")
    stats_df.to_csv(stats_csv_path, index=False)
    print(f"\nSaved stats: {stats_csv_path}  ({len(stats_df)} rows)")


# --- Main ---

def main():
    parser = argparse.ArgumentParser(description="NB-specific target-space analysis")
    parser.add_argument("--pairs-csv", required=True, help="CSV with miRNA pair columns")
    parser.add_argument("--targetscan-tsv", required=True, help="TargetScan v7.2 predictions TSV")
    parser.add_argument("--sig-dir", required=True, help="Directory with NB signature gene lists")
    parser.add_argument("--outdir", default="nb_specific_results", help="Output directory")
    args = parser.parse_args()

    print("Loading TargetScan data...")
    df_ts = load_targets(args.targetscan_tsv, TARGETSCAN_SPEC)
    universe = get_universe(df_ts, TARGETSCAN_SPEC)
    print(f"  Universe: {len(universe)} genes")

    print("Loading NB-specific modules...")
    module_meta, module_gene_sets = load_modules(args.sig_dir)
    for m in module_meta:
        in_universe = len(module_gene_sets[m["module"]] & universe)
        print(f"  {m['display']}: {m['n_genes']} genes ({in_universe} in TargetScan universe)")

    print("\nRunning batch analysis...")
    os.makedirs(args.outdir, exist_ok=True)
    combined = run_batch(args.pairs_csv, df_ts, TARGETSCAN_SPEC,
                          module_meta, module_gene_sets, universe, args.outdir)

    print("\nGenerating summary figure...")
    generate_summary_figure(combined, args.outdir)

    print(f"\nDone. Results in {args.outdir}/")


if __name__ == "__main__":
    main()
