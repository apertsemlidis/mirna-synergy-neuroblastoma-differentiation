#!/usr/bin/env python3
"""
miRNA synergy figure builder (v7): curated gene sets + correct null + complementarity/coverage metrics

What v7 adds on top of v6b:
- Keeps curated gene sets via Enrichr + module-conditioned hypergeometric p-values (correct baseline)
- Fixes Hallmark module term matching (robust includes)
- Adds metrics that support a "complementarity / dilution" narrative even when overlap isn't enriched:
    * coverage_A = |A∩M| / |M|
    * coverage_B = |B∩M| / |M|
    * coverage_union = |(A∪B)∩M| / |M|      (phenotype-relevant coverage)
    * incremental_coverage = coverage_union - max(coverage_A, coverage_B)
    * complementarity_index = 1 - (|A∩B∩M| / |(A∪B)∩M|)  (0=all shared, 1=fully non-overlapping)
    * off_target_burden = |(A∪B)∩M_off| / |A∪B|  for user-defined OFF modules grouping
- Adds summary plots:
    * coverage_bars.png (coverage_A/B/union per module)
    * incremental_coverage.png
    * complementarity.png

This lets you say (if supported by data):
- The combination increases *coverage* of neurite-length modules (complementary targeting),
  even if direct overlap is not enriched.
- Off-target modules show lower incremental coverage and/or higher complementarity (dilution) vs on-target.

Deps:
  pip install pandas numpy matplotlib matplotlib-venn requests
Optional:
  pip install scipy   (faster hypergeom)
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import json
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import requests

# ----------------------------
# Input table specs
# ----------------------------

@dataclass(frozen=True)
class TargetSource:
    name: str
    mir_col: str
    gene_col: str

MIRTARBASE_SPEC = TargetSource("miRTarBase", "miRNA", "Target Gene")
TARGETSCAN_SPEC = TargetSource("TargetScan", "miRNA", "GeneSymbol")


def load_table(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df


def extract_targets(df: pd.DataFrame, spec: TargetSource, mirna: str) -> Set[str]:
    mirna = str(mirna).strip()
    if spec.mir_col not in df.columns or spec.gene_col not in df.columns:
        raise ValueError(
            f"{spec.name} table must have columns '{spec.mir_col}' and '{spec.gene_col}'. "
            f"Found: {list(df.columns)}"
        )
    sub = df[df[spec.mir_col].astype(str).str.strip() == mirna]
    genes = set(
        g.strip() for g in sub[spec.gene_col].dropna().astype(str).tolist()
        if g.strip() and g.strip().upper() not in {"NA", "N/A", "NULL"}
    )
    return genes


def extract_universe(df: pd.DataFrame, spec: TargetSource) -> Set[str]:
    if spec.gene_col not in df.columns:
        raise ValueError(f"Missing gene column '{spec.gene_col}' in {spec.name} table")
    genes = set(
        g.strip() for g in df[spec.gene_col].dropna().astype(str).tolist()
        if g.strip() and g.strip().upper() not in {"NA", "N/A", "NULL"}
    )
    return genes


# ----------------------------
# Enrichr library retrieval
# ----------------------------

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"


def _safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def fetch_enrichr_library(library_name: str, cache_dir: str) -> Dict[str, List[str]]:
    _safe_mkdir(cache_dir)
    cache_path = os.path.join(cache_dir, f"{library_name}.txt")

    if os.path.exists(cache_path):
        raw = open(cache_path, "r", encoding="utf-8").read()
    else:
        url = f"{ENRICHR_BASE}/geneSetLibrary"
        params = {"mode": "text", "libraryName": library_name}
        resp = requests.get(url, params=params, timeout=60)
        if resp.status_code != 200:
            raise RuntimeError(f"Failed to fetch Enrichr library '{library_name}' (HTTP {resp.status_code}).")
        raw = resp.text
        open(cache_path, "w", encoding="utf-8").write(raw)

    lib: Dict[str, List[str]] = {}
    for line in raw.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        term = parts[0].strip()
        genes = [g.strip() for g in parts[1:] if g.strip()]
        if term and genes:
            lib[term] = genes
    if not lib:
        raise RuntimeError(f"Parsed 0 gene sets from Enrichr library '{library_name}'.")
    return lib


def pick_terms(lib: Dict[str, List[str]], include: List[str], exclude: Optional[List[str]] = None,
              max_terms: int = 25) -> List[str]:
    exclude = exclude or []
    inc_res = [re.compile(p, re.IGNORECASE) for p in include]
    exc_res = [re.compile(p, re.IGNORECASE) for p in exclude]

    hits = []
    for term in lib.keys():
        if any(r.search(term) for r in inc_res) and not any(r.search(term) for r in exc_res):
            hits.append(term)
    return sorted(hits)[:max_terms]


def union_gene_sets(lib: Dict[str, List[str]], terms: List[str]) -> Set[str]:
    genes: Set[str] = set()
    for t in terms:
        for g in lib.get(t, []):
            if g and g.upper() not in {"NA", "N/A", "NULL"}:
                genes.add(g)
    return genes


DEFAULT_LIB_CANDIDATES = {
    "reactome": ["Reactome_2022", "Reactome_2023", "Reactome_2016"],
    "hallmark": ["MSigDB_Hallmark_2023", "MSigDB_Hallmark_2022", "MSigDB_Hallmark_2020"],
    "go_bp": ["GO_Biological_Process_2023", "GO_Biological_Process_2021"],
}


# ----------------------------
# v7 module definitions
# ----------------------------
# We keep Reactome for neuronal pathways and use GO_BP as a more specific neurite/axon/synapse option.
# Hallmark matching is made robust by matching either HALLMARK_* or plain phrases.

MODULE_DEFS = [
    # ON-target modules (neurite-length phenotype)
    {
        "module": "Neurite outgrowth & projection morphogenesis",
        "library_group": "go_bp",
        "include": [
            r"neurite outgrowth",
            r"neuron projection (development|morphogenesis|organization)",
            r"axonogenesis",
            r"axon (development|extension|guidance|outgrowth|regeneration|pathfinding)",
            r"dendrite (development|morphogenesis|growth)",
            r"growth cone",
        ],
        "exclude": [r"immune", r"viral", r"interferon"],
        "max_terms": 40,
        "display": "Neurite outgrowth (on-target coverage)",
        "group": "on_target",
    },
    {
        "module": "Synapse formation & neuronal maturation",
        "library_group": "go_bp",
        "include": [
            r"synapse (assembly|organization)",
            r"synaptic signaling",
            r"chemical synaptic transmission",
            r"neurotransmitter (secretion|transport|release)",
        ],
        "exclude": [],
        "max_terms": 40,
        "display": "Synapse/maturation (on-target coverage)",
        "group": "on_target",
    },

    # OFF-target / controls
    {
        "module": "Translation / ribosome (housekeeping)",
        "library_group": "reactome",
        "include": [r"\bTranslation\b", r"Peptide chain elongation", r"Translation initiation"],
        "exclude": [r"selenocysteine"],
        "max_terms": 25,
        "display": "Translation/ribosome (off-target burden)",
        "group": "off_target",
    },
    {
        "module": "RNA processing / splicing (housekeeping)",
        "library_group": "reactome",
        "include": [r"splic", r"mRNA processing", r"Processing of Capped", r"RNA Polymerase II.*processing"],
        "exclude": [r"tRNA"],
        "max_terms": 25,
        "display": "RNA processing/splicing (off-target burden)",
        "group": "off_target",
    },
    {
        "module": "UPR / ER stress (liability)",
        "library_group": "hallmark",
        "include": [r"UNFOLDED_PROTEIN_RESPONSE", r"unfolded protein response"],
        "exclude": [],
        "max_terms": 10,
        "display": "UPR/ER stress (off-target liability)",
        "group": "off_target",
    },
    {
        "module": "Apoptosis / p53 / DNA repair (liability)",
        "library_group": "hallmark",
        "include": [r"APOPTOSIS", r"P53_PATHWAY", r"DNA_REPAIR", r"apoptosis", r"dna repair", r"p53"],
        "exclude": [],
        "max_terms": 20,
        "display": "Apoptosis/DNA repair (off-target liability)",
        "group": "off_target",
    },
    {
        "module": "EMT (alternate morphology)",
        "library_group": "hallmark",
        "include": [r"EPITHELIAL_MESENCHYMAL_TRANSITION", r"epithelial.*mesenchymal"],
        "exclude": [],
        "max_terms": 10,
        "display": "EMT (off-target morphology)",
        "group": "off_target",
    },
    {
        "module": "Trafficking / cytoskeleton (non-specific morphology)",
        "library_group": "go_bp",
        "include": [
            r"vesicle[- ]mediated transport",
            r"endocyt",
            r"exocyt",
            r"golgi",
            r"endosome",
            r"lysosome",
            r"actin cytoskeleton",
            r"microtubule",
        ],
        "exclude": [r"axon", r"synapse", r"neurite", r"neuron projection"],
        "max_terms": 40,
        "display": "Trafficking/cytoskeleton (off-target morphology)",
        "group": "off_target",
    },
]


def resolve_library(group: str, cache_dir: str) -> Tuple[str, Dict[str, List[str]]]:
    candidates = DEFAULT_LIB_CANDIDATES.get(group, [])
    if not candidates:
        raise ValueError(f"Unknown library group: {group}")
    last_err = None
    for libname in candidates:
        try:
            lib = fetch_enrichr_library(libname, cache_dir=cache_dir)
            return libname, lib
        except Exception as e:
            last_err = e
    raise RuntimeError(f"Could not fetch any Enrichr library for group '{group}'. Last error: {last_err}")


def build_modules(cache_dir: str) -> Tuple[List[dict], Dict[str, Set[str]]]:
    libs: Dict[str, Tuple[str, Dict[str, List[str]]]] = {}
    module_gene_sets: Dict[str, Set[str]] = {}
    module_meta: List[dict] = []

    for md in MODULE_DEFS:
        group = md["library_group"]
        if group not in libs:
            libs[group] = resolve_library(group, cache_dir=cache_dir)
        libname, lib = libs[group]

        terms = pick_terms(lib, md["include"], md.get("exclude", []), md.get("max_terms", 25))
        genes = union_gene_sets(lib, terms)

        module_gene_sets[md["module"]] = genes
        module_meta.append({
            "module": md["module"],
            "display": md["display"],
            "category": md.get("group", "uncategorized"),
            "library_group": group,
            "library_name": libname,
            "n_terms": len(terms),
            "terms": terms,
            "n_genes": len(genes),
        })

    return module_meta, module_gene_sets


# ----------------------------
# Stats
# ----------------------------

def jaccard(a: Set[str], b: Set[str]) -> float:
    u = a | b
    return len(a & b) / len(u) if u else float("nan")


def overlap_coeff(a: Set[str], b: Set[str]) -> float:
    denom = min(len(a), len(b))
    return len(a & b) / denom if denom else float("nan")


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    # P(X >= k) for Hypergeometric(N, K, n)
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(k - 1, N, K, n))
    except Exception:
        import math
        if k <= 0:
            return 1.0
        kmax = min(K, n)
        if k > kmax:
            return 0.0

        def logC(a: int, b: int) -> float:
            return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

        denom = logC(N, n)
        s = 0.0
        for x in range(k, kmax + 1):
            s += math.exp(logC(K, x) + logC(N - K, n - x) - denom)
        return float(min(max(s, 0.0), 1.0))


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.empty(n, dtype=float)
    for i, idx in enumerate(order):
        ranked[idx] = pvals[idx] * n / (i + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(ranked, 0, 1)


# ----------------------------
# Plotting
# ----------------------------

def venn_or_text(ax, a: Set[str], b: Set[str], title: str, labelA: str, labelB: str):
    a_only, b_only, inter = len(a - b), len(b - a), len(a & b)
    uni = len(a | b)
    J = jaccard(a, b)
    OC = overlap_coeff(a, b)

    ax.set_title(title, fontsize=12, pad=10)

    if uni == 0:
        ax.text(0.5, 0.5, "No module targets\nfor either miRNA",
                ha="center", va="center", fontsize=11)
        ax.set_axis_off()
        return

    venn2(subsets=(a_only, b_only, inter), set_labels=(labelA, labelB), ax=ax)
    ax.text(
        0.5, 1.01,
        f"overlap={inter}, union={uni}, J={J:.3f}, OC={OC:.3f}",
        ha="center", va="bottom", transform=ax.transAxes, fontsize=10
    )
    ax.set_axis_off()


def save_venn_grid(module_sets: Dict[str, Tuple[Set[str], Set[str]]],
                   title_map: Dict[str, str],
                   outpath: str,
                   labelA: str, labelB: str) -> None:
    n = len(module_sets)
    cols = 2
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(12, 5.0 * rows))
    if rows == 1:
        axes = np.array(axes).reshape(1, -1)

    for i, (module, (a, b)) in enumerate(module_sets.items()):
        r, c = divmod(i, cols)
        venn_or_text(axes[r, c], a, b, title_map.get(module, module), labelA, labelB)

    for i in range(n, rows * cols):
        r, c = divmod(i, cols)
        axes[r, c].axis("off")

    fig.suptitle("Curated gene-set overlap stratified by neurite-length modules (v7)",
                 fontsize=14, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_grouped_coverage(metrics: pd.DataFrame, outpath: str) -> None:
    df = metrics.copy()
    x = np.arange(len(df))
    w = 0.25
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(x - w, df["coverage_A"], width=w, label="miRNA A coverage")
    ax.bar(x,     df["coverage_B"], width=w, label="miRNA B coverage")
    ax.bar(x + w, df["coverage_union"], width=w, label="Union coverage (A∪B)")
    ax.set_ylim(0, 1)
    ax.set_xticks(x)
    ax.set_xticklabels(df["module_display"], rotation=22, ha="right")
    ax.set_ylabel("Fraction of module genes targeted")
    ax.set_title("Module coverage: complementarity vs redundancy")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_simple_bar(metrics: pd.DataFrame, col: str, title: str, ylabel: str, outpath: str) -> None:
    df = metrics.copy().sort_values(col, ascending=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df["module_display"], df[col])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=22)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mirA", default="hsa-miR-137", help="miRNA A name")
    ap.add_argument("--mirB", default="hsa-miR-449b-5p", help="miRNA B name")

    ap.add_argument("--mode", choices=["validated", "predicted", "union"], default="union")
    ap.add_argument("--mirtarbase_tsv", default=None)
    ap.add_argument("--targetscan_tsv", default=None)

    ap.add_argument("--cache_dir", default=os.path.join(os.path.expanduser("~"), ".cache", "enrichr_libraries"))
    ap.add_argument("--outdir", default="out_v7")

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sources: List[Tuple[pd.DataFrame, TargetSource]] = []
    universes: List[Set[str]] = []

    if args.mode in ["validated", "union"]:
        if not args.mirtarbase_tsv:
            print("ERROR: --mirtarbase_tsv required for validated/union mode", file=sys.stderr)
            sys.exit(2)
        df_mt = load_table(args.mirtarbase_tsv)
        sources.append((df_mt, MIRTARBASE_SPEC))
        universes.append(extract_universe(df_mt, MIRTARBASE_SPEC))

    if args.mode in ["predicted", "union"]:
        if not args.targetscan_tsv:
            print("ERROR: --targetscan_tsv required for predicted/union mode", file=sys.stderr)
            sys.exit(2)
        df_ts = load_table(args.targetscan_tsv)
        sources.append((df_ts, TARGETSCAN_SPEC))
        universes.append(extract_universe(df_ts, TARGETSCAN_SPEC))

    universe = set().union(*universes)

    genes_A = set().union(*(extract_targets(df, spec, args.mirA) for df, spec in sources))
    genes_B = set().union(*(extract_targets(df, spec, args.mirB) for df, spec in sources))
    genes_union = genes_A | genes_B

    print(f"{args.mirA}: {len(genes_A)} targets")
    print(f"{args.mirB}: {len(genes_B)} targets")
    print(f"Universe (all genes in input DBs): {len(universe)}")

    print("Fetching curated gene set libraries and building modules...")
    module_meta, module_gene_sets = build_modules(cache_dir=args.cache_dir)

    meta_path = os.path.join(args.outdir, "module_terms.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(module_meta, f, indent=2)

    title_map = {m["module"]: m["display"] for m in module_meta}

    labelA = "miR-137" if "137" in args.mirA else "miRNA A"
    labelB = "miR-449b" if "449b" in args.mirB.lower() else "miRNA B"

    rows = []
    module_sets_for_plot: Dict[str, Tuple[Set[str], Set[str]]] = {}

    # For aggregate on/off burden stats
    on_union_hits = 0
    off_union_hits = 0

    for m in module_meta:
        module = m["module"]
        category = m["category"]
        gset = (module_gene_sets[module] & universe)
        N = len(gset)

        mA = genes_A & gset
        mB = genes_B & gset
        mU = genes_union & gset

        inter = len(mA & mB)
        uni = len(mA | mB)

        jac = jaccard(mA, mB)
        oc = overlap_coeff(mA, mB)

        if N == 0 or len(mA) == 0 or len(mB) == 0:
            p = float("nan")
        else:
            p = hypergeom_sf(inter, N=N, K=len(mA), n=len(mB))

        covA = len(mA) / N if N else float("nan")
        covB = len(mB) / N if N else float("nan")
        covU = len(mU) / N if N else float("nan")
        inc_cov = covU - max(covA, covB) if (not np.isnan(covU) and not np.isnan(covA) and not np.isnan(covB)) else float("nan")
        comp = 1.0 - (inter / len(mU)) if len(mU) else float("nan")

        if category == "on_target":
            on_union_hits += len(mU)
        elif category == "off_target":
            off_union_hits += len(mU)

        module_sets_for_plot[module] = (mA, mB)

        rows.append({
            "module": module,
            "module_display": m["display"],
            "category": category,
            "library_group": m["library_group"],
            "library_name": m["library_name"],
            "n_terms": m["n_terms"],
            "module_genes_total": m["n_genes"],
            "module_genes_in_universe": N,
            "mirA": args.mirA,
            "mirB": args.mirB,
            "mirA_module_targets": len(mA),
            "mirB_module_targets": len(mB),
            "union_module_targets": len(mU),
            "overlap": inter,
            "union": uni,
            "jaccard": jac,
            "overlap_coeff": oc,
            "p_hypergeom_ge": p,
            "coverage_A": covA,
            "coverage_B": covB,
            "coverage_union": covU,
            "incremental_coverage": inc_cov,
            "complementarity_index": comp,
        })

        p_txt = f"{p:.4g}" if not np.isnan(p) else "nan"
        print(f"[{module}] N={N} |A|={len(mA)} |B|={len(mB)} U={len(mU)} "
              f"ov={inter} J={jac:.3f} OC={oc:.3f} p={p_txt} "
              f"covU={covU:.3f}" if N else f"[{module}] N=0")

    metrics = pd.DataFrame(rows)

    pvals = metrics["p_hypergeom_ge"].values.astype(float)
    mask = ~np.isnan(pvals)
    fdr = np.full_like(pvals, np.nan, dtype=float)
    if mask.sum() > 0:
        fdr[mask] = bh_fdr(pvals[mask])
    metrics["fdr_bh"] = fdr

    # Aggregate on/off burden fraction of union targets that land in these modules
    total_union = len(genes_union) if len(genes_union) else np.nan
    metrics.attrs["on_target_union_hits"] = on_union_hits
    metrics.attrs["off_target_union_hits"] = off_union_hits
    metrics.attrs["on_target_burden_fraction"] = (on_union_hits / total_union) if total_union else np.nan
    metrics.attrs["off_target_burden_fraction"] = (off_union_hits / total_union) if total_union else np.nan

    metrics_path = os.path.join(args.outdir, "metrics.csv")
    metrics.to_csv(metrics_path, index=False)

    # Plots
    venn_path = os.path.join(args.outdir, "venn_panels.png")
    save_venn_grid(module_sets_for_plot, title_map, venn_path, labelA, labelB)

    cov_path = os.path.join(args.outdir, "coverage_bars.png")
    save_grouped_coverage(metrics, cov_path)

    inc_path = os.path.join(args.outdir, "incremental_coverage.png")
    save_simple_bar(metrics.fillna(0), "incremental_coverage",
                    "Incremental module coverage gained by combining miRNAs",
                    "Incremental coverage (A∪B − max(A,B))", inc_path)

    comp_path = os.path.join(args.outdir, "complementarity.png")
    save_simple_bar(metrics.fillna(0), "complementarity_index",
                    "Complementarity of targeting within modules",
                    "Complementarity index (1 − overlap / union_module_targets)", comp_path)

    print("\nWrote:")
    print(f"  {metrics_path}")
    print(f"  {venn_path}")
    print(f"  {cov_path}")
    print(f"  {inc_path}")
    print(f"  {comp_path}")
    print(f"  {meta_path}  (terms used per module)")

    # Print aggregate burden numbers for quick interpretation
    print("\nAggregate (counting union targets that fall into each curated module):")
    print(f"  union targets total: {len(genes_union)}")
    print(f"  on-target union hits:  {on_union_hits}")
    print(f"  off-target union hits: {off_union_hits}")
    if len(genes_union):
        print(f"  on-target fraction:  {on_union_hits/len(genes_union):.3f}")
        print(f"  off-target fraction: {off_union_hits/len(genes_union):.3f}")


if __name__ == "__main__":
    main()
