#!/usr/bin/env python3
"""
miRNA synergy figure builder (v8): coverage + complementarity ONLY

Purpose
-------
Generate manuscript-ready plots supporting a "complementary / distributed regulation"
mechanism for synergistic miRNAs using curated gene sets and pathway/module-level metrics.

Outputs
-------
- metrics.csv
- coverage_bars.png
- incremental_coverage.png
- complementarity.png
- module_terms.json (audit trail: gene-set terms per module)

Notes
-----
- Gene sets are retrieved from Enrichr libraries (cached locally).
- Universe = genes present in the input target database(s); each module is intersected with the universe.
- Metrics quantify module coverage and complementarity; they do NOT claim enriched overlap.

Dependencies
------------
pip install pandas numpy matplotlib requests
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
              max_terms: int = 60) -> List[str]:
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
# Modules (neurite-length phenotype)
# ----------------------------

MODULE_DEFS = [
    # ON-target modules
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
        "max_terms": 80,
        "display": "Neurite outgrowth (on-target)",
        "category": "on_target",
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
        "max_terms": 80,
        "display": "Synapse/maturation (on-target)",
        "category": "on_target",
    },

    # Housekeeping / liabilities / alternative programs
    {
        "module": "Translation / ribosome (housekeeping)",
        "library_group": "reactome",
        "include": [r"\bTranslation\b", r"Peptide chain elongation", r"Translation initiation"],
        "exclude": [r"selenocysteine"],
        "max_terms": 40,
        "display": "Translation/ribosome (housekeeping)",
        "category": "housekeeping",
    },
    {
        "module": "RNA processing / splicing (housekeeping)",
        "library_group": "reactome",
        "include": [r"splic", r"mRNA processing", r"Processing of Capped", r"RNA Polymerase II.*processing"],
        "exclude": [r"tRNA"],
        "max_terms": 40,
        "display": "RNA processing/splicing (housekeeping)",
        "category": "housekeeping",
    },
    {
        "module": "UPR / ER stress (liability)",
        "library_group": "hallmark",
        "include": [r"UNFOLDED_PROTEIN_RESPONSE", r"unfolded protein response"],
        "exclude": [],
        "max_terms": 50,
        "display": "UPR/ER stress (liability)",
        "category": "liability",
    },
    {
        "module": "Apoptosis / p53 / DNA repair (liability)",
        "library_group": "hallmark",
        "include": [r"APOPTOSIS", r"P53_PATHWAY", r"DNA_REPAIR", r"apoptosis", r"dna repair", r"p53"],
        "exclude": [],
        "max_terms": 120,
        "display": "Apoptosis/DNA repair (liability)",
        "category": "liability",
    },
    {
        "module": "EMT (alternate morphology)",
        "library_group": "hallmark",
        "include": [r"EPITHELIAL_MESENCHYMAL_TRANSITION", r"epithelial.*mesenchymal"],
        "exclude": [],
        "max_terms": 50,
        "display": "EMT (alternate morphology)",
        "category": "alternate_morphology",
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
        "max_terms": 80,
        "display": "Trafficking/cytoskeleton (non-specific)",
        "category": "non_specific_morphology",
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

        terms = pick_terms(lib, md["include"], md.get("exclude", []), md.get("max_terms", 60))
        genes = union_gene_sets(lib, terms)

        module_gene_sets[md["module"]] = genes
        module_meta.append({
            "module": md["module"],
            "display": md["display"],
            "category": md["category"],
            "library_group": group,
            "library_name": libname,
            "n_terms": len(terms),
            "terms": terms,
            "n_genes": len(genes),
        })

    return module_meta, module_gene_sets


def safe_div(num: float, den: float) -> float:
    return float(num) / float(den) if den else float("nan")


CATEGORY_ORDER = ["on_target", "liability", "housekeeping", "non_specific_morphology", "alternate_morphology"]


def save_coverage_bars(metrics: pd.DataFrame, outpath: str) -> None:
    df = metrics.copy()
    df["category"] = pd.Categorical(df["category"], categories=CATEGORY_ORDER, ordered=True)
    df = df.sort_values(["category", "module_display"])

    x = np.arange(len(df))
    w = 0.25
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(x - w, df["coverage_A"], width=w, label="miRNA A")
    ax.bar(x,     df["coverage_B"], width=w, label="miRNA B")
    ax.bar(x + w, df["coverage_union"], width=w, label="Union (A∪B)")
    ax.set_ylim(0, min(1.0, max(0.25, float(np.nanmax(df["coverage_union"].values) * 1.15))))
    ax.set_xticks(x)
    ax.set_xticklabels(df["module_display"], rotation=22, ha="right")
    ax.set_ylabel("Fraction of module genes targeted")
    ax.set_title("Pathway/module coverage by each miRNA and their union")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_simple_bar(metrics: pd.DataFrame, col: str, title: str, ylabel: str, outpath: str) -> None:
    df = metrics.copy()
    df["category"] = pd.Categorical(df["category"], categories=CATEGORY_ORDER, ordered=True)
    df = df.sort_values(["category", col], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.bar(df["module_display"], df[col])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=22)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mirA", default="hsa-miR-137", help="miRNA A name")
    ap.add_argument("--mirB", default="hsa-miR-449b-5p", help="miRNA B name")

    ap.add_argument("--mode", choices=["validated", "predicted", "union"], default="union")
    ap.add_argument("--mirtarbase_tsv", default=None)
    ap.add_argument("--targetscan_tsv", default=None)

    ap.add_argument("--cache_dir", default=os.path.join(os.path.expanduser("~"), ".cache", "enrichr_libraries"))
    ap.add_argument("--outdir", default="out_v8")

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
    print(f"Universe (genes in input DBs): {len(universe)}")
    print("Fetching curated gene set libraries and building modules...")

    module_meta, module_gene_sets = build_modules(cache_dir=args.cache_dir)
    meta_path = os.path.join(args.outdir, "module_terms.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(module_meta, f, indent=2)

    rows = []
    for m in module_meta:
        module = m["module"]
        gset = module_gene_sets[module] & universe
        N = len(gset)

        mA = genes_A & gset
        mB = genes_B & gset
        mU = genes_union & gset

        overlap = len(mA & mB)
        union_in_module = len(mU)

        covA = safe_div(len(mA), N)
        covB = safe_div(len(mB), N)
        covU = safe_div(union_in_module, N)

        inc_cov = covU - max(covA, covB) if (not np.isnan(covU) and not np.isnan(covA) and not np.isnan(covB)) else float("nan")
        comp = 1.0 - safe_div(overlap, union_in_module) if union_in_module else float("nan")

        rows.append({
            "module": module,
            "module_display": m["display"],
            "category": m["category"],
            "library_group": m["library_group"],
            "library_name": m["library_name"],
            "n_terms": m["n_terms"],
            "module_genes_total": m["n_genes"],
            "module_genes_in_universe": N,
            "mirA": args.mirA,
            "mirB": args.mirB,
            "mirA_targets_total": len(genes_A),
            "mirB_targets_total": len(genes_B),
            "union_targets_total": len(genes_union),
            "mirA_module_targets": len(mA),
            "mirB_module_targets": len(mB),
            "union_module_targets": union_in_module,
            "overlap": overlap,
            "coverage_A": covA,
            "coverage_B": covB,
            "coverage_union": covU,
            "incremental_coverage": inc_cov,
            "complementarity_index": comp,
        })

        if N:
            print(f"[{module}] N={N} |A|={len(mA)} |B|={len(mB)} U={union_in_module} ov={overlap} "
                  f"covA={covA:.3f} covB={covB:.3f} covU={covU:.3f} inc={inc_cov:.3f} comp={comp:.3f}")
        else:
            print(f"[{module}] N=0")

    metrics = pd.DataFrame(rows)
    metrics_path = os.path.join(args.outdir, "metrics.csv")
    metrics.to_csv(metrics_path, index=False)

    cov_path = os.path.join(args.outdir, "coverage_bars.png")
    save_coverage_bars(metrics, cov_path)

    inc_path = os.path.join(args.outdir, "incremental_coverage.png")
    save_simple_bar(metrics.fillna(0), "incremental_coverage",
                    "Incremental module coverage gained by combining miRNAs",
                    "Incremental coverage (A∪B − max(A,B))", inc_path)

    comp_path = os.path.join(args.outdir, "complementarity.png")
    save_simple_bar(metrics.fillna(0), "complementarity_index",
                    "Complementarity of targeting within modules",
                    "Complementarity (1 − overlap / union_module_targets)", comp_path)

  # print("\nWrote:")
  # print(f"  {metrics_path}")
  # print(f"  {cov_path}")
  # print(f"  {inc_path}")
  # print(f"  {comp_path}")
  # print(f"  {meta_path}  (terms used per module)")


if __name__ == "__main__":
    main()
