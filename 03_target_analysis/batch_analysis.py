#!/usr/bin/env python3
"""
Batch analysis for miRNA-pair "coverage + complementarity" using miRTarBase validated targets.

Adapted from mirna_synergy_batch_v1_2.py to work with miRTarBase instead of TargetScan.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import json
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Optional, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

# ----------------------------
# Enrichr library retrieval + caching (same as TargetScan version)
# ----------------------------

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

DEFAULT_LIB_CANDIDATES = {
    "reactome": ["Reactome_2022", "Reactome_2023", "Reactome_2016"],
    "hallmark": ["MSigDB_Hallmark_2023", "MSigDB_Hallmark_2022", "MSigDB_Hallmark_2020"],
    "go_bp": ["GO_Biological_Process_2023", "GO_Biological_Process_2021"],
}

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

# ----------------------------
# Modules (same as TargetScan version)
# ----------------------------

MODULE_DEFS = [
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
    {
        "module": "Translation / ribosome (housekeeping)",
        "library_group": "reactome",
        "include": [r"\bTranslation\b", r"Peptide chain elongation", r"Translation initiation"],
        "exclude": [r"selenocysteine"],
        "max_terms": 40,
        "display": "Translation/ribosome (housekeeping)",
        "category": "control",
    },
    {
        "module": "RNA processing / splicing (housekeeping)",
        "library_group": "reactome",
        "include": [r"splic", r"mRNA processing", r"Processing of Capped", r"RNA Polymerase II.*processing"],
        "exclude": [r"tRNA"],
        "max_terms": 40,
        "display": "RNA processing/splicing (housekeeping)",
        "category": "control",
    },
    {
        "module": "UPR / ER stress (liability)",
        "library_group": "hallmark",
        "include": [r"UNFOLDED_PROTEIN_RESPONSE", r"unfolded protein response"],
        "exclude": [],
        "max_terms": 50,
        "display": "UPR/ER stress (control)",
        "category": "control",
    },
    {
        "module": "Apoptosis / p53 / DNA repair (liability)",
        "library_group": "hallmark",
        "include": [r"APOPTOSIS", r"P53_PATHWAY", r"DNA_REPAIR", r"apoptosis", r"dna repair", r"p53"],
        "exclude": [],
        "max_terms": 120,
        "display": "Apoptosis/DNA repair (control)",
        "category": "control",
    },
    {
        "module": "EMT (alternate morphology)",
        "library_group": "hallmark",
        "include": [r"EPITHELIAL_MESENCHYMAL_TRANSITION", r"epithelial.*mesenchymal"],
        "exclude": [],
        "max_terms": 50,
        "display": "EMT (control)",
        "category": "control",
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
        "display": "Trafficking/cytoskeleton (control)",
        "category": "control",
    },
]

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

# ----------------------------
# miRTarBase-specific helpers
# ----------------------------

def load_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df

def extract_universe(df: pd.DataFrame, gene_col: str) -> Set[str]:
    genes = set(
        g.strip() for g in df[gene_col].dropna().astype(str).tolist()
        if g.strip() and g.strip().upper() not in {"NA", "N/A", "NULL"}
    )
    return genes

def safe_div(num: float, den: float) -> float:
    return float(num) / float(den) if den else float("nan")

# miRNA normalization
def normalize_mirna(name: str) -> str:
    n = str(name).strip()
    if not n:
        return n
    n = n.replace("–", "-").replace("—", "-")
    n = re.sub(r"^hsa[-_]", "hsa-", n, flags=re.IGNORECASE)
    n = re.sub(r"^mir", "miR", n, flags=re.IGNORECASE)
    if not re.search(r"miR-", n, flags=re.IGNORECASE):
        n = re.sub(r"^miR", "miR-", n, flags=re.IGNORECASE)
    if not re.match(r"^[a-z]{3}-", n, flags=re.IGNORECASE):
        n = "hsa-" + n
    return n

DEFAULT_ALIASES = {
    "hsa-miR-124": "hsa-miR-124-3p",
    "hsa-miR-34b": "hsa-miR-34b-5p",
    "hsa-miR-450b": "hsa-miR-450b-5p",
    "hsa-miR-449b": "hsa-miR-449b-5p",
    "hsa-miR-137": "hsa-miR-137-3p",
}

def apply_alias(m: str, aliases: Dict[str, str]) -> str:
    return aliases.get(m, m)

def extract_targets(df: pd.DataFrame, mir_col: str, gene_col: str, mirna_name: str) -> Set[str]:
    """Extract targets for a given miRNA from miRTarBase table"""
    mirna_name = str(mirna_name).strip()
    sub = df[df[mir_col].astype(str).str.strip() == mirna_name]
    genes = set(
        g.strip() for g in sub[gene_col].dropna().astype(str).tolist()
        if g.strip() and g.strip().upper() not in {"NA", "N/A", "NULL"}
    )
    return genes

def slug_pair(a: str, b: str) -> str:
    s = f"{a}__{b}"
    s = re.sub(r"[^A-Za-z0-9_\-\.]+", "_", s)
    return s

# ----------------------------
# Per-pair analysis
# ----------------------------

def compute_pair_metrics(
    mirA_name: str,
    mirB_name: str,
    label: str,
    df_mt: pd.DataFrame,
    mir_col: str,
    gene_col: str,
    module_meta: List[dict],
    module_gene_sets: Dict[str, Set[str]],
    universe: Set[str],
    outdir_pair: str,
) -> Tuple[pd.DataFrame, dict]:
    os.makedirs(outdir_pair, exist_ok=True)

    genes_A = extract_targets(df_mt, mir_col, gene_col, mirA_name)
    genes_B = extract_targets(df_mt, mir_col, gene_col, mirB_name)
    genes_union = genes_A | genes_B

    with open(os.path.join(outdir_pair, "module_terms.json"), "w", encoding="utf-8") as f:
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
        inc = covU - max(covA, covB) if (not np.isnan(covU) and not np.isnan(covA) and not np.isnan(covB)) else float("nan")
        comp = 1.0 - safe_div(overlap, union_in_module) if union_in_module else float("nan")

        rows.append({
            "label": label,
            "mirA_name": mirA_name,
            "mirB_name": mirB_name,
            "module": module,
            "module_display": m["display"],
            "category": m["category"],
            "module_genes_in_universe": N,
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
            "incremental_coverage": inc,
            "complementarity_index": comp,
        })

    metrics = pd.DataFrame(rows)
    metrics.to_csv(os.path.join(outdir_pair, "metrics.csv"), index=False)

    on = metrics[metrics["category"] == "on_target"]
    ctrl = metrics[metrics["category"] == "control"]

    summary = {
        "label": label,
        "mirA_name": mirA_name,
        "mirB_name": mirB_name,
        "mirA_targets_total": int(metrics["mirA_targets_total"].iloc[0]),
        "mirB_targets_total": int(metrics["mirB_targets_total"].iloc[0]),
        "union_targets_total": int(metrics["union_targets_total"].iloc[0]),
        "Inc_on": float(on["incremental_coverage"].mean(skipna=True)) if len(on) else float("nan"),
        "Inc_ctrl": float(ctrl["incremental_coverage"].mean(skipna=True)) if len(ctrl) else float("nan"),
        "Comp_on": float(on["complementarity_index"].mean(skipna=True)) if len(on) else float("nan"),
        "Comp_ctrl": float(ctrl["complementarity_index"].mean(skipna=True)) if len(ctrl) else float("nan"),
    }
    summary["NBS"] = summary["Inc_on"] - summary["Inc_ctrl"] if (not np.isnan(summary["Inc_on"]) and not np.isnan(summary["Inc_ctrl"])) else float("nan")
    return metrics, summary

# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs_csv", required=True, help="CSV with mirA, mirB, label(optional)")
    ap.add_argument("--mirtarbase_tsv", required=True, help="miRTarBase TSV")
    ap.add_argument("--outdir", default="batch_mirtarbase")
    ap.add_argument("--cache_dir", default=os.path.join(os.path.expanduser("~"), ".cache", "enrichr_libraries"))
    ap.add_argument("--aliases_json", default=None, help="Optional JSON dict mapping miRNA aliases")
    ap.add_argument("--mir_col", default="miRNA", help="miRNA column name in miRTarBase TSV")
    ap.add_argument("--gene_col", default="Target Gene", help="gene symbol column name")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    per_pair_dir = os.path.join(args.outdir, "per_pair")
    os.makedirs(per_pair_dir, exist_ok=True)

    aliases = dict(DEFAULT_ALIASES)
    if args.aliases_json:
        with open(args.aliases_json, "r", encoding="utf-8") as f:
            user_aliases = json.load(f)
        if not isinstance(user_aliases, dict):
            raise ValueError("--aliases_json must contain a JSON object/dict")
        aliases.update({str(k): str(v) for k, v in user_aliases.items()})

    pairs = pd.read_csv(args.pairs_csv, dtype=str).fillna("")
    if "mirA" not in pairs.columns or "mirB" not in pairs.columns:
        raise ValueError(f"--pairs_csv must contain mirA, mirB. Found: {list(pairs.columns)}")
    if "label" not in pairs.columns:
        pairs["label"] = ""

    df_mt = load_table(args.mirtarbase_tsv)
    if args.mir_col not in df_mt.columns or args.gene_col not in df_mt.columns:
        raise ValueError(f"miRTarBase TSV missing required columns. Need '{args.mir_col}' and '{args.gene_col}'. Found: {list(df_mt.columns)}")

    # Filter to human only
    df_mt = df_mt[df_mt[args.mir_col].astype(str).str.lower().str.startswith("hsa-")].copy()

    present_ids = set(df_mt[args.mir_col].dropna().astype(str).str.strip().unique().tolist())
    universe = extract_universe(df_mt, args.gene_col)

    print(f"Loaded {len(present_ids)} unique human miRNAs from miRTarBase")
    print(f"Universe: {len(universe)} genes\n")

    module_meta, module_gene_sets = build_modules(cache_dir=args.cache_dir)

    summaries = []
    unmatched = []

    for _, row in pairs.iterrows():
        rawA, rawB = row["mirA"], row["mirB"]
        label = row.get("label", "")

        normA = apply_alias(normalize_mirna(rawA), aliases)
        normB = apply_alias(normalize_mirna(rawB), aliases)

        if normA not in present_ids:
            unmatched.append({"raw": rawA, "normalized": normA, "reason": "not in miRTarBase"})
            print(f"WARNING: miRNA '{rawA}' (normalized '{normA}') not found in miRTarBase", file=sys.stderr)
            continue
        if normB not in present_ids:
            unmatched.append({"raw": rawB, "normalized": normB, "reason": "not in miRTarBase"})
            print(f"WARNING: miRNA '{rawB}' (normalized '{normB}') not found in miRTarBase", file=sys.stderr)
            continue

        print(f"Processing: {normA} + {normB}")

        pair_slug = slug_pair(normA, normB)
        out_pair = os.path.join(per_pair_dir, pair_slug)

        _, summary = compute_pair_metrics(
            mirA_name=normA,
            mirB_name=normB,
            label=label,
            df_mt=df_mt,
            mir_col=args.mir_col,
            gene_col=args.gene_col,
            module_meta=module_meta,
            module_gene_sets=module_gene_sets,
            universe=universe,
            outdir_pair=out_pair,
        )
        summary["pair_slug"] = pair_slug
        summary["mirA_raw"] = rawA
        summary["mirB_raw"] = rawB
        summaries.append(summary)

    if unmatched:
        pd.DataFrame(unmatched).to_csv(os.path.join(args.outdir, "unmatched_mirnas.csv"), index=False)

    summary_df = pd.DataFrame(summaries)
    summary_path = os.path.join(args.outdir, "pairs_summary.csv")
    summary_df.to_csv(summary_path, index=False)

    print(f"\nWrote:")
    print(f"  {summary_path}")
    print(f"  Successfully processed {len(summary_df)}/{len(pairs)} pairs")
    if unmatched:
        print(f"  {len(unmatched)} miRNAs not matched (see unmatched_mirnas.csv)")

if __name__ == "__main__":
    main()
