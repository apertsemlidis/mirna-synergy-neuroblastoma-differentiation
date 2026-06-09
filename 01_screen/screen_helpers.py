#!/usr/bin/env python3
"""
Shared helpers for the screen analysis tier.

Used by the three screen-processing steps (plate loading, heatmap tables,
and superhit calling).

ORDERED_MIRS is inlined here; HSA_df / abs_heatmap_df are self-contained
with the same float32 cast, singlet/combo split, and min vs max for the
HSA branch.
"""

import numpy as np
import pandas as pd

# Canonical miRNA ordering — clustermap-derived dendrogram order.
ORDERED_MIRS = [
    "hsa-miR-340-5p",
    "hsa-miR-124-3p",
    "hsa-miR-506-3p",
    "hsa-miR-19a-3p",
    "hsa-miR-19b-3p",
    "hsa-miR-106b-5p",
    "hsa-miR-20a-5p",
    "hsa-miR-106a-5p",
    "hsa-miR-17-5p",
    "hsa-miR-20b-5p",
    "hsa-miR-93-5p",
    "hsa-miR-27a-3p",
    "hsa-miR-27b-3p",
    "hsa-miR-137",
    "hsa-miR-429",
    "hsa-miR-200b-3p",
    "hsa-miR-200c-3p",
    "hsa-miR-92a-3p",
    "hsa-miR-25-3p",
    "hsa-miR-363-3p",
    "hsa-miR-128-3p",
    "hsa-miR-3714",
    "hsa-miR-204-5p",
    "hsa-miR-211-5p",
    "hsa-miR-200a-3p",
    "hsa-miR-135b-5p",
    "hsa-miR-449b-5p",
    "hsa-miR-34a-5p",
    "hsa-miR-449a",
    "hsa-miR-103a-3p",
    "hsa-miR-107",
    "hsa-miR-2110",
    "hsa-miR-452-5p",
    "hsa-miR-34b-5p",
    "hsa-miR-143-5p",
    "hsa-miR-18a-5p",
    "hsa-miR-18b-5p",
    "hsa-miR-873-5p",
    "hsa-miR-10a-5p",
    "hsa-miR-10b-5p",
    "hsa-miR-196a-5p",
    "hsa-miR-193b-3p",
    "hsa-miR-3937",
    "hsa-miR-450b-3p",
]


def load_complete(parquet_path):
    """Inverse of the plate-loading step's save_complete.

    Returns dict[metric] -> DataFrame (condition index, time-point columns
    as int). The parquet schema is:
        index columns: [metric, condition]
        data columns:  time points stored as str ("0", "6", ...)
    """
    raw = pd.read_parquet(parquet_path)
    out = {}
    for metric, sub in raw.groupby("metric", sort=False):
        df = sub.drop(columns=["metric", "condition"]).copy()
        df.index = sub["condition"].values
        df.columns = [int(c) for c in df.columns]
        out[metric] = df
    return out


def _combos_singlets(s):
    """Split a Series-of-conditions index into combos (contain '+') and
    the de-duplicated single-miRNA components those combos reference."""
    combos = s.index[s.index.str.contains(r"\+")]
    singlets = []
    for combo in combos:
        a, b = combo.split(" + ")
        if a not in singlets:
            singlets.append(a)
        if b not in singlets:
            singlets.append(b)
    return combos, singlets


def abs_heatmap_df(complete, measure, time_slice, ordered=ORDERED_MIRS):
    """44x44 of mean(measure) per single agent (diagonal) and per
    combination (off-diagonal). Symmetric. Float32 to match v14.
    """
    s = complete[measure].loc[:, time_slice]
    combos, singlets = _combos_singlets(s)
    df = pd.DataFrame(index=ordered, columns=ordered)
    for single in singlets:
        df.loc[single, single] = s.loc[single].mean().mean()
    for comb in combos:
        one, two = comb.split(" + ")
        v = s.loc[comb].mean().mean()
        df.loc[one, two] = v
        df.loc[two, one] = v
    return df.astype(np.float32)


def hsa_df(complete, measure, time_slice, ordered=ORDERED_MIRS):
    """44x44 of combo-mean minus the better single-agent mean. Diagonal
    is 0. Direction depends on metric: cell-body metrics use min (we want
    fewer cell bodies than the best single agent); everything else uses
    max (we want longer neurites than the best single agent).
    """
    s = complete[measure].loc[:, time_slice]
    combos, singlets = _combos_singlets(s)
    df = pd.DataFrame(index=ordered, columns=ordered)
    for single in singlets:
        df.loc[single, single] = 0
    use_min = measure in ("cell body cluster area", "cell body clusters")
    agg = min if use_min else max
    for comb in combos:
        one, two = comb.split(" + ")
        v = s.loc[comb].mean().mean() - agg(
            s.loc[one].mean().mean(), s.loc[two].mean().mean()
        )
        df.loc[one, two] = v
        df.loc[two, one] = v
    return df.astype(np.float32)
