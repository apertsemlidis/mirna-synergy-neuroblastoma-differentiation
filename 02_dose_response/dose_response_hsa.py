#!/usr/bin/env python
"""Dose-response HSA synergy with replicate bootstrap CIs — the statistics behind Figure 7.

For each miRNA pair in a SynergyFinder-format dose-response table (columns
PairIndex, Response, Drug1, Drug2, Conc1, Conc2, ConcUnit; replicate-level rows),
computes the interior-mean Highest Single Agent (HSA) synergy and a replicate
bootstrap 95% CI + one-sided p-value.

HSA synergy is scale-free (combination - best single agent), so combo>best-single
significance is invariant to the response normalization. Per cell:
    HSA delta = mean(combination) - max(mean(single_A), mean(single_B))
the figure-7 score is the mean over the 4x4 non-zero-dose ("interior") cells; the CI
is from resampling the replicates within each cell.

Input : dose_response_maxnlnorm.csv  (per-plate-maximum-normalized neurite length @120 h)
Output: dose_response_hsa_stats.csv  (pair, hsa_synergy, ci_low, ci_high, p_value, n_validated)

A pair is HSA-synergistic if the CI lower bound > 0 (one-sided p < 0.025). On the
manuscript input only the two miR-124-3p pairs pass (124+363, 124+34b; both p < 0.001).

Usage: python dose_response_hsa.py [input.csv] [output.csv]
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INP = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "dose_response_maxnlnorm.csv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "dose_response_hsa_stats.csv"
B = 5000
rng = np.random.default_rng(1)

df = pd.read_csv(INP)
doses = sorted(set(df["Conc1"]) | set(df["Conc2"]))
inter = [(i, j) for i in doses[1:] for j in doses[1:]]

rows = []
for pi, g in df.groupby("PairIndex"):
    pair = f"{g['Drug1'].iloc[0]} + {g['Drug2'].iloc[0]}"
    cells = {
        (i, j): g[(g.Conc1 == i) & (g.Conc2 == j)].Response.values
        for i in doses
        for j in doses
        if len(g[(g.Conc1 == i) & (g.Conc2 == j)])
    }

    def score(means):
        return float(
            np.mean(
                [
                    means[(i, j)] - max(means[(i, 0.0)], means[(0.0, j)])
                    for (i, j) in inter
                ]
            )
        )

    obs = score({k: v.mean() for k, v in cells.items()})
    bs = np.array(
        [
            score(
                {
                    k: rng.choice(v, size=len(v), replace=True).mean()
                    for k, v in cells.items()
                }
            )
            for _ in range(B)
        ]
    )
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p1 = float(np.mean(bs <= 0))
    rows.append(
        dict(
            pair=pair,
            hsa_synergy=round(obs, 2),
            ci_low=round(lo, 2),
            ci_high=round(hi, 2),
            p_value=p1,
            hsa_synergistic=bool(lo > 0),
        )
    )

out = pd.DataFrame(rows).sort_values("hsa_synergy", ascending=False)
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
print(
    f"\nwrote {OUT.name} — {int(out['hsa_synergistic'].sum())}/{len(out)} pairs HSA-synergistic (CI excludes 0)"
)
