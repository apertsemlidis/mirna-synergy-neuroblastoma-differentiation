#!/usr/bin/env python
"""Assemble the candidate-disposition supplementary table (Additional file 4).

Derives the 34 dual-positive screen hits from the full 946-combination master
(analysis/decision_artifacts/candidate_disposition_table_946.csv) and adds the dose-response
validation outcome, reconciled with the v19 synergy analysis (analysis/synergyfinder_version_repro/):
of the six pairs taken to dose-response, only the two miR-124-3p pairs validate as HSA-synergistic
(replicate bootstrap, p<0.001); the others are additive-to-antagonistic. HSA is the validation metric
because it is scale-free and implementation-robust (the Bliss/ZIP scalar is tool- and
normalization-dependent — see analysis/synergyfinder_version_repro/README.md).

Rows are ordered by combination index ascending (strongest synergy first) across all 34 — neutral and
reproducible; disposition is carried by the selected/outcome columns, so the order transparently shows
that several non-pursued hits have stronger CI than the six pursued pairs.

Outputs (figures/, for the manuscript supplement): Additional file 4 v19.csv + .xlsx
The full 946-combination master also ships as a repo/Zenodo data file
(github_repo/01_screen/candidate_disposition_all_946.csv), not as a typeset supplement.
"""

import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# Source of truth = the full 946-combination disposition table (tracked). The 34 dual-positive
# hits are derived from it here, so the supplement and the repo/Zenodo data file stay consistent.
SRC = ROOT / "analysis" / "decision_artifacts" / "candidate_disposition_table_946.csv"

# dose-response outcome, reconciled with the v19 HSA bootstrap (run_hsa_bootstrap.py)
DR_OUTCOME = {
    "hsa-miR-124-3p + hsa-miR-363-3p": "Validated — HSA synergy 17.3 (95% CI 13.1–21.2), p<0.001",
    "hsa-miR-124-3p + hsa-miR-34b-5p": "Validated — HSA synergy 7.6 (95% CI 3.4–11.9), p<0.001",
    "hsa-miR-137 + hsa-miR-450b-3p": "Not validated — additive (HSA NS)",
    "hsa-miR-137 + hsa-miR-449b-5p": "Not validated — additive (HSA NS)",
    "hsa-miR-137 + hsa-miR-17-5p": "Not validated — antagonistic trend (HSA NS)",
    "hsa-miR-19b-3p + hsa-miR-2110": "Not validated — antagonistic (HSA NS)",
}

df = pd.read_csv(SRC)
# normalize boolean-ish columns (the master may carry TRUE/FALSE or True/False)
for c in ("dual_positive", "selected_for_dose_response", "same_family"):
    df[c] = df[c].astype(str).str.strip().str.upper() == "TRUE"
df = df[df["dual_positive"]].copy()  # 34 dual-positive hits
df["dose_response_outcome"] = (
    df["combination"].map(DR_OUTCOME).fillna("Not tested (screen-level hit)")
)

# tidy presentation: round screen metrics, friendlier headers, sort hits by synergy strength
df = df.rename(
    columns={
        "combination": "miRNA pair",
        "nl_synergy": "NL synergy (combo − HSA)",
        "combination_index": "combination index (CI)",
        "nl_synergy_pvalue": "NL synergy p",
        "cbca_improvement": "CBCA improvement vs ATRA",
        "same_family": "same seed family",
        "dual_positive": "dual-positive",
        "selected_for_dose_response": "selected for dose-response",
        "dose_response_outcome": "dose-response outcome",
    }
)
for c in [
    "NL synergy (combo − HSA)",
    "combination index (CI)",
    "CBCA improvement vs ATRA",
]:
    df[c] = df[c].round(3)
df["NL synergy p"] = df["NL synergy p"].map(lambda v: f"{v:.2g}")
# Order by synergy strength (combination index ascending = strongest synergy first) across all 34.
# Neutral/reproducible; disposition is carried by the selected/outcome columns, not the row order —
# this transparently shows that several non-pursued hits have stronger CI than the six pursued pairs.
df = df.sort_values("combination index (CI)", ascending=True)

out_csv = ROOT / "figures" / "Additional file 4 v19.csv"
out_xlsx = ROOT / "figures" / "Additional file 4 v19.xlsx"
df.to_csv(out_csv, index=False)
try:
    df.to_excel(out_xlsx, index=False)
    print(f"wrote {out_csv.name} + {out_xlsx.name}  ({len(df)} hits)")
except Exception as e:
    print(f"wrote {out_csv.name}  ({len(df)} hits); xlsx skipped: {e}")
