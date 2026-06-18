#!/usr/bin/env python
"""Assemble the candidate-disposition supplementary table (Additional file 4).

Derives the 34 dual-positive screen hits from the full 946-combination disposition table
and adds the dose-response validation outcome: of the six pairs taken to dose-response, only
the two miR-124-3p pairs validate as HSA-synergistic (replicate bootstrap, p<0.001); the
others are additive-to-antagonistic. HSA is the validation metric because it is scale-free and
implementation-robust, whereas the Bliss/ZIP scalars are tool- and normalization-dependent.

Rows are ordered by combination index ascending (strongest synergy first) across all 34 —
neutral and reproducible; disposition is carried by the selected/outcome columns, so the order
transparently shows that several non-pursued hits have stronger CI than the six pursued pairs.

Usage:  python build_additional_file_4.py [input_946.csv] [output_af4.csv]
Defaults (repo layout): reads candidate_disposition_all_946.csv next to this script and writes
additional_file_4.csv alongside it.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
INP = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else HERE / "candidate_disposition_all_946.csv"
)
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "additional_file_4.csv"

# dose-response outcome, reconciled with the HSA bootstrap (dose_response_hsa.py)
DR_OUTCOME = {
    "hsa-miR-124-3p + hsa-miR-363-3p": "Validated — HSA synergy 17.3 (95% CI 13.1–21.2), p<0.001",
    "hsa-miR-124-3p + hsa-miR-34b-5p": "Validated — HSA synergy 7.6 (95% CI 3.4–11.9), p<0.001",
    "hsa-miR-137 + hsa-miR-450b-3p": "Not validated — additive (HSA NS)",
    "hsa-miR-137 + hsa-miR-449b-5p": "Not validated — additive (HSA NS)",
    "hsa-miR-137 + hsa-miR-17-5p": "Not validated — antagonistic trend (HSA NS)",
    "hsa-miR-19b-3p + hsa-miR-2110": "Not validated — antagonistic (HSA NS)",
}

df = pd.read_csv(INP)
# normalize boolean-ish columns (the master may carry TRUE/FALSE or True/False)
for c in ("dual_positive", "selected_for_dose_response", "same_family"):
    df[c] = df[c].astype(str).str.strip().str.upper() == "TRUE"
df = df[df["dual_positive"]].copy()  # 34 dual-positive hits
df["dose_response_outcome"] = (
    df["combination"].map(DR_OUTCOME).fillna("Not tested (screen-level hit)")
)

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
df = df.sort_values("combination index (CI)", ascending=True)

df.to_csv(OUT, index=False)
print(f"wrote {OUT.name}  ({len(df)} hits)")
try:
    df.to_excel(OUT.with_suffix(".xlsx"), index=False)
    print(f"wrote {OUT.with_suffix('.xlsx').name}")
except Exception as e:
    print(f"(xlsx skipped: {e})")
