# miRNA Synergy in Neuroblastoma Differentiation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Code and processed data for:

> Lawson SA, Zhang Y, Kosti A, Hart MJ, Penalva LOF, Pertsemlidis A. "MicroRNA combinations function as synergistic network regulators of neuroblastoma differentiation." *Journal of Biomedical Science* (2026).
>
> Preprint: [bioRxiv](https://doi.org/10.64898/2026.03.04.709166)

## Repository Structure

```
├── 01_screen/                              # Primary combinatorial screen (Figures 3, 5)
│   ├── screen_analysis.py                  # HSA synergy scoring from raw Incucyte data
│   ├── screen_helpers.py                   # Shared functions for screen_analysis.py
│   ├── volcano_plot.py                     # Dual-phenotype volcano (standalone Panel A of Figure 5)
│   ├── nl_cbca_correlation.py              # NL synergy vs CBCA improvement (standalone Panel B of Figure 5)
│   ├── qqplots.py                          # NL/CBCA p-value QQ plots (supplementary)
│   ├── nl_hsa_scores.csv                   # Neurite length HSA scores, 946 combinations
│   ├── cbca_scores.csv                     # Cell body cluster area scores, 946 combinations
│   ├── mirna_family_info.csv               # miRNA family annotations (same-family pair filtering)
│   ├── superhits.csv                       # 33 superhits passing NL + CBCA significance gates
│   └── heatmap_data/                       # 96–126 h time-window matrices (Figure 3)
│       ├── nl_hsa.csv                      # NL relative to highest single agent
│       ├── nl_absolute.csv                 # NL in absolute units
│       ├── cbca_hsa.csv                    # CBCA relative to lowest single agent
│       └── cbca_absolute.csv               # CBCA in absolute units
│
├── 02_dose_response/                       # Dose-response validation (Methods/Discussion only; no figure)
│   ├── process_dose_response.py            # Processing pipeline (requires raw plate data, not shipped)
│   └── output/                             # Processed dose-response artifacts
│
├── 03_target_analysis/                     # Target-space complementarity (Figures 4, 7, 8)
│   ├── target_complementarity.py           # Per-pair coverage & complementarity analysis (Figure 4)
│   ├── batch_analysis.py                   # Run target_complementarity across all 31 pairs
│   ├── statistical_tests.py                # NBS statistics + manuscript text generation
│   ├── nb_specific_analysis.py             # ADRN/MES/MYCN/retinoid module coverage (Figure 8)
│   ├── synergy_features.py                 # Distinguishing features of synergistic pairs (Figure 7)
│   ├── compare_databases.py                # TargetScan vs miRTarBase comparison
│   ├── ontology_venn.py                    # Pathway-overlap Venn diagrams (supplementary)
│   ├── synergistic_pairs.csv               # 31 synergistic miRNA pairs (input)
│   ├── external/                           # External resources (TargetScan, miRTarBase, NB signatures)
│   │   ├── targetscan72_hsa.tsv
│   │   ├── mirtarbase_mti.csv
│   │   └── nb_signatures/
│   └── outputs/                            # Pre-computed metrics for the figure composites
│       ├── pairs_summary.csv               # Per-pair Inc_on, Inc_ctrl, NBS
│       ├── per_pair/                       # Per-pair module metrics + term lists (Figure 4)
│       ├── all_pairs_nb_metrics.csv        # NB-specific incremental coverage (Figure 8)
│       └── all_features.csv                # Synergy-feature comparison data (Figure 7)
│
├── 04_survival/                            # Patient survival (Figure 6, Additional files 2, 3)
│   ├── km_3group.py                        # 3-group KM curves (Figure 6, per pair)
│   ├── km_3group_screen.py                 # Sweep all 31 pairs
│   ├── km_mycn_stratified.py               # MYCN-amplified/non-amplified KM (Additional file 3)
│   ├── km_mycn_stratified_screen.py        # Sweep all 31 pairs, MYCN-stratified
│   ├── cox_forest_per_pair.py              # Per-pair Cox forest plots
│   ├── cox_forest_combined.py              # Combined Cox forest, all four pairs (Additional file 2)
│   ├── cox_multivariate.py                 # Cox PH adjusting for MYCN + age stratification
│   ├── cox_time_split.py                   # Time-split HR sensitivity table (Additional file 1)
│   ├── cox_diagnostic.py                   # Single-pair Cox diagnostic (137 + 450b)
│   └── data/                               # GSE155945-derived patient data
│       ├── miRNA_expression_data.csv
│       └── survival_data.csv
│
├── figures/                                # Composite-assembly scripts for manuscript figures
│   ├── figure3_screen_heatmaps.py          # Figure 3: NL + CBCA heatmaps with hit overlays
│   ├── figure4_target_complementarity.py   # Figure 4 wrapper (uses figure4_panels.load_batch_data)
│   ├── figure4_panels.py                   # Figure 4 plotting logic (panel A boxplot + B/C scatter)
│   ├── figure5_volcano_correlation.py      # Figure 5: dual-phenotype volcano + NL/CBCA correlation
│   ├── figure6_km_3group.py                # Figure 6: 3-group KMs for the four synergistic pairs
│   ├── figure7_synergy_features.py         # Figure 7: features distinguishing synergistic pairs
│   ├── figure8_nb_specific.py              # Figure 8: NB-specific module incremental coverage
│   ├── additional_file_2_cox_forest.py     # Additional file 2 wrapper (delegates to cox_forest_combined.py)
│   ├── additional_file_3_km_mycn.py        # Additional file 3: MYCN-stratified KMs
│   └── figure_style.py                     # Shared matplotlib style settings
│
├── POPULATE.sh                             # Source-of-truth script for repopulating from the project tree
├── LICENSE
└── README.md
```

## Figure → Script → Data Mapping

| Manuscript artifact | Composite script | Data dependency |
|---|---|---|
| Figure 3 | `figures/figure3_screen_heatmaps.py` | `01_screen/heatmap_data/`, `01_screen/superhits.csv` |
| Figure 4 | `figures/figure4_target_complementarity.py` | `03_target_analysis/outputs/per_pair/` |
| Figure 5 | `figures/figure5_volcano_correlation.py` | `01_screen/{nl_hsa_scores,cbca_scores,mirna_family_info}.csv` |
| Figure 6 | `figures/figure6_km_3group.py` | `04_survival/data/{miRNA_expression_data,survival_data}.csv` |
| Figure 7 | `figures/figure7_synergy_features.py` | `03_target_analysis/outputs/all_features.csv` |
| Figure 8 | `figures/figure8_nb_specific.py` | `03_target_analysis/outputs/all_pairs_nb_metrics.csv` |
| Additional file 1 | `04_survival/cox_time_split.py` (pre-rendered copy at `figures/Additional file 1 v15.csv`) | `04_survival/data/` |
| Additional file 2 | `figures/additional_file_2_cox_forest.py` → `04_survival/cox_forest_combined.py` | `04_survival/data/` |
| Additional file 3 | `figures/additional_file_3_km_mycn.py` | `04_survival/data/` |

Each composite is self-contained: it loads the shipped CSVs, runs the analysis, and saves both PNG and PDF outputs into `figures/` under the canonical `Figure N v15.{png,pdf}` / `Additional file N v15.{png,pdf}` name. Running any composite from the repo root works without arguments:

```bash
python figures/figure3_screen_heatmaps.py
python figures/figure5_volcano_correlation.py
python figures/figure6_km_3group.py
# ... etc.
```

Figures 1 and 2 are vector schematics / microscopy panels assembled outside this codebase; they are not regenerated by any script in this repo.

## What Each Subdirectory Does

### 01_screen — Combinatorial miRNA Screen

946 pairwise miRNA combinations were screened in SK-N-BE(2)-C neuroblastoma cells for synergistic effects on neurite length and growth arrest. HSA (Highest Single Agent) synergy scores identify combinations exceeding the strongest single agent. Dual-phenotype filtering requires both NL synergy (p < 0.05, CI < 1) and CBCA improvement over ATRA (p < 0.05).

- `screen_analysis.py` is the original notebook-derived pipeline that processes raw Incucyte NeuroTrack data (not in this repo — 152 MB; available on request) into the summary CSVs and `heatmap_data/`. Included here for reproducibility and audit. Self-contained scripts below operate on the shipped CSV outputs.
- `volcano_plot.py`, `nl_cbca_correlation.py` produce the original standalone Figure 5 panels and also emit the manuscript's summary statistics (Fisher's combined test + BH, permutation tests).

### 02_dose_response — Synergy Validation

Dose-response interaction modeling using SynergyFinder 3.0. Validation results for the four synergistic pairs are discussed in the manuscript prose (no main-text figure). Raw IncuCyte plate data is not included; the `output/` folder ships the processed dose-response curves and synergy reports.

### 03_target_analysis — Target-Space Complementarity

Computes incremental pathway coverage and target complementarity using TargetScan v7.2 predictions. Three analyses:

1. **Generic-module coverage (Figure 4)** — `target_complementarity.py` + `batch_analysis.py`. On-target (neurite outgrowth) vs liability (e.g., apoptosis) vs housekeeping pathways.
2. **NB-specific modules (Figure 8)** — `nb_specific_analysis.py`. Coverage of adrenergic/mesenchymal/MYCN/retinoid signatures from van Groningen et al., Wei et al., and GO:BP retinoid response.
3. **Synergy-feature comparison (Figure 7)** — `synergy_features.py`. Jaccard, combined target set size, individual potency, and tumor expression correlation, comparing synergistic vs non-synergistic pairs.

```bash
# Single pair
python 03_target_analysis/target_complementarity.py \
  --mode predicted \
  --targetscan_tsv 03_target_analysis/external/targetscan72_hsa.tsv \
  --mirA hsa-miR-124-3p --mirB hsa-miR-363-3p \
  --outdir output_example

# All 31 pairs
python 03_target_analysis/batch_analysis.py \
  --pairs_csv 03_target_analysis/synergistic_pairs.csv \
  --targetscan_tsv 03_target_analysis/external/targetscan72_hsa.tsv \
  --outdir 03_target_analysis/outputs

# Statistics
python 03_target_analysis/statistical_tests.py 03_target_analysis/outputs
```

### 04_survival — Patient Survival

Kaplan-Meier and Cox proportional hazards analyses stratifying 96 neuroblastoma patients by coordinated miRNA expression. Source: GSE155945 (Misiak et al., 2021). For each synergistic pair, patients are tertile-split into "both high", "mixed", and "both low" groups. Cox models adjust for MYCN amplification with age stratification; the time-split sensitivity table tests the proportional hazards assumption per pair.

Penalized (Firth-like) Cox estimation is used for `miR-137-3p + miR-450b-5p`, where the "both high" group has zero events.

## External Data (not included)

| Resource | Source | Place at |
|---|---|---|
| Raw Incucyte NeuroTrack imaging (~152 MB) | Available upon request | needed only by `01_screen/screen_analysis.py` |
| Raw IncuCyte dose-response plates | Available upon request | needed only by `02_dose_response/process_dose_response.py` |
| GSE155945 raw miRNA + clinical data | [GSE155945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155945) | derived CSVs already ship in `04_survival/data/` |
| TargetScan v7.2 (human predictions, 4.1 MB) | [targetscan.org](https://www.targetscan.org/) | `03_target_analysis/external/targetscan72_hsa.tsv` |
| miRTarBase MTI table (375 MB) | [mirtarbase.cuhk.edu.cn](https://mirtarbase.cuhk.edu.cn/) | `03_target_analysis/external/mirtarbase_mti.csv` (needed for `compare_databases.py` only) |

External databases (TargetScan, miRTarBase) are not redistributed here — fetch them from their canonical sources to keep licensing and provenance clean. The `03_target_analysis/external/nb_signatures/` folder, by contrast, ships the project-curated gene-set inputs derived from the cited papers (van Groningen et al., Wei et al., GO:BP).

## Requirements

```
python >= 3.9
pandas >= 2.0
numpy >= 1.26
scipy >= 1.11
matplotlib >= 3.9
seaborn >= 0.13
lifelines >= 0.27      # 04_survival/ only
adjustText >= 0.8      # volcano_plot.py only
requests >= 2.28       # ontology lookups in 03_target_analysis/
```

## Provenance

This repository is populated from the project tree by `POPULATE.sh`. Script filenames in the project carry `_v14` / `_v15` version suffixes (the project's never-overwrite versioning scheme); those suffixes are dropped here because the Zenodo-archived release is itself a fixed version. CSV filenames are likewise renamed from their project-internal forms (e.g., `NL_allcombinations_vs_HSA.csv` → `nl_hsa_scores.csv`); the populator script patches the shipped `read_csv` calls to match.

## License

MIT — see [LICENSE](LICENSE).

## Citation

```bibtex
@article{lawson2026mirna,
  title={MicroRNA combinations function as synergistic network regulators
         of neuroblastoma differentiation},
  author={Lawson, Seth A and Zhang, Yiqiang and Kosti, Adam
          and Hart, Matthew J and Penalva, Luiz O F
          and Pertsemlidis, Alexander},
  journal={Journal of Biomedical Science},
  year={2026}
}
```
