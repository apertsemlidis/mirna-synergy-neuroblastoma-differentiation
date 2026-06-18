# miRNA Synergy in Neuroblastoma Differentiation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Code and processed data for:

> Lawson SA, Zhang Y, Kosti A, Hart MJ, Penalva LOF, Pertsemlidis A. "MicroRNA combinations function as synergistic network regulators of neuroblastoma differentiation." *Journal of Biomedical Science* (2026).
>
> Preprint: [bioRxiv](https://doi.org/10.64898/2026.03.04.709166)

## What this repository ships

This is a **reproducibility archive**, not a turnkey figure-rendering pipeline. For each analysis domain it provides:

- the **analysis scripts** and their **pre-computed output CSVs**,
- the **processed input data** (the raw instrument data is too large / licensed separately — see *External Data*),
- the **final manuscript figures** as static `figures/figure_N.{pdf,png}` and `figures/additional_file_N.*`.

The figure-assembly composites (the scripts that lay out panels into the publication figures) are **not** included; they live in the private project tree. To reproduce a figure's underlying numbers, re-run the listed analysis script — it regenerates the CSVs under its `outputs/` (or domain) folder, which match the shipped copies. The static figure files are the exact images in the paper.

## Repository structure

```
├── 01_screen/                         # Primary combinatorial screen (Figures 3, 5, 6; Additional file 4)
│   ├── 01_load_plates.py              # Load raw Incucyte plates → normalized table (needs raw data)
│   ├── 02_compute_heatmaps.py         # 44×44 HSA + absolute heatmap tables (Figure 3)
│   ├── 03_compute_superhits.py        # Superhit calls (NL + CBCA significance gates)
│   ├── volcano_analysis.py            # Dual-phenotype NL volcano + stats (Figure 5A)
│   ├── qqplots.py                     # NL/CBCA p-value QQ plots (supplementary)
│   ├── pvalue_histograms.py           # NL/CBCA p-value histograms (supplementary; companion to qqplots)
│   ├── screen_helpers.py             # Shared screen functions
│   ├── nl_hsa_scores.csv              # Neurite-length HSA scores, 946 combinations
│   ├── cbca_scores.csv                # Cell-body cluster-area scores, 946 combinations
│   ├── nl_volcano_stats.csv           # Per-combination volcano statistics (Figure 5A)
│   ├── mirna_family_info.csv          # miRNA family annotations (same-family pair filtering)
│   ├── superhits.csv                  # Superhits passing NL + CBCA gates
│   ├── candidate_disposition_all_946.csv  # Full 946-combination disposition (Additional file 4 data; 34 hits typeset)
│   └── heatmap_data/                  # 96–126 h time-window matrices (Figure 3)
│
├── 02_dose_response/                  # Dose-response validation (Figure 7: HSA synergy surfaces)
│   ├── process_dose_response.py       # Raw plate exports → SynergyFinder-format CSVs (needs raw data)
│   ├── dose_response_hsa.py           # HSA synergy + replicate bootstrap (Figure 7 statistics)
│   ├── dose_response_maxnlnorm.csv    # Per-plate-max-normalized dose matrices (HSA input)
│   ├── dose_response_hsa_stats.csv    # Per-pair HSA synergy, 95% CI, p (2/6 pairs significant)
│   └── output/                        # Processed dose-response matrices, combined CSVs, synergy reports
│
├── 03_target_analysis/                # Target-space complementarity (Figures 4, 9, 10)
│   ├── target_complementarity.py      # Per-pair module coverage & complementarity (Figure 4)
│   ├── batch_analysis.py              # miRTarBase (validated-target) batch over the pair set
│   ├── statistical_tests.py           # Reinforcement-vs-dilution statistics
│   ├── synergy_features.py            # Synergistic vs non-synergistic feature comparison (Figure 9)
│   ├── nb_specific_analysis.py        # ADRN/MES/MYCN/retinoid module coverage (Figure 10)
│   ├── compare_databases.py           # TargetScan vs miRTarBase comparison
│   ├── ontology_venn.py               # Pathway-overlap Venn diagrams (supplementary)
│   ├── synergistic_pairs.csv          # 31 synergistic miRNA pairs (master input)
│   ├── external/nb_signatures/        # Curated ADRN/MES/MYCN/retinoid gene sets (shipped)
│   └── outputs/                       # Pre-computed metrics consumed by the figures
│       ├── pairs_summary.csv          # Per-pair incremental coverage summary
│       ├── per_pair/                  # Per-pair module metrics + term lists (Figure 4; 20 pairs)
│       ├── all_pairs_nb_metrics.csv   # NB-specific incremental coverage (Figure 10)
│       ├── nb_specific_stats.csv      # NB-module Mann-Whitney statistics (Figure 10)
│       ├── all_features.csv           # Synergy-feature comparison data (Figure 9)
│       └── synergy_features_stats.csv # Synergy-feature statistics (Figure 9)
│
├── 04_survival/                       # Patient survival (Figure 8, Additional files 1–3)
│   ├── km_3group.py                   # 3-group KM curves, per pair (Figure 8)
│   ├── km_3group_screen.py            # Sweep all 31 pairs
│   ├── km_mycn_stratified.py          # MYCN-amplified / non-amplified KM (Additional file 3)
│   ├── km_mycn_stratified_screen.py   # Sweep all 31 pairs, MYCN-stratified
│   ├── cox_forest_per_pair.py         # Per-pair Cox forest plots
│   ├── cox_forest_combined.py         # Combined Cox forest, six dose-response pairs (Additional file 2)
│   ├── cox_multivariate.py            # Cox PH adjusting for MYCN + age stratification
│   ├── cox_time_split.py              # Time-split HR sensitivity table (Additional file 1)
│   ├── cox_diagnostic.py              # Single-pair Cox diagnostic (miR-137-3p + miR-450b-5p)
│   ├── *_stats.csv                    # Pre-computed KM / Cox statistics
│   └── data/                          # GSE155945-derived patient data (expression + survival)
│
├── figures/                           # Final manuscript figures (static)
│   ├── figure_1.{pdf,png} … figure_10.{pdf,png}
│   ├── additional_file_1.csv          # Time-split Cox HR table
│   ├── additional_file_2.{pdf,png}, additional_file_3.{pdf,png}
│   └── additional_file_4.csv          # Candidate-disposition table (34 dual-positive hits)
│
├── LICENSE
└── README.md
```

## Figure → analysis script → data mapping

Each figure ships as a static `figures/figure_N.{pdf,png}`. The table below lists the analysis script and data that produce its underlying numbers (composites that lay out the panels are not shipped — see *What this repository ships*).

| Manuscript figure | Underlying analysis script | Data |
|---|---|---|
| Figure 1 (concept schematic) | — (illustration, no code) | — |
| Figure 2 (screening method) | — (microscopy / plate map) | raw Incucyte (not shipped) |
| Figure 3 (NL + CBCA heatmaps) | `01_screen/02_compute_heatmaps.py`, `03_compute_superhits.py` | `01_screen/heatmap_data/`, `superhits.csv` |
| Figure 4 (target complementarity) | `03_target_analysis/target_complementarity.py` + `batch_analysis.py` | `03_target_analysis/outputs/per_pair/` |
| Figure 5 (volcano + NL/CBCA correlation) | `01_screen/volcano_analysis.py` | `01_screen/{nl_hsa_scores,cbca_scores,nl_volcano_stats,mirna_family_info}.csv` |
| Figure 6 (screen NL/CBCA time course) | screen time-course composite (not shipped) | `01_screen/` screen data |
| Figure 7 (dose-response HSA synergy surfaces) | `02_dose_response/dose_response_hsa.py` (HSA + bootstrap) | `02_dose_response/dose_response_maxnlnorm.csv`, `output/` |
| Figure 8 (3-group KM survival) | `04_survival/km_3group.py` | `04_survival/data/{miRNA_expression_data,survival_data}.csv` |
| Figure 9 (synergy features) | `03_target_analysis/synergy_features.py` | `03_target_analysis/outputs/all_features.csv` |
| Figure 10 (NB-specific module coverage) | `03_target_analysis/nb_specific_analysis.py` | `03_target_analysis/outputs/all_pairs_nb_metrics.csv` |
| Additional file 1 (time-split Cox table) | `04_survival/cox_time_split.py` | `04_survival/data/` |
| Additional file 2 (combined Cox forest) | `04_survival/cox_forest_combined.py` | `04_survival/data/` |
| Additional file 3 (MYCN-stratified KM) | `04_survival/km_mycn_stratified.py` | `04_survival/data/` |
| Additional file 4 (candidate-disposition table) | `01_screen/build_additional_file_4.py` | `01_screen/candidate_disposition_all_946.csv` (full 946; 34 dual-positive hits typeset) |

Figures 1 and 2 are a vector schematic and a microscopy/plate-map panel assembled outside this codebase.

## What each subdirectory does

### 01_screen — combinatorial miRNA screen

946 pairwise miRNA combinations were screened in SK-N-BE(2)-C neuroblastoma cells for synergistic effects on neurite length (NL) and growth arrest (cell-body cluster area, CBCA). HSA (Highest Single Agent) scores identify combinations exceeding the strongest single agent; dual-phenotype filtering requires both NL synergy (p < 0.05, CI < 1) and CBCA improvement over ATRA (p < 0.05).

- `01_load_plates.py` → `02_compute_heatmaps.py` → `03_compute_superhits.py` is the screen-processing pipeline (split from the original monolith). Step 1 needs the raw Incucyte NeuroTrack data (≈152 MB; available on request); steps 2–3 produce the shipped `heatmap_data/` and `superhits.csv`.
- `volcano_analysis.py` produces the Figure 5A dual-phenotype volcano and emits `nl_volcano_stats.csv` (Fisher's combined test + BH correction).
- `candidate_disposition_all_946.csv` is the full disposition of all 946 screened combinations (screen metrics, dual-positive status, dose-response selection + outcome); the 34 dual-positive hits form the typeset **Additional file 4**. The screen-based NL/CBCA time courses for the six dose-response pairs underlie **Figure 6**.

### 02_dose_response — dose-response validation

Dose-response interaction modeling (SynergyFinder 3.0) for the six selected pairs across a 5×5 dose matrix. `process_dose_response.py` converts raw IncuCyte exports into the SynergyFinder-format CSVs in `output/`; the per-cell HSA synergy surfaces derived from these matrices are **Figure 7**. `dose_response_hsa.py` computes the interior-mean HSA synergy and a replicate bootstrap 95% CI / p-value per pair from `dose_response_maxnlnorm.csv`, writing `dose_response_hsa_stats.csv` — only the two miR-124-3p pairs reach significant HSA synergy (124+363: 17.3 [13.1–21.2]; 124+34b: 7.6 [3.4–11.9]; both p < 0.001). HSA is used because it is scale-free and reproducible across implementations, whereas the Bliss/ZIP scalars are tool- and normalization-dependent. Raw plate data is not included.

### 03_target_analysis — target-space complementarity

Incremental pathway coverage and target complementarity from TargetScan v7.2 predictions, across three analyses:

1. **Generic-module coverage (Figure 4)** — `target_complementarity.py`: on-target (neurite outgrowth) vs liability (apoptosis / ER stress) vs housekeeping pathways, computed per pair. The shipped `outputs/per_pair/` metrics are the pre-computed batch over all pairs.
2. **Synergy-feature comparison (Figure 9)** — `synergy_features.py`: Jaccard overlap, combined target-set size, individual potency, and tumor expression correlation, synergistic vs non-synergistic pairs.
3. **NB-specific modules (Figure 10)** — `nb_specific_analysis.py`: coverage of adrenergic (ADRN) / mesenchymal (MES) / MYCN-target / retinoid-response signatures (van Groningen et al., Wei et al. `WEI_MYCN_TARGETS_WITH_E_BOX`, GO:BP `GOBP_RESPONSE_TO_RETINOIC_ACID`).

`batch_analysis.py` and `compare_databases.py` provide the parallel miRTarBase (experimentally-validated target) comparison.

**miR-34b-5p handling (matches the published figures).** miR-34b-5p is absent from the conserved-default TargetScan v7.2 predictions; because it shares its seed family (AGGCAGU) with miR-449b-5p, the shipped target-space outputs assign it that partner's conserved target set (see the manuscript Methods, "Systematic feature comparison"). miR-450b-5p and miR-2110, which also lack conserved-default predictions, carry no target set. To regenerate the outputs from a fresh TargetScan download, apply the same seed-family assignment before running the scripts below.

```bash
# Per-pair target-space coverage (Figure 4) — single-pair example
python 03_target_analysis/target_complementarity.py \
  --mode predicted \
  --targetscan_tsv 03_target_analysis/external/targetscan72_hsa.tsv \
  --mirA hsa-miR-124-3p --mirB hsa-miR-363-3p \
  --outdir example_out

# NB-specific module coverage (Figure 10)
python 03_target_analysis/nb_specific_analysis.py \
  --pairs-csv 03_target_analysis/synergistic_pairs.csv \
  --targetscan-tsv 03_target_analysis/external/targetscan72_hsa.tsv \
  --sig-dir 03_target_analysis/external/nb_signatures \
  --outdir 03_target_analysis/outputs
```

### 04_survival — patient survival

Kaplan-Meier and Cox proportional-hazards analyses of 96 neuroblastoma patients (GSE155945; Misiak et al., 2021), stratified by coordinated miRNA expression. For each pair, patients are grouped by the number of the two miRNAs expressed above the cohort median (0, 1, or 2 "high"). Cox models adjust for MYCN amplification with age stratification; the time-split table tests the proportional-hazards assumption per pair. Survival figures cover the **six** dose-response pairs. Penalized (Firth-like) Cox estimation is used where the "both high" group has zero events (e.g., miR-137-3p + miR-450b-5p).

## External data (not included)

| Resource | Source | Place at |
|---|---|---|
| Raw Incucyte NeuroTrack imaging (~152 MB) | available on request | needed only by `01_screen/01_load_plates.py` |
| Raw IncuCyte dose-response plates | available on request | needed only by `02_dose_response/process_dose_response.py` |
| GSE155945 raw miRNA + clinical data | [GSE155945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155945) | derived CSVs already ship in `04_survival/data/` |
| TargetScan v7.2 (human predictions, 4.1 MB) | [targetscan.org](https://www.targetscan.org/) | `03_target_analysis/external/targetscan72_hsa.tsv` |
| miRTarBase MTI table | [mirtarbase.cuhk.edu.cn](https://mirtarbase.cuhk.edu.cn/) | `03_target_analysis/external/` (needed for `compare_databases.py` only) |

TargetScan and miRTarBase are not redistributed here — fetch them from their canonical sources to keep licensing and provenance clean. The `03_target_analysis/external/nb_signatures/` folder ships the project-curated gene-set inputs derived from the cited papers (van Groningen et al., Wei et al., GO:BP).

## Requirements

```
python >= 3.9
pandas >= 2.0
numpy >= 1.26
scipy >= 1.11
matplotlib >= 3.9
seaborn >= 0.13
lifelines >= 0.27      # 04_survival/ only
adjustText >= 0.8      # 01_screen/volcano_analysis.py only
requests >= 2.28       # ontology lookups in 03_target_analysis/ontology_venn.py
```

## Provenance

This repository is a curated export of the analysis scripts, processed data, and static figures behind the manuscript. The shipped figures correspond to the manuscript's final figure set (Figures 1–10, Additional files 1–4). Some processed-data CSVs are renamed from their analysis-internal forms for readability (e.g., `NL_allcombinations_vs_HSA.csv` → `nl_hsa_scores.csv`); the shipped scripts' `read_csv` paths point to these public names.

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
