#!/bin/bash
# Populate the GitHub repo with manuscript-relevant files only.
# Run from the project root:
#   bash github_repo/POPULATE.sh
#
# Path 4 model (v16):
#   Ships analysis scripts + their stats CSVs + static figure outputs.
#   Does NOT ship figure-rendering composites — those stay private in the
#   project under scripts/figures/. A reviewer reproducing the work can:
#     (a) read the shipped statistics CSVs alongside the figure PDFs, or
#     (b) re-run the analysis scripts to regenerate the statistics CSVs.
#
# Version conventions:
#   - Project carries versioned filenames (_v14.py, _v16.py).
#   - Public repo strips version suffixes — the GitHub release tag
#     (v1.0.0 etc.) carries the snapshot identity.
#
# Updated for v19 (2026-06-17): Figure 6 = screen time-course, Figure 7 = dose-response
# HSA synergy surfaces (new), old Figures 7/8/9 renumbered to 8/9/10; Additional file 4
# (candidate-disposition table) added + full 946-combination master shipped as data.
# Static figures (Fig 1-10 + Add files 1-4) shipped as PDF/PNG; figure composites excluded.

set -u

REPO="github_repo"
PAPER="$(pwd)"

if [ ! -d "$PAPER/$REPO" ]; then
    echo "ERROR: $REPO/ not found. Run from the project root."
    exit 1
fi

echo "Source: $PAPER"
echo "Target: $REPO/"
echo ""

# Wipe previously populated subdirs to avoid leaving stale files.
# Preserves: README.md, LICENSE, POPULATE.sh, .gitignore, .git, .DS_Store
for sub in 01_screen 02_dose_response 03_target_analysis 04_survival figures; do
    if [ -d "$PAPER/$REPO/$sub" ]; then
        rm -rf "$PAPER/$REPO/$sub"
    fi
done

# Helper: copy with explicit error message on missing source.
COPIES=0
MISSING=0
cpy() {
    local src="$1"
    local dst="$2"
    if [ -e "$src" ]; then
        mkdir -p "$(dirname "$dst")"
        cp -R "$src" "$dst"
        COPIES=$((COPIES + 1))
    else
        echo "  MISSING: $src"
        MISSING=$((MISSING + 1))
    fi
}

# ---------------------------------------------------------------------
# 01_screen — Primary combinatorial screen (Figures 3, 5)
# ---------------------------------------------------------------------
echo "[01_screen]"
# Analysis scripts — current three-step tier (split out from the archived
# v14 monolith) plus the volcano + qqplots producers. The NL/CBCA
# correlation panel script (generate_nl_cbca_correlation_panel.py) is kept
# in the project tree but not shipped publicly: Figure 5 v17's composite
# is self-contained, so the standalone panel script isn't needed for
# manuscript reproducibility.
cpy "$PAPER/scripts/screen/load_screen.py"                         "$REPO/01_screen/01_load_plates.py"
cpy "$PAPER/scripts/screen/compute_hsa_abs_tables.py"              "$REPO/01_screen/02_compute_heatmaps.py"
cpy "$PAPER/scripts/screen/compute_superhits.py"                   "$REPO/01_screen/03_compute_superhits.py"
cpy "$PAPER/scripts/screen/screen_helpers.py"                      "$REPO/01_screen/screen_helpers.py"
cpy "$PAPER/scripts/screen/generate_nl_volcano_dual.py"            "$REPO/01_screen/volcano_analysis.py"
cpy "$PAPER/scripts/screen/generate_qqplots.py"                    "$REPO/01_screen/qqplots.py"
cpy "$PAPER/scripts/screen/generate_pvalue_histograms.py"          "$REPO/01_screen/pvalue_histograms.py"
# Screen scores + family annotations (inputs)
cpy "$PAPER/data/screen/NL_allcombinations_vs_HSA.csv"             "$REPO/01_screen/nl_hsa_scores.csv"
cpy "$PAPER/data/screen/CBCA_allcombinations_vs_ATRA.csv"          "$REPO/01_screen/cbca_scores.csv"
cpy "$PAPER/data/screen/miR_Family_Info.csv"                       "$REPO/01_screen/mirna_family_info.csv"
# Heatmap inputs (Figure 3) — pre-processed CSVs; raw Incucyte data not included (152 MB)
cpy "$PAPER/data/screen/HSA_dfs/nl_slice(96, 126, None).csv"       "$REPO/01_screen/heatmap_data/nl_hsa.csv"
cpy "$PAPER/data/screen/HSA_dfs/cbca_slice(96, 126, None).csv"     "$REPO/01_screen/heatmap_data/cbca_hsa.csv"
cpy "$PAPER/data/screen/ABS_dfs/nl_slice(96, 126, None).csv"       "$REPO/01_screen/heatmap_data/nl_absolute.csv"
cpy "$PAPER/data/screen/ABS_dfs/cbca_slice(96, 126, None).csv"     "$REPO/01_screen/heatmap_data/cbca_absolute.csv"
cpy "$PAPER/data/screen/HSAhits_p05_cytostatic.csv"                "$REPO/01_screen/superhits.csv"
# Stats CSV from analysis (v16)
cpy "$PAPER/data/screen/nl_volcano_stats.csv"                      "$REPO/01_screen/nl_volcano_stats.csv"
# Full 946-combination candidate-disposition table (the Zenodo data file behind Additional file 4;
# the 34 dual-positive hits are the typeset supplement, derived from this master)
cpy "$PAPER/analysis/decision_artifacts/candidate_disposition_table_946.csv" "$REPO/01_screen/candidate_disposition_all_946.csv"

# ---------------------------------------------------------------------
# 02_dose_response — Dose-response validation
# ---------------------------------------------------------------------
echo "[02_dose_response]"
cpy "$PAPER/scripts/dose_response/process_dose_response.py"        "$REPO/02_dose_response/process_dose_response.py"
# Raw plate data omitted (large); processed outputs included
cpy "$PAPER/data/dose_response/output"                              "$REPO/02_dose_response/output"

# ---------------------------------------------------------------------
# 03_target_analysis — Target-space complementarity (Figures 4, 8, 9)
# ---------------------------------------------------------------------
echo "[03_target_analysis]"
# Analysis scripts (v16 where applicable; v14 for utility)
cpy "$PAPER/scripts/target_analysis/mirna_synergy_coverage_complementarity.py" "$REPO/03_target_analysis/target_complementarity.py"
cpy "$PAPER/scripts/target_analysis/mirna_synergy_batch_mirtarbase.py"          "$REPO/03_target_analysis/batch_analysis.py"
cpy "$PAPER/scripts/target_analysis/analyze_reinforce_dilute.py"                "$REPO/03_target_analysis/statistical_tests.py"
cpy "$PAPER/scripts/target_analysis/nb_specific_analysis.py"                        "$REPO/03_target_analysis/nb_specific_analysis.py"
cpy "$PAPER/scripts/target_analysis/synergy_features.py"                            "$REPO/03_target_analysis/synergy_features.py"
cpy "$PAPER/scripts/target_analysis/compare_databases.py"                           "$REPO/03_target_analysis/compare_databases.py"
cpy "$PAPER/scripts/target_analysis/mirna_synergy_ontology_venn.py"                 "$REPO/03_target_analysis/ontology_venn.py"
# Inputs
cpy "$PAPER/mirna_pairs.csv"                                                         "$REPO/03_target_analysis/synergistic_pairs.csv"
# External databases NOT redistributed (licensing/provenance; see .gitignore + README).
# (The published v18 figures use the "B-partial" TargetScan input, in which miR-34b-5p is
# assigned its seed-family partner miR-449b-5p's conserved target set; this is described in
# the manuscript Methods. The shipped outputs below already reflect it.)
cpy "$PAPER/data/external/nb_signatures"                                             "$REPO/03_target_analysis/external/nb_signatures"
# Pre-computed outputs (per-pair metrics + module term lists) — v18 B-partial, matching the paper
cpy "$PAPER/data/target_analysis/batch_ts_v18/pairs_summary.csv"                      "$REPO/03_target_analysis/outputs/pairs_summary.csv"
cpy "$PAPER/data/target_analysis/batch_ts_nb_specific_v18/all_pairs_nb_metrics.csv"   "$REPO/03_target_analysis/outputs/all_pairs_nb_metrics.csv"
cpy "$PAPER/data/target_analysis/synergy_features_v18/all_features.csv"               "$REPO/03_target_analysis/outputs/all_features.csv"
cpy "$PAPER/data/target_analysis/batch_ts_v18/per_pair"                               "$REPO/03_target_analysis/outputs/per_pair"
# Stats CSVs from analysis (v18 B-partial)
cpy "$PAPER/data/target_analysis/synergy_features_v18/synergy_features_stats.csv"     "$REPO/03_target_analysis/outputs/synergy_features_stats.csv"
cpy "$PAPER/data/target_analysis/batch_ts_nb_specific_v18/nb_specific_stats.csv"      "$REPO/03_target_analysis/outputs/nb_specific_stats.csv"

# ---------------------------------------------------------------------
# 04_survival — Patient survival (Figure 7, Additional files 1, 2, 3)
# Source: GSE155945 (Misiak et al. 2021), 96 NB tumors
# ---------------------------------------------------------------------
echo "[04_survival]"
# Analysis scripts (v16 master pair set + stats emission)
cpy "$PAPER/survival/km_3group/km_3group_all.py"                    "$REPO/04_survival/km_3group.py"
cpy "$PAPER/survival/km_3group/km_3group_screen.py"                 "$REPO/04_survival/km_3group_screen.py"
cpy "$PAPER/survival/km_mycn_stratified/km_mycn_stratified_all.py"       "$REPO/04_survival/km_mycn_stratified.py"
cpy "$PAPER/survival/km_mycn_stratified/km_mycn_stratified_screen.py"    "$REPO/04_survival/km_mycn_stratified_screen.py"
cpy "$PAPER/survival/cox_forest/cox_forest_all.py"                  "$REPO/04_survival/cox_forest_per_pair.py"
cpy "$PAPER/survival/cox_forest/cox_forest_combined.py"             "$REPO/04_survival/cox_forest_combined.py"
cpy "$PAPER/survival/cox_multivariate/multivariate_survival.py"     "$REPO/04_survival/cox_multivariate.py"
cpy "$PAPER/survival/cox_time_split/cox_time_split_table.py"        "$REPO/04_survival/cox_time_split.py"
cpy "$PAPER/survival/cox_diagnostic/cox_diagnostic_137+450b.py"     "$REPO/04_survival/cox_diagnostic.py"
# Data (GSE155945-derived)
cpy "$PAPER/survival/miRNA_expression_data.csv"                     "$REPO/04_survival/data/miRNA_expression_data.csv"
cpy "$PAPER/survival/survival_data.csv"                             "$REPO/04_survival/data/survival_data.csv"
# Stats CSVs from analysis (v16)
cpy "$PAPER/survival/km_3group/km_3group_stats.csv"                 "$REPO/04_survival/km_3group_stats.csv"
cpy "$PAPER/survival/km_mycn_stratified/km_mycn_stratified_stats.csv"   "$REPO/04_survival/km_mycn_stratified_stats.csv"
cpy "$PAPER/survival/cox_forest/cox_forest_combined_stats.csv"      "$REPO/04_survival/cox_forest_combined_stats.csv"
cpy "$PAPER/survival/cox_time_split/cox_time_split_table.csv"       "$REPO/04_survival/cox_time_split_table.csv"

# ---------------------------------------------------------------------
# figures — Static deliverables only (PDF / PNG); composites NOT shipped
# ---------------------------------------------------------------------
echo "[figures]"
# Figures 1-10. Prefer the latest available vintage per figure: v19 → v18 → v17 → v16 → v15.
# (v19 submission set: Figure 6 = screen time-course, Figure 7 = dose-response HSA surfaces (new),
# old Figures 7/8/9 renumbered to 8/9/10. All v19 figure files are labelled v19; older vintages in
# the fallback chain are kept for backward compatibility and live in figures/archive/.)
for n in 1 2 3 4 5 6 7 8 9 10; do
    for v in v19 v18 v17 v16 v15; do
        if [ -f "$PAPER/figures/Figure $n $v.pdf" ]; then
            cpy "$PAPER/figures/Figure $n $v.pdf"   "$REPO/figures/figure_$n.pdf"
            break
        fi
    done
    for v in v19 v18 v17 v16 v15; do
        if [ -f "$PAPER/figures/Figure $n $v.png" ]; then
            cpy "$PAPER/figures/Figure $n $v.png"   "$REPO/figures/figure_$n.png"
            break
        fi
    done
done
# Additional file 1 (survival Cox table, CSV)
for v in v19 v18 v16 v15; do
    if [ -f "$PAPER/figures/Additional file 1 $v.csv" ]; then
        cpy "$PAPER/figures/Additional file 1 $v.csv" "$REPO/figures/additional_file_1.csv"
        break
    fi
done
# Additional files 2-3 (figures)
for n in 2 3; do
    for v in v19 v18 v16 v15; do
        if [ -f "$PAPER/figures/Additional file $n $v.pdf" ]; then
            cpy "$PAPER/figures/Additional file $n $v.pdf" "$REPO/figures/additional_file_$n.pdf"
            break
        fi
    done
    for v in v19 v18 v16 v15; do
        if [ -f "$PAPER/figures/Additional file $n $v.png" ]; then
            cpy "$PAPER/figures/Additional file $n $v.png" "$REPO/figures/additional_file_$n.png"
            break
        fi
    done
done
# Additional file 4 (candidate-disposition table, CSV; new in v19)
for v in v19 v18; do
    if [ -f "$PAPER/figures/Additional file 4 $v.csv" ]; then
        cpy "$PAPER/figures/Additional file 4 $v.csv" "$REPO/figures/additional_file_4.csv"
        break
    fi
done

# ---------------------------------------------------------------------
# Patch filename + path references in the populated analysis scripts
#
# The public repo ships CSVs under cleaner public-facing names than the
# project uses internally. Analysis scripts copied above still reference
# the original names; patch them. Also fix path resolution where the
# project's nested layout (data/screen, survival/<subdir>) collapses into
# the flat 0N_*/ structure in the public repo.
# ---------------------------------------------------------------------
echo ""
echo "Patching filename + path references..."

# Use a portable in-place sed (BSD/macOS-friendly: -i ''; GNU: -i).
SED_INPLACE() {
    if sed --version >/dev/null 2>&1; then
        sed -i "$@"          # GNU
    else
        sed -i '' "$@"       # BSD/macOS
    fi
}

# Global literal renames applied across every .py file in the public repo.
RENAMES=(
    "s|NL_allcombinations_vs_HSA\.csv|nl_hsa_scores.csv|g"
    "s|CBCA_allcombinations_vs_ATRA\.csv|cbca_scores.csv|g"
    "s|miR_Family_Info\.csv|mirna_family_info.csv|g"
    "s|HSAhits_p05_cytostatic\.csv|superhits.csv|g"
    "s|miRTarBase_MTI\.csv|mirtarbase_mti.csv|g"
)
# Per-pair output filenames in survival scripts are now unsuffixed
# upstream (post-2026-05 convention sweep); no public-side rename needed.
while IFS= read -r -d '' pyfile; do
    for expr in "${RENAMES[@]}"; do
        SED_INPLACE -e "$expr" "$pyfile"
    done
done < <(find "$REPO" -name "*.py" -print0)

# 01_screen/: scripts that os.chdir into data/screen now chdir to their own folder.
for f in "$REPO/01_screen/volcano_analysis.py" \
         "$REPO/01_screen/qqplots.py" \
         "$REPO/01_screen/pvalue_histograms.py"; do
    [ -f "$f" ] && SED_INPLACE \
        -e 's|os\.chdir(Path(__file__)\.resolve()\.parents\[2\] / "data" / "screen")|os.chdir(Path(__file__).resolve().parent)|g' \
        "$f"
done

# 04_survival/: data CSVs live in 04_survival/data/, not survival/ parent.
for f in "$REPO/04_survival"/*.py; do
    [ -f "$f" ] && SED_INPLACE \
        -e 's|"\.\./miRNA_expression_data\.csv"|"data/miRNA_expression_data.csv"|g' \
        -e 's|"\.\./survival_data\.csv"|"data/survival_data.csv"|g' \
        -e 's|"\.\./\.\./mirna_pairs\.csv"|"../03_target_analysis/synergistic_pairs.csv"|g' \
        "$f"
done

# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------
TOTAL=$(find "$REPO" -type f \( ! -path "*/.git/*" \) | wc -l | tr -d ' ')
echo ""
echo "Done. Copied $COPIES files; $MISSING source(s) missing; $TOTAL files total in $REPO/"
echo ""
echo "REVIEW before committing:"
echo "  1. Check for hardcoded paths in copied scripts (os.chdir, /Users/...)"
echo "       grep -rn 'os.chdir\\|/Users/' $REPO --include='*.py'"
echo "  2. Confirm processed-data filenames in scripts match the public-facing"
echo "     copies after the rename pass."
echo "  3. Verify survival scripts find their CSVs at the expected relative path."
echo "  4. Spot-check the static figures shipped under figures/ — selection prefers"
echo "     v19 → v18 → v17 → v16 → v15 per figure (the v19 submission set ships all figures at v19)."
echo "  5. Update the Zenodo DOI placeholder in README.md once minted."
echo ""
echo "Then initialize / update + publish:"
echo "  cd $REPO"
echo "  git add -A && git commit -m 'v19 release for JBS submission'"
echo "  # Push, tag v1.0.0, enable Zenodo integration"
