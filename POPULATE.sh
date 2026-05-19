#!/bin/bash
# Populate the GitHub repo with manuscript-relevant files only.
# Run from the project root:
#   bash github_repo/POPULATE.sh
#
# Source naming: project uses versioned filenames (_v14.py, _v15.py).
# Public-facing names drop the version suffix — the repo represents a
# single submission snapshot.
#
# Updated for v15 (2026-05-18): Cox forest + multivariate + time-split
# survival; NB-specific module analysis; synergy features; native re-
# render composites for Figs 3, 5, 6, 7, 8; Additional files 2, 3.

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

# Wipe previously populated subdirs to avoid leaving stale v4-era files.
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
cpy "$PAPER/scripts/screen/screen_analysis_v14.py"                  "$REPO/01_screen/screen_analysis.py"
cpy "$PAPER/scripts/screen/fxns.py"                                  "$REPO/01_screen/screen_helpers.py"
cpy "$PAPER/scripts/screen/generate_nl_volcano_dual_v14.py"          "$REPO/01_screen/volcano_plot.py"
cpy "$PAPER/scripts/screen/generate_nl_cbca_correlation_panel_v14.py" "$REPO/01_screen/nl_cbca_correlation.py"
cpy "$PAPER/scripts/screen/generate_qqplots_v14.py"                  "$REPO/01_screen/qqplots.py"
# Screen scores + family annotations (input to volcano + family-pair filtering)
cpy "$PAPER/data/screen/NL_allcombinations_vs_HSA.csv"               "$REPO/01_screen/nl_hsa_scores.csv"
cpy "$PAPER/data/screen/CBCA_allcombinations_vs_ATRA.csv"            "$REPO/01_screen/cbca_scores.csv"
cpy "$PAPER/data/screen/miR_Family_Info.csv"                         "$REPO/01_screen/mirna_family_info.csv"
# Heatmap inputs (Figure 3) — already-processed CSVs; raw Incucyte data not included (152 MB)
cpy "$PAPER/data/screen/HSA_dfs/nl_slice(96, 126, None).csv"         "$REPO/01_screen/heatmap_data/nl_hsa.csv"
cpy "$PAPER/data/screen/HSA_dfs/cbca_slice(96, 126, None).csv"       "$REPO/01_screen/heatmap_data/cbca_hsa.csv"
cpy "$PAPER/data/screen/ABS_dfs/nl_slice(96, 126, None).csv"         "$REPO/01_screen/heatmap_data/nl_absolute.csv"
cpy "$PAPER/data/screen/ABS_dfs/cbca_slice(96, 126, None).csv"       "$REPO/01_screen/heatmap_data/cbca_absolute.csv"
cpy "$PAPER/data/screen/HSAhits_p05_cytostatic.csv"                  "$REPO/01_screen/superhits.csv"

# ---------------------------------------------------------------------
# 02_dose_response — Dose-response validation (discussed in prose; no figure)
# ---------------------------------------------------------------------
echo "[02_dose_response]"
cpy "$PAPER/scripts/dose_response/process_dose_response_v14.py"      "$REPO/02_dose_response/process_dose_response.py"
# Raw plate data omitted (large); processed outputs included
cpy "$PAPER/data/dose_response/output"                                "$REPO/02_dose_response/output"

# ---------------------------------------------------------------------
# 03_target_analysis — Target-space complementarity (Figures 4, 7, 8)
# ---------------------------------------------------------------------
echo "[03_target_analysis]"
cpy "$PAPER/scripts/target_analysis/mirna_synergy_coverage_complementarity_v14.py" "$REPO/03_target_analysis/target_complementarity.py"
cpy "$PAPER/scripts/target_analysis/mirna_synergy_batch_mirtarbase_v14.py"          "$REPO/03_target_analysis/batch_analysis.py"
cpy "$PAPER/scripts/target_analysis/analyze_reinforce_dilute_v14.py"                "$REPO/03_target_analysis/statistical_tests.py"
cpy "$PAPER/scripts/target_analysis/nb_specific_analysis_v14.py"                    "$REPO/03_target_analysis/nb_specific_analysis.py"
cpy "$PAPER/scripts/target_analysis/synergy_features_v14.py"                        "$REPO/03_target_analysis/synergy_features.py"
cpy "$PAPER/scripts/target_analysis/compare_databases_v14.py"                       "$REPO/03_target_analysis/compare_databases.py"
cpy "$PAPER/scripts/target_analysis/mirna_synergy_ontology_venn_v14.py"             "$REPO/03_target_analysis/ontology_venn.py"
# Inputs
cpy "$PAPER/mirna_pairs.csv"                                                         "$REPO/03_target_analysis/synergistic_pairs.csv"
# External databases intentionally NOT copied — users download from canonical
# sources to keep licensing/provenance clean. Both are listed in .gitignore.
#   - TargetScan v7.2: https://www.targetscan.org/  -> 03_target_analysis/external/targetscan72_hsa.tsv
#   - miRTarBase MTI:  https://mirtarbase.cuhk.edu.cn/ -> 03_target_analysis/external/mirtarbase_mti.csv
cpy "$PAPER/data/external/nb_signatures"                                             "$REPO/03_target_analysis/external/nb_signatures"
# Outputs
cpy "$PAPER/data/target_analysis/batch_ts/pairs_summary.csv"                         "$REPO/03_target_analysis/outputs/pairs_summary.csv"
cpy "$PAPER/data/target_analysis/batch_ts_nb_specific/all_pairs_nb_metrics.csv"      "$REPO/03_target_analysis/outputs/all_pairs_nb_metrics.csv"
cpy "$PAPER/data/target_analysis/synergy_features/all_features.csv"                  "$REPO/03_target_analysis/outputs/all_features.csv"
# Per-pair metrics + module-term lists (input to figure4 composite)
cpy "$PAPER/data/target_analysis/batch_ts/per_pair"                                   "$REPO/03_target_analysis/outputs/per_pair"

# ---------------------------------------------------------------------
# 04_survival — Patient survival (Figure 6, Additional files 2, 3)
# Source: GSE155945 (Misiak et al. 2021), 96 NB tumors
# ---------------------------------------------------------------------
echo "[04_survival]"
cpy "$PAPER/survival/km_3group/km_3group_all_v14.py"                     "$REPO/04_survival/km_3group.py"
cpy "$PAPER/survival/km_3group/km_3group_screen_v14.py"                  "$REPO/04_survival/km_3group_screen.py"
cpy "$PAPER/survival/km_mycn_stratified/km_mycn_stratified_all_v14.py"   "$REPO/04_survival/km_mycn_stratified.py"
cpy "$PAPER/survival/km_mycn_stratified/km_mycn_stratified_screen_v14.py" "$REPO/04_survival/km_mycn_stratified_screen.py"
cpy "$PAPER/survival/cox_forest/cox_forest_all_v14.py"                   "$REPO/04_survival/cox_forest_per_pair.py"
cpy "$PAPER/survival/cox_forest/cox_forest_combined_v14.py"              "$REPO/04_survival/cox_forest_combined.py"
cpy "$PAPER/survival/cox_multivariate/multivariate_survival_v14.py"      "$REPO/04_survival/cox_multivariate.py"
cpy "$PAPER/survival/cox_time_split/cox_time_split_table_v14.py"         "$REPO/04_survival/cox_time_split.py"
cpy "$PAPER/survival/cox_diagnostic/cox_diagnostic_137+450b_v14.py"      "$REPO/04_survival/cox_diagnostic.py"
# Data (GSE155945 derived)
cpy "$PAPER/survival/miRNA_expression_data.csv"                          "$REPO/04_survival/data/miRNA_expression_data.csv"
cpy "$PAPER/survival/survival_data.csv"                                  "$REPO/04_survival/data/survival_data.csv"

# ---------------------------------------------------------------------
# figures — Composite-assembly scripts (v15) for manuscript figures
# ---------------------------------------------------------------------
echo "[figures]"
cpy "$PAPER/scripts/figures/screen_heatmaps_composite_v15.py"            "$REPO/figures/figure3_screen_heatmaps.py"
cpy "$PAPER/scripts/figures/target_complementarity_composite_v15.py"     "$REPO/figures/figure4_target_complementarity.py"
cpy "$PAPER/scripts/figures/create_final_figure_v14.py"                  "$REPO/figures/figure4_panels.py"
cpy "$PAPER/scripts/figures/nl_volcano_correlation_composite_v15.py"     "$REPO/figures/figure5_volcano_correlation.py"
cpy "$PAPER/scripts/figures/km_3group_composite_v15.py"                  "$REPO/figures/figure6_km_3group.py"
cpy "$PAPER/scripts/figures/synergy_features_composite_v15.py"           "$REPO/figures/figure7_synergy_features.py"
cpy "$PAPER/scripts/figures/nb_specific_modules_composite_v15.py"        "$REPO/figures/figure8_nb_specific.py"
cpy "$PAPER/scripts/figures/cox_forest_combined_v15.py"                  "$REPO/figures/additional_file_2_cox_forest.py"
cpy "$PAPER/scripts/figures/km_mycn_stratified_composite_v15.py"         "$REPO/figures/additional_file_3_km_mycn.py"
cpy "$PAPER/scripts/figures/figure_style.py"                             "$REPO/figures/figure_style.py"

# ---------------------------------------------------------------------
# Patch filename + path references in the populated scripts
#
# The repo ships CSVs under cleaner public-facing names than the project
# uses internally. Scripts copied above still reference the original
# names; patch them to use the public names, and rewrite the path
# resolution so scripts find their data when run from the public repo.
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

# Global literal renames applied across every .py file in the repo.
RENAMES=(
    "s|NL_allcombinations_vs_HSA\.csv|nl_hsa_scores.csv|g"
    "s|CBCA_allcombinations_vs_ATRA\.csv|cbca_scores.csv|g"
    "s|miR_Family_Info\.csv|mirna_family_info.csv|g"
    "s|HSAhits_p05_cytostatic\.csv|superhits.csv|g"
    "s|miRTarBase_MTI\.csv|mirtarbase_mti.csv|g"
)
while IFS= read -r -d '' pyfile; do
    for expr in "${RENAMES[@]}"; do
        SED_INPLACE -e "$expr" "$pyfile"
    done
done < <(find "$REPO" -name "*.py" -print0)

# Path-resolution patches: per-script, because the project's nested
# `data/screen`, `survival/<subdir>` layout collapses into a flat
# `01_screen/`, `04_survival/` in the public repo.

# 01_screen/: scripts that os.chdir into data/screen now chdir to their own folder
for f in "$REPO/01_screen/volcano_plot.py" \
         "$REPO/01_screen/nl_cbca_correlation.py" \
         "$REPO/01_screen/qqplots.py"; do
    [ -f "$f" ] && SED_INPLACE \
        -e 's|os\.chdir(Path(__file__)\.resolve()\.parents\[2\] / "data" / "screen")|os.chdir(Path(__file__).resolve().parent)|g' \
        "$f"
done

# 04_survival/: data CSVs live in 04_survival/data/, not survival/ parent
for f in "$REPO/04_survival"/*.py; do
    [ -f "$f" ] && SED_INPLACE \
        -e 's|"\.\./miRNA_expression_data\.csv"|"data/miRNA_expression_data.csv"|g' \
        -e 's|"\.\./survival_data\.csv"|"data/survival_data.csv"|g' \
        -e 's|"\.\./\.\./mirna_pairs\.csv"|"../03_target_analysis/synergistic_pairs.csv"|g' \
        "$f"
done

# cox_forest_combined.py: drop `_v14` suffix from its savefig outputs so the
# wrapper (Additional file 2) finds them under the public-facing name.
F_COX="$REPO/04_survival/cox_forest_combined.py"
if [ -f "$F_COX" ]; then
    SED_INPLACE \
        -e 's|cox_forest_combined_v14\.png|cox_forest_combined.png|g' \
        -e 's|cox_forest_combined_v14\.svg|cox_forest_combined.svg|g' \
        -e 's|cox_forest_combined_v14\.pdf|cox_forest_combined.pdf|g' \
        "$F_COX"
fi

# figures/figure3_screen_heatmaps.py: ROOT/SCREEN repointed; SLICE_TAG removed
# (heatmap CSVs are now named like `nl_hsa.csv`, `cbca_absolute.csv` etc., not by time-slice)
F3="$REPO/figures/figure3_screen_heatmaps.py"
if [ -f "$F3" ]; then
    SED_INPLACE \
        -e 's|ROOT = Path(__file__)\.resolve()\.parents\[2\]|ROOT = Path(__file__).resolve().parents[1]|g' \
        -e 's|SCREEN = ROOT / "data" / "screen"|SCREEN = ROOT / "01_screen"|g' \
        -e 's|SLICE_TAG = "slice(96, 126, None)"||g' \
        -e 's|SCREEN / "ABS_dfs" / f"{ini}_{SLICE_TAG}\.csv"|SCREEN / "heatmap_data" / f"{ini}_absolute.csv"|g' \
        -e 's|SCREEN / "HSA_dfs" / f"{ini}_{SLICE_TAG}\.csv"|SCREEN / "heatmap_data" / f"{ini}_hsa.csv"|g' \
        -e 's|SCREEN / "HSAhits_p05_cytostatic\.csv"|SCREEN / "superhits.csv"|g' \
        "$F3"
fi

# figures/figure5_volcano_correlation.py: ROOT/SCREEN repointed
F5="$REPO/figures/figure5_volcano_correlation.py"
if [ -f "$F5" ]; then
    SED_INPLACE \
        -e 's|ROOT = Path(__file__)\.resolve()\.parents\[2\]|ROOT = Path(__file__).resolve().parents[1]|g' \
        -e 's|SCREEN = ROOT / "data" / "screen"|SCREEN = ROOT / "01_screen"|g' \
        "$F5"
fi

# Every other figure composite uses `parents[2]` — collapse to `parents[1]`
# (the public repo is one level shallower than the project layout).
for f in "$REPO"/figures/figure*.py "$REPO"/figures/additional_file_*.py; do
    [ -f "$f" ] && SED_INPLACE \
        -e 's|ROOT = Path(__file__)\.resolve()\.parents\[2\]|ROOT = Path(__file__).resolve().parents[1]|g' \
        "$f"
done

# Survival figure composites: DATA_DIR points at the survival data folder
for f in "$REPO/figures/figure6_km_3group.py" \
         "$REPO/figures/additional_file_3_km_mycn.py"; do
    [ -f "$f" ] && SED_INPLACE \
        -e 's|DATA_DIR = ROOT / "survival"|DATA_DIR = ROOT / "04_survival" / "data"|g' \
        "$f"
done

# Target-analysis figure composites: rewrite the project data path
# component to the public repo's outputs folder. Use perl (multiline)
# because some scripts split the path across multiple lines inside a
# `DATA = ( ... )` block, which line-oriented sed can't match.
for f in "$REPO"/figures/figure*.py; do
    [ -f "$f" ] && perl -0777 -i -pe '
        s|ROOT\s*/\s*"data"\s*/\s*"target_analysis"\s*/\s*"batch_ts_nb_specific"\s*/\s*"all_pairs_nb_metrics\.csv"|ROOT / "03_target_analysis" / "outputs" / "all_pairs_nb_metrics.csv"|gs;
        s|ROOT\s*/\s*"data"\s*/\s*"target_analysis"\s*/\s*"synergy_features"\s*/\s*"all_features\.csv"|ROOT / "03_target_analysis" / "outputs" / "all_features.csv"|gs;
        s|ROOT\s*/\s*"data"\s*/\s*"target_analysis"\s*/\s*"batch_ts"\s*/\s*"pairs_summary\.csv"|ROOT / "03_target_analysis" / "outputs" / "pairs_summary.csv"|gs;
        # Collapse a `DATA = ( <single-segment-ROOT>... )` into a one-liner if all that remains is one path
        s|DATA\s*=\s*\(\s*(ROOT[^()]+?)\s*\)|DATA = $1|gs;
    ' "$f"
done

# Figure 4 wrapper: imports from project's scripts/figures; redirect to public figures/
F4="$REPO/figures/figure4_target_complementarity.py"
if [ -f "$F4" ]; then
    SED_INPLACE \
        -e 's|SCRIPTS_FIGURES = ROOT / "scripts" / "figures"|SCRIPTS_FIGURES = ROOT / "figures"|g' \
        -e 's|from create_final_figure_v14 import|from figure4_panels import|g' \
        -e 's|BATCH_DIR = ROOT / "data" / "target_analysis" / "batch_ts"|BATCH_DIR = ROOT / "03_target_analysis" / "outputs"|g' \
        "$F4"
fi

# Additional file 2 wrapper: source-side path points back at the project's
# survival/cox_forest folder, which doesn't exist in the public repo.
# Redirect to 04_survival/ and drop the `_v14` suffix from the source script
# + output filenames (those are renamed in the public repo).
A2="$REPO/figures/additional_file_2_cox_forest.py"
if [ -f "$A2" ]; then
    SED_INPLACE \
        -e 's|ROOT / "survival" / "cox_forest"|ROOT / "04_survival"|g' \
        -e 's|cox_forest_combined_v14\.py|cox_forest_combined.py|g' \
        -e 's|cox_forest_combined_v14\.png|cox_forest_combined.png|g' \
        -e 's|cox_forest_combined_v14\.svg|cox_forest_combined.svg|g' \
        -e 's|cox_forest_combined_v14\.pdf|cox_forest_combined.pdf|g' \
        "$A2"
fi

# figure3 heatmap composite: drop now-orphaned SLICE_TAG assignment line
# (its references were already rewritten to use the renamed heatmap_data CSVs)
F3="$REPO/figures/figure3_screen_heatmaps.py"
if [ -f "$F3" ]; then
    SED_INPLACE \
        -e '/^SLICE_TAG = ""\s*$/d' \
        -e '/^SLICE_TAG = $/d' \
        "$F3"
fi

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
echo "  2. Confirm processed-data filenames in scripts still match the renamed copies"
echo "     (e.g., screen scripts expect 'NL_allcombinations_vs_HSA.csv'; we ship as"
echo "     'nl_hsa_scores.csv' — either rename in scripts or rename in copies)."
echo "  3. Verify survival scripts find their CSVs at the expected relative path."
echo "  4. Update the Zenodo DOI placeholder in README.md once minted."
echo ""
echo "Then initialize and publish:"
echo "  cd $REPO"
echo "  git init && git add -A && git commit -m 'Initial release for JBS submission'"
echo "  # Create private GitHub repo, push, tag v1.0.0, enable Zenodo integration"
