#!/usr/bin/env python3
"""
Generate enhanced NL vs CBCA correlation plot highlighting dual-positive hits.
Panel B for the volcano plot figure.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import os

os.chdir(Path(__file__).resolve().parent)

# Set font to match volcano plot
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10

def main():
    print("Generating enhanced NL vs CBCA correlation plot...")

    # Load data
    nl_df = pd.read_csv('nl_hsa_scores.csv', index_col=0)
    cbca_df = pd.read_csv('cbca_scores.csv', index_col=0)

    # Compute ATRA CBCA reference from data (columns 0-95 are ATRA replicates)
    atra_cols = [str(i) for i in range(96)]
    atra_cbca = cbca_df.iloc[0][atra_cols].astype(float).values.mean()
    print(f"ATRA CBCA reference (from data): {atra_cbca:.4f}")

    # Calculate NL synergy and CBCA improvement for each combination
    nl_synergy_values = []
    cbca_improvement_values = []
    dual_positive_combos = []
    combo_names = []

    for combo in nl_df.index:
        # NL synergy calculation
        combo_nl_reps = nl_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
        hsa_nl_reps = nl_df.loc[combo, ['HSA_rep_1', 'HSA_rep_2', 'HSA_rep_3']].astype(float).values

        combo_nl = combo_nl_reps.mean()
        hsa_nl = hsa_nl_reps.mean()

        # NL synergy = combo - HSA (positive values = synergy, i.e., combo better than HSA)
        nl_synergy = combo_nl - hsa_nl

        # Two-sided p-value for NL
        t, p_nl = stats.ttest_ind(combo_nl_reps, hsa_nl_reps, equal_var=False)

        # CI for NL
        ci = hsa_nl / combo_nl if combo_nl > 0 else np.nan

        # CBCA improvement calculation
        try:
            combo_cbca_reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            combo_cbca = combo_cbca_reps.mean()

            # CBCA improvement = ATRA - combo (positive values = improvement)
            cbca_improvement = atra_cbca - combo_cbca

            # One-sided p-value for CBCA
            t, p_cbca_twosided = stats.ttest_1samp(combo_cbca_reps, atra_cbca)
            if combo_cbca < atra_cbca:
                p_cbca = p_cbca_twosided / 2
            else:
                p_cbca = 1 - (p_cbca_twosided / 2)

            cbca_pass = (combo_cbca < atra_cbca and p_cbca < 0.05)
        except:
            cbca_improvement = np.nan
            cbca_pass = False

        # Check if dual-positive
        is_dual_positive = (ci < 1.0 and p_nl < 0.05) and cbca_pass

        nl_synergy_values.append(nl_synergy)
        cbca_improvement_values.append(cbca_improvement)
        combo_names.append(combo)
        dual_positive_combos.append(is_dual_positive)

    # Convert to arrays
    nl_synergy_values = np.array(nl_synergy_values)
    cbca_improvement_values = np.array(cbca_improvement_values)
    dual_positive_combos = np.array(dual_positive_combos)

    # Remove NaN values
    valid_idx = ~(np.isnan(nl_synergy_values) | np.isnan(cbca_improvement_values))
    nl_synergy_values = nl_synergy_values[valid_idx]
    cbca_improvement_values = cbca_improvement_values[valid_idx]
    dual_positive_combos = dual_positive_combos[valid_idx]

    # Calculate correlation
    r, p = stats.pearsonr(nl_synergy_values, cbca_improvement_values)

    print("\nCorrelation statistics:")
    print(f"  Pearson r = {r:.3f}")
    print(f"  p-value = {p:.4f}")
    print(f"  R² = {r**2:.4f}")
    print(f"  Independent variance = {(1 - r**2)*100:.1f}%")
    print(f"  Dual-positive hits = {dual_positive_combos.sum()}")

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot background points (non-dual-positive) in gray
    non_dual_idx = ~dual_positive_combos
    ax.scatter(nl_synergy_values[non_dual_idx],
              cbca_improvement_values[non_dual_idx],
              color='gray', alpha=0.4, s=50, zorder=1,
              label='Non-significant')

    # Plot dual-positive points in orange
    dual_idx = dual_positive_combos
    ax.scatter(nl_synergy_values[dual_idx],
              cbca_improvement_values[dual_idx],
              color='#FF8C00', s=80, zorder=3,
              edgecolors='black', linewidths=0.5,
              label=f'Dual-positive (n={dual_positive_combos.sum()})')

    # Add quadrant lines
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.7, zorder=2)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.7, zorder=2)

    # Add regression line (optional, to show lack of correlation)
    # Calculate best fit line
    slope, intercept = np.polyfit(nl_synergy_values, cbca_improvement_values, 1)
    x_line = np.array([nl_synergy_values.min(), nl_synergy_values.max()])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'b-', alpha=0.3, linewidth=2, zorder=2)

    # Labels (simplified, no subtitles)
    ax.set_xlabel('NL Synergy (combo - HSA)', fontsize=24)
    ax.set_ylabel('CBCA Improvement (ATRA - combo)', fontsize=24)
    ax.tick_params(axis='both', labelsize=16)

    # Legend
    ax.legend(frameon=True, loc='lower right', fontsize=13,
             fancybox=False, edgecolor='black', framealpha=1)

    # Keep spines visible
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    # Grid
    ax.grid(True, alpha=0.3, zorder=0)

    # Save
    outdir = Path('QC_images')
    outdir.mkdir(exist_ok=True)

    fig.tight_layout()
    fig.savefig(outdir / 'NL_vs_CBCA_correlation_ENHANCED.png', dpi=300, bbox_inches='tight')
    fig.savefig(outdir / 'NL_vs_CBCA_correlation_ENHANCED.svg', bbox_inches='tight')
    plt.close()

    print(f"\n✓ Created: {outdir / 'NL_vs_CBCA_correlation_ENHANCED.png'}")
    print(f"✓ Created: {outdir / 'NL_vs_CBCA_correlation_ENHANCED.svg'}")

    # Print quadrant distribution
    q1 = (nl_synergy_values > 0) & (cbca_improvement_values > 0)  # Upper-right: dual benefit (DESIRED)
    q2 = (nl_synergy_values < 0) & (cbca_improvement_values > 0)  # Upper-left: antagonism NL, improved CBCA
    q3 = (nl_synergy_values > 0) & (cbca_improvement_values < 0)  # Lower-right: synergy NL, worse CBCA
    q4 = (nl_synergy_values < 0) & (cbca_improvement_values < 0)  # Lower-left: both worse

    print("\nQuadrant distribution:")
    print(f"  Q1 (NL synergy, CBCA improved): {q1.sum()} ({100*q1.sum()/len(nl_synergy_values):.1f}%)")
    print(f"  Q2 (NL antagonism, CBCA improved): {q2.sum()} ({100*q2.sum()/len(nl_synergy_values):.1f}%)")
    print(f"  Q3 (NL synergy, CBCA worse): {q3.sum()} ({100*q3.sum()/len(nl_synergy_values):.1f}%)")
    print(f"  Q4 (NL antagonism, CBCA worse): {q4.sum()} ({100*q4.sum()/len(nl_synergy_values):.1f}%)")
    print(f"\nDual-positive in Q1: {(dual_positive_combos & q1).sum()}/{dual_positive_combos.sum()}")

if __name__ == '__main__':
    main()
