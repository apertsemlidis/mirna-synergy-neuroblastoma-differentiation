#!/usr/bin/env python3
"""
Generate QQ plots for NL and CBCA p-value distributions.
Two-panel figure for supplementary materials.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import os

os.chdir(Path(__file__).resolve().parent)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10

def qq_plot(ax, pvals, title, color='black', highlight_color='#FF8C00'):
    """QQ plot of observed vs expected -log10(p-values) with confidence band."""
    pvals = np.array(pvals)
    pvals = pvals[~np.isnan(pvals)]
    n = len(pvals)
    
    # Sort p-values
    sorted_p = np.sort(pvals)
    
    # Expected uniform quantiles
    expected = (np.arange(1, n + 1) - 0.5) / n
    
    # -log10 transform
    obs = -np.log10(sorted_p)
    exp = -np.log10(expected)
    
    # 95% confidence band under the null (order statistics of uniform)
    ci_lower = np.zeros(n)
    ci_upper = np.zeros(n)
    for i in range(n):
        # Beta distribution for i-th order statistic of Uniform(0,1)
        a = i + 1
        b = n - i
        ci_lower[i] = -np.log10(stats.beta.ppf(0.975, a, b))
        ci_upper[i] = -np.log10(stats.beta.ppf(0.025, a, b))
    
    # Fill confidence band
    ax.fill_between(exp, ci_lower, ci_upper, alpha=0.15, color='gray', 
                     label='95% CI (null)')
    
    # Determine significance threshold
    sig_threshold = -np.log10(0.05)
    sig_mask = obs > sig_threshold
    
    # Plot non-significant points
    ax.scatter(exp[~sig_mask], obs[~sig_mask], color='gray', s=15, alpha=0.5, zorder=2)
    
    # Plot significant points
    ax.scatter(exp[sig_mask], obs[sig_mask], color=highlight_color, s=20, alpha=0.8, 
               zorder=3, label=f'p < 0.05 (n={sig_mask.sum()})')
    
    # Diagonal
    max_val = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], 'k--', linewidth=1, alpha=0.7, label='Expected (null)')
    
    # Significance threshold line
    ax.axhline(y=sig_threshold, color='red', ls=':', linewidth=0.8, alpha=0.5)
    
    ax.set_xlabel('Expected -log₁₀(p-value)', fontsize=16)
    ax.set_ylabel('Observed -log₁₀(p-value)', fontsize=16)
    ax.set_title(title, fontsize=18)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(fontsize=11, frameon=True, fancybox=False, edgecolor='black')
    ax.set_xlim(-0.1, max_val)
    ax.set_ylim(-0.1, max_val)
    ax.set_aspect('equal')
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

def main():
    print("Generating QQ plots...")
    
    # Load data
    nl_df = pd.read_csv('nl_hsa_scores.csv', index_col=0)
    cbca_df = pd.read_csv('cbca_scores.csv', index_col=0)
    
    # ATRA reference
    atra_cols = [str(i) for i in range(96)]
    atra_cbca = cbca_df.iloc[0][atra_cols].astype(float).values.mean()
    
    # Compute NL p-values (two-sided Welch's t-test vs HSA)
    nl_pvals = []
    for combo in nl_df.index:
        combo_reps = nl_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
        hsa_reps = nl_df.loc[combo, ['HSA_rep_1', 'HSA_rep_2', 'HSA_rep_3']].astype(float).values
        _, p = stats.ttest_ind(combo_reps, hsa_reps, equal_var=False)
        nl_pvals.append(p)
    
    # Compute CBCA p-values (one-sided t-test vs ATRA)
    cbca_pvals = []
    for combo in cbca_df.index:
        try:
            combo_reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            combo_mean = combo_reps.mean()
            _, p_two = stats.ttest_1samp(combo_reps, atra_cbca)
            if combo_mean < atra_cbca:
                p_one = p_two / 2
            else:
                p_one = 1 - (p_two / 2)
            cbca_pvals.append(p_one)
        except:
            pass
    
    print(f"NL p-values: {len(nl_pvals)} (sig: {sum(1 for p in nl_pvals if p < 0.05)})")
    print(f"CBCA p-values: {len(cbca_pvals)} (sig: {sum(1 for p in cbca_pvals if p < 0.05)})")
    
    # Create two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    qq_plot(ax1, nl_pvals, 'Neurite Length (two-sided)')
    qq_plot(ax2, cbca_pvals, 'Cell Body Cluster Area (one-sided)')
    
    fig.tight_layout(w_pad=3)
    
    outdir = Path('QC_images')
    outdir.mkdir(exist_ok=True)
    fig.savefig(outdir / 'qqplots_NL_CBCA.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✓ Created: {outdir / 'qqplots_NL_CBCA.png'}")

if __name__ == '__main__':
    main()
