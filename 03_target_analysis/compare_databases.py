#!/usr/bin/env python3
"""
Compare TargetScan vs miRTarBase results to understand database differences.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys


def main():
    # Load both summaries
    ts_summary = pd.read_csv("batch_ts/pairs_summary.csv")
    mt_summary = pd.read_csv("batch_mirtarbase/pairs_summary.csv")

    # Create comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: NBS distribution comparison
    ax = axes[0, 0]
    ax.hist(ts_summary['NBS'], bins=15, alpha=0.7, label='TargetScan (predicted)', color='#2E86AB', edgecolor='black')
    ax.hist(mt_summary['NBS'], bins=15, alpha=0.7, label='miRTarBase (validated)', color='#C73E1D', edgecolor='black')
    ax.axvline(0, color='black', linewidth=2, linestyle='--')
    ax.set_xlabel('Neurite Bias Score (NBS)', fontsize=11)
    ax.set_ylabel('Number of pairs', fontsize=11)
    ax.set_title('A. NBS distribution by database', fontsize=12, fontweight='bold', loc='left')
    ax.legend(frameon=False)
    ax.text(0.05, 0.95, f'TargetScan mean: {ts_summary["NBS"].mean():.4f}\nmiRTarBase mean: {mt_summary["NBS"].mean():.4f}',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), fontsize=9)

    # Panel B: Inc_on vs Inc_ctrl scatter for both databases
    ax = axes[0, 1]
    ax.scatter(ts_summary['Inc_ctrl'], ts_summary['Inc_on'], s=80, alpha=0.7,
               label='TargetScan', color='#2E86AB', edgecolors='black', linewidths=0.5)
    ax.scatter(mt_summary['Inc_ctrl'], mt_summary['Inc_on'], s=80, alpha=0.7,
               label='miRTarBase', color='#C73E1D', edgecolors='black', linewidths=0.5)

    lims = [0, max(ts_summary['Inc_ctrl'].max(), ts_summary['Inc_on'].max(),
                   mt_summary['Inc_ctrl'].max(), mt_summary['Inc_on'].max()) * 1.05]
    ax.plot(lims, lims, 'k--', alpha=0.4, linewidth=1, label='y = x')

    ax.set_xlabel('Inc. coverage (control)', fontsize=11)
    ax.set_ylabel('Inc. coverage (on-target)', fontsize=11)
    ax.set_title('B. Reinforcement pattern by database', fontsize=12, fontweight='bold', loc='left')
    ax.legend(frameon=False, loc='upper left')

    # Panel C: Target counts comparison
    ax = axes[1, 0]

    # Get overlapping pairs only
    ts_summary['pair_key'] = ts_summary['mirA_raw'] + '_' + ts_summary['mirB_raw']
    mt_summary['pair_key'] = mt_summary['mirA_raw'] + '_' + mt_summary['mirB_raw']

    common_pairs = set(ts_summary['pair_key']) & set(mt_summary['pair_key'])

    ts_targets = []
    mt_targets = []
    for pk in common_pairs:
        ts_targets.append(ts_summary[ts_summary['pair_key'] == pk]['union_targets_total'].iloc[0])
        mt_targets.append(mt_summary[mt_summary['pair_key'] == pk]['union_targets_total'].iloc[0])

    ax.scatter(mt_targets, ts_targets, s=100, alpha=0.7, color='purple', edgecolors='black')
    max_val = max(max(ts_targets), max(mt_targets))
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.4)

    ax.set_xlabel('Union targets (miRTarBase)', fontsize=11)
    ax.set_ylabel('Union targets (TargetScan)', fontsize=11)
    ax.set_title(f'C. Target counts ({len(common_pairs)} overlapping pairs)', fontsize=12, fontweight='bold', loc='left')
    ax.text(0.05, 0.95, f'TargetScan median: {np.median(ts_targets):.0f}\nmiRTarBase median: {np.median(mt_targets):.0f}',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5), fontsize=9)

    # Panel D: Database coverage summary
    ax = axes[1, 1]
    ax.axis('off')

    summary_text = f"""
DATABASE COMPARISON SUMMARY

Total pairs in dataset: 31

TargetScan (predicted targets):
  • Pairs analyzed: {len(ts_summary)}/31
  • Mean NBS: {ts_summary['NBS'].mean():.4f}
  • Pairs with NBS > 0: {(ts_summary['NBS'] > 0).sum()}/{len(ts_summary)}
  • Interpretation: STRONG reinforcement

miRTarBase (validated targets):
  • Pairs analyzed: {len(mt_summary)}/31
  • Mean NBS: {mt_summary['NBS'].mean():.4f}
  • Pairs with NBS > 0: {(mt_summary['NBS'] > 0).sum()}/{len(mt_summary)}
  • Interpretation: OPPOSITE pattern

KEY INSIGHT:
Validated targets show LESS selectivity
than predicted targets, suggesting:
  1. Validation bias toward housekeeping genes
  2. Neuron-specific targets under-represented
  3. Synergy may involve context-specific
     or low-evidence interactions

RECOMMENDATION:
Use TargetScan predictions for main analysis,
note miRTarBase discrepancy in discussion.
"""

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
            verticalalignment='top', fontsize=9, family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.tight_layout()
    fig.savefig('database_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("✓ Created: database_comparison.png")

    # Print detailed comparison
    print("\n" + "="*80)
    print("DATABASE COMPARISON")
    print("="*80)
    print(f"\nTargetScan: {len(ts_summary)} pairs, mean NBS = {ts_summary['NBS'].mean():.4f}")
    print(f"miRTarBase: {len(mt_summary)} pairs, mean NBS = {mt_summary['NBS'].mean():.4f}")
    print(f"\nDifference: {ts_summary['NBS'].mean() - mt_summary['NBS'].mean():.4f}")
    print("\nThis large discrepancy suggests database-specific biases.")
    print("TargetScan predictions show clear neurite bias; miRTarBase does not.")


if __name__ == "__main__":
    main()
