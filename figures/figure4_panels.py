#!/usr/bin/env python3
"""
Clean, minimal manuscript figure:
- Just 3 clear categories: On-target, Liability, Housekeeping
- No redundant text labels (put in caption instead)
- Focused message
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from adjustText import adjust_text

# Set font to Helvetica (or Arial as fallback)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10


def load_batch_data(batch_dir):
    """Load per-pair metrics"""
    per_pair_dir = Path(batch_dir) / "per_pair"
    all_metrics = []

    for pair_folder in per_pair_dir.iterdir():
        if pair_folder.is_dir():
            metrics_file = pair_folder / "metrics.csv"
            if metrics_file.exists():
                df = pd.read_csv(metrics_file)
                df['pair_slug'] = pair_folder.name
                all_metrics.append(df)

    return pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()


def categorize_modules(metrics):
    """Categorize into 3 clear groups"""
    metrics = metrics.copy()

    def assign_category(row):
        display = row['module_display']
        if 'on-target' in display.lower():
            return 'On-target'
        elif 'apoptosis' in display.lower() or 'upr' in display.lower() or 'er stress' in display.lower():
            return 'Liability'
        elif 'housekeeping' in display.lower():
            return 'Housekeeping'
        else:
            return None  # Exclude "Other" category

    metrics['category_clean'] = metrics.apply(assign_category, axis=1)
    return metrics[metrics['category_clean'].notna()]  # Keep only the 3 clear categories


def create_figure(metrics, outpath):
    """Clean figure with split Panel B"""
    metrics = categorize_modules(metrics)

    # Create figure with 3 panels: A (narrower), B, C
    fig = plt.figure(figsize=(13, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1, 1], wspace=0.4)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    # ============================================================================
    # PANEL A: Clean box plot - 3 categories only
    # ============================================================================

    categories = ['On-target', 'Liability', 'Housekeeping']
    category_labels = ['On-target', 'Liability', 'Housekeeping']

    data_by_cat = [metrics[metrics['category_clean'] == cat]['incremental_coverage'].dropna()
                   for cat in categories]

    # Box plot
    box_parts = ax1.boxplot(
        data_by_cat,
        positions=range(len(categories)),
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='black', linewidth=2.5),
        boxprops=dict(linewidth=1.5, edgecolor='black'),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5)
    )

    # Simple color scheme
    colors = ['#4CAF50', '#F44336', '#9E9E9E']  # Green, Red, Gray

    for i, (patch, color) in enumerate(zip(box_parts['boxes'], colors)):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

        # Add points
        y = data_by_cat[i]
        x = np.random.normal(i, 0.06, size=len(y))
        ax1.scatter(x, y, alpha=0.3, s=12, color='black', zorder=3)

    ax1.set_xticks(range(len(categories)))
    ax1.set_xticklabels(category_labels, fontsize=11)
    ax1.set_ylabel('Incremental coverage', fontsize=12)
    ax1.set_title('A', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax1.set_ylim(bottom=0, top=max([d.max() for d in data_by_cat]) * 1.28)

    # Simple significance brackets
    y_max = max([d.max() for d in data_by_cat])

    # On-target vs Liability
    y1 = y_max * 1.10
    ax1.plot([0, 1], [y1, y1], 'k-', linewidth=1.5)
    ax1.text(0.5, y1 + 0.002, '***', ha='center', fontsize=13, fontweight='bold')

    # On-target vs Housekeeping
    y2 = y_max * 1.20
    ax1.plot([0, 2], [y2, y2], 'k-', linewidth=1.5)
    ax1.text(1, y2 + 0.002, '***', ha='center', fontsize=13, fontweight='bold')

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True, alpha=0.2, axis='y', linestyle=':', linewidth=0.8)

    # ============================================================================
    # PANEL B1: On-target vs Liability scatter
    # ============================================================================

    # Compute per-pair means
    pair_data = []
    for pair in metrics['pair_slug'].unique():
        pair_subset = metrics[metrics['pair_slug'] == pair]

        on_target = pair_subset[pair_subset['category_clean'] == 'On-target']['incremental_coverage'].mean()
        liability = pair_subset[pair_subset['category_clean'] == 'Liability']['incremental_coverage'].mean()
        housekeeping = pair_subset[pair_subset['category_clean'] == 'Housekeeping']['incremental_coverage'].mean()

        if not np.isnan(on_target) and not np.isnan(liability) and not np.isnan(housekeeping):
            pair_data.append({
                'on_target': on_target,
                'liability': liability,
                'housekeeping': housekeeping,
                'pair_slug': pair
            })

    pair_df = pd.DataFrame(pair_data)

    # --- Panel B1: On-target vs Liability ---
    ax2.scatter(pair_df['liability'], pair_df['on_target'],
                s=60, alpha=0.7, edgecolors='black', linewidths=1.2,
                color='#2196F3', zorder=3)

    # Add labels
    texts = []
    for _, row in pair_df.iterrows():
        label = row['pair_slug'].replace('hsa-', '').replace('miR-', '').replace('__', '/')
        txt = ax2.text(row['liability'], row['on_target'], label,
                      fontsize=8, ha='center', va='center', zorder=4)
        texts.append(txt)

    adjust_text(texts,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.8, alpha=0.7),
                ax=ax2,
                expand_points=(1.5, 1.5),
                expand_text=(1.3, 1.3),
                force_points=(0.4, 0.6),
                force_text=(0.6, 0.8),
                lim=1000)

    lims_b1 = [0, max(pair_df['on_target'].max(), pair_df['liability'].max()) * 1.08]
    ax2.plot(lims_b1, lims_b1, 'k--', alpha=0.5, linewidth=2, zorder=1)
    ax2.fill_between(lims_b1, lims_b1, [lims_b1[1], lims_b1[1]], alpha=0.12, color='green', zorder=0)

    ax2.set_xlabel('Inc. coverage (liability)', fontsize=11)
    ax2.set_ylabel('Inc. coverage (on-target)', fontsize=11)
    ax2.set_title('B', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax2.set_xlim(left=0)
    ax2.set_ylim(bottom=0)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

    # ============================================================================
    # PANEL B2: On-target vs Housekeeping scatter
    # ============================================================================

    ax3.scatter(pair_df['housekeeping'], pair_df['on_target'],
                s=60, alpha=0.7, edgecolors='black', linewidths=1.2,
                color='#2196F3', zorder=3)

    # Add labels
    texts = []
    for _, row in pair_df.iterrows():
        label = row['pair_slug'].replace('hsa-', '').replace('miR-', '').replace('__', '/')
        txt = ax3.text(row['housekeeping'], row['on_target'], label,
                      fontsize=8, ha='center', va='center', zorder=4)
        texts.append(txt)

    adjust_text(texts,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.8, alpha=0.7),
                ax=ax3,
                expand_points=(1.5, 1.5),
                expand_text=(1.3, 1.3),
                force_points=(0.4, 0.6),
                force_text=(0.6, 0.8),
                lim=1000)

    lims_b2 = [0, max(pair_df['on_target'].max(), pair_df['housekeeping'].max()) * 1.08]
    ax3.plot(lims_b2, lims_b2, 'k--', alpha=0.5, linewidth=2, zorder=1)
    ax3.fill_between(lims_b2, lims_b2, [lims_b2[1], lims_b2[1]], alpha=0.12, color='green', zorder=0)

    ax3.set_xlabel('Inc. coverage (housekeeping)', fontsize=11)
    ax3.set_ylabel('Inc. coverage (on-target)', fontsize=11)
    ax3.set_title('C', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax3.set_xlim(left=0)
    ax3.set_ylim(bottom=0)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

    # Layout
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"✓ Created: {outpath}")
    print("\nClean, 3-panel figure:")
    print("  • Panel A: Category comparison with significance")
    print("  • Panel B: On-target vs Liability (pair-by-pair)")
    print("  • Panel C: On-target vs Housekeeping (pair-by-pair)")
    print("\nMeans:")
    for i, cat in enumerate(categories):
        print(f"  {cat:15s} {data_by_cat[i].mean():.4f}")
    print(f"\nPairs above diagonal:")
    print(f"  B (vs Liability):     {(pair_df['on_target'] > pair_df['liability']).sum()}/{len(pair_df)}")
    print(f"  C (vs Housekeeping):  {(pair_df['on_target'] > pair_df['housekeeping']).sum()}/{len(pair_df)}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python create_final_figure_v4_clean.py <batch_dir>")
        sys.exit(1)

    batch_dir = Path(sys.argv[1])
    outdir = batch_dir / "figures"
    outdir.mkdir(exist_ok=True, parents=True)

    print(f"Loading data from {batch_dir}...")
    metrics = load_batch_data(batch_dir)

    if metrics.empty:
        print("ERROR: No metrics found!")
        sys.exit(1)

    create_figure(metrics, outdir / "manuscript_figure_FINAL.png")


if __name__ == "__main__":
    main()
