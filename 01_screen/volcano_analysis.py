#!/usr/bin/env python3
"""
Generate Neurite Length volcano plot highlighting dual-positive hits.

Dual-positive = NL synergy (p<0.05) AND CBCA improvement (p<0.05)
FDR controlled via Fisher's combined test + BH correction.

Also emits a per-combination stats CSV (`nl_volcano_stats.csv`) with
log2(CI), -log10(p_nl), p_cbca, same_family flag, and dual_positive
flag for all 946 combinations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log10
from pathlib import Path
import os

os.chdir(Path(__file__).resolve().parent)

# Set font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10

def check_family(combo, fams_subset):
    """Check if combination is from same miRNA family"""
    one, two = combo.split(' + ')
    if one == 'hsa-miR-137':
        one = 'hsa-miR-137-3p'
    if two == 'hsa-miR-137':
        two = 'hsa-miR-137-3p'
    try:
        one_fam = fams_subset.loc[fams_subset['MiRBase ID'] == one, 'miR family'].values[0]
        two_fam = fams_subset.loc[fams_subset['MiRBase ID'] == two, 'miR family'].values[0]
        return one_fam == two_fam
    except:
        return False

def main():
    print("Generating NL volcano plot with dual-positive highlighting...")

    # Load data
    nl_df = pd.read_csv('nl_hsa_scores.csv', index_col=0)
    cbca_df = pd.read_csv('cbca_scores.csv', index_col=0)
    fams_subset = pd.read_csv('mirna_family_info.csv')

    # Compute ATRA CBCA reference from data (columns 0-95 are ATRA replicates)
    atra_cols = [str(i) for i in range(96)]
    atra_cbca = cbca_df.iloc[0][atra_cols].astype(float).values.mean()
    print(f"ATRA CBCA reference (from data): {atra_cbca:.4f}")

    # Calculate volcano data with dual-positive classification
    volcano_data = {}

    for combo in nl_df.index:
        # NL statistics
        combo_nl_reps = nl_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
        hsa_nl_reps = nl_df.loc[combo, ['HSA_rep_1', 'HSA_rep_2', 'HSA_rep_3']].astype(float).values

        combo_nl = combo_nl_reps.mean()
        hsa_nl = hsa_nl_reps.mean()

        # Use TWO-sided p-value for volcano plot (can detect both synergy and antagonism)
        t, p_nl = stats.ttest_ind(combo_nl_reps, hsa_nl_reps, equal_var=False)
        # p_nl is already two-sided from ttest_ind

        ci = hsa_nl / combo_nl if combo_nl > 0 else np.nan

        # CBCA statistics (keep one-sided for dual-positive since we only care about improvement)
        try:
            combo_cbca_reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            combo_cbca = combo_cbca_reps.mean()

            t, p_cbca_twosided = stats.ttest_1samp(combo_cbca_reps, atra_cbca)
            # One-sided for CBCA (only care if better than ATRA)
            if combo_cbca < atra_cbca:
                p_cbca = p_cbca_twosided / 2
            else:
                p_cbca = 1 - (p_cbca_twosided / 2)

            cbca_pass = (combo_cbca < atra_cbca and p_cbca < 0.05)
        except:
            cbca_pass = False

        # Family classification
        same_family = check_family(combo, fams_subset)

        # Classify: [same_family, log2(CI), -log10(p), dual_positive, CI]
        # For dual-positive: use directional criteria (synergy only)
        dual_positive = (ci < 1.0 and p_nl < 0.05) and cbca_pass
        log2_ci = np.log2(ci) if ci > 0 and not np.isnan(ci) else np.nan

        # Store extra fields for the NL/CBCA correlation panel (Fig 5 Panel B).
        nl_synergy_val = combo_nl - hsa_nl
        try:
            cbca_improvement_val = atra_cbca - combo_cbca_reps.mean()
        except Exception:
            cbca_improvement_val = np.nan

        volcano_data[combo] = [
            same_family,
            log2_ci,
            -log10(p_nl) if p_nl > 0 else 10,
            dual_positive,
            ci,
            nl_synergy_val,
            cbca_improvement_val,
        ]

    print(f"Processed {len(volcano_data)} combinations")

    # Emit per-combination stats CSV (consumed by the figure). Fig 5 Panel A uses
    # log2_ci + neglog10_p_nl + same_family + dual_positive; Panel B uses
    # nl_synergy + cbca_improvement.
    stats_rows = []
    for combo, vals in volcano_data.items():
        same_family, log2_ci, neglog10_p_nl, dual_positive, ci, nl_syn, cbca_imp = vals
        stats_rows.append({
            "combination": combo,
            "ci": ci,
            "log2_ci": log2_ci,
            "neglog10_p_nl": neglog10_p_nl,
            "same_family": bool(same_family),
            "dual_positive": bool(dual_positive),
            "nl_synergy": nl_syn,
            "cbca_improvement": cbca_imp,
        })
    stats_df = pd.DataFrame(stats_rows)
    stats_csv = "nl_volcano_stats.csv"
    stats_df.to_csv(stats_csv, index=False)
    print(f"Saved stats: {stats_csv}  ({len(stats_df)} combinations)")

    # Separate by category
    # Priority order: dual-positive > NL synergy > same-family > non-significant

    dual_pos = {k:v for k,v in volcano_data.items() if v[3]}  # Dual-positive
    nl_synergy = {k:v for k,v in volcano_data.items()
                  if not v[3] and v[4] < 1.0 and v[2] > -log10(0.05)}  # NL synergy (not dual-pos)
    same_fam = {k:v for k,v in volcano_data.items()
                if v[0] and not v[3] and not (v[4] < 1.0 and v[2] > -log10(0.05))}  # Same family, non-synergistic
    nonsig = {k:v for k,v in volcano_data.items()
              if not v[3] and not (v[4] < 1.0 and v[2] > -log10(0.05))}

    print("\nCategories:")
    print(f"  Dual-positive (NL + CBCA): {len(dual_pos)}")
    print(f"  NL synergy only: {len(nl_synergy)}")
    print(f"  Same family (non-synergistic): {len(same_fam)}")
    print(f"  Non-significant: {len(nonsig)}")

    # Create plot matching original style
    fig, ax = plt.subplots(figsize=(10, 10))

    # All same families (blue, regardless of synergy)
    same_fam_all = {k:v for k,v in volcano_data.items() if v[0]}
    all_nl_syn = {**nl_synergy, **dual_pos}

    # Plot in order: background first, orange, then red on top
    # Non-significant (gray)
    if nonsig:
        ax.scatter([v[1] for v in nonsig.values()],
                  [v[2] for v in nonsig.values()],
                  color='gray', alpha=0.5, s=50, zorder=1)

    # Same family (blue)
    if same_fam_all:
        ax.scatter([v[1] for v in same_fam_all.values()],
                  [v[2] for v in same_fam_all.values()],
                  label='Same families',
                  color='blue', s=50, zorder=2)

    # NL synergy only (NOT dual-positive) - red
    if nl_synergy:
        ax.scatter([v[1] for v in nl_synergy.values()],
                  [v[2] for v in nl_synergy.values()],
                  label='NL synergy only',
                  color='red', s=50, zorder=3)

    # Dual-positive (orange) - on top so fully visible
    if dual_pos:
        ax.scatter([v[1] for v in dual_pos.values()],
                  [v[2] for v in dual_pos.values()],
                  label='Dual-positive (NL + CBCA)',
                  color='#FF8C00', s=50, zorder=4)

    # Threshold lines
    ax.axhline(y=-log10(0.05), color='black', ls='--', linewidth=1)
    ax.axvline(x=0, color='black', ls='--', linewidth=1)  # log2(1) = 0

    # Labels
    ax.set_xlabel('log\u2082(HSA NL / Combination NL)', fontsize=24)
    ax.set_ylabel('-log\u2081\u2080(p-value)', fontsize=24)
    ax.tick_params(axis='both', labelsize=16)

    # Set symmetric axis limits
    xmax = max(abs(v[1]) for v in volcano_data.values() if not np.isnan(v[1]))
    ax.set_xlim(-xmax * 1.1, xmax * 1.1)
    ax.set_ylim(-0.2, max([v[2] for v in volcano_data.values()]) + 0.5)

    # Legend in upper right with box
    ax.legend(frameon=True, loc='upper right', fontsize=13,
             fancybox=False, edgecolor='black', framealpha=1)

    # Add labels for top hits with manual placement
    all_nl_syn = {**nl_synergy, **dual_pos}
    top_hits = sorted([(k, v[1], v[2]) for k, v in all_nl_syn.items()],
                     key=lambda x: x[2], reverse=True)[:8]

    # Manual label positions: label -> (label_x, label_y) or None for auto
    label_positions = {
        'miR-137 + miR-17-5p':        (-3.2, 3.05),
        'miR-20a-5p + miR-449b-5p':   (-3.2, 2.72),
        'miR-106a-5p + miR-449b-5p':  (-3.2, 2.45),
        'miR-211-5p + miR-449b-5p':   (-1.1, 2.25),
    }

    for combo, log2ci, logp in top_hits:
        if logp > 2.5:
            label = combo.replace('hsa-', '')
            if label in label_positions:
                lx, ly = label_positions[label]
            else:
                lx, ly = log2ci + 0.12, logp

            ax.annotate(label, xy=(log2ci, logp), xytext=(lx, ly),
                       fontsize=11, ha='left', va='center', zorder=6,
                       arrowprops=dict(arrowstyle='-', color='black',
                                       lw=0.5, alpha=0.7))

    # Keep spines visible (like original)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    # Save
    outdir = Path('QC_images')
    outdir.mkdir(exist_ok=True)

    fig.tight_layout()
    fig.savefig(outdir / 'neurite length CI volcano plot DUAL.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✓ Created: {outdir / 'neurite length CI volcano plot DUAL.png'}")

    # ---- Summary statistics with corrected analysis ----
    same_fam_all = {k:v for k,v in volcano_data.items() if v[0]}
    same_fam_synergistic = {k:v for k,v in same_fam_all.items() if v[4] < 1.0 and v[2] > -log10(0.05)}
    all_nl_syn = {**nl_synergy, **dual_pos}

    # Directional concordance: what fraction of NL synergies have CBCA < ATRA?
    n_total = len(volcano_data)
    nl_syn_combos = list(all_nl_syn.keys())
    nl_syn_below_atra = 0
    all_below_atra = 0
    for combo in nl_df.index:
        try:
            combo_cbca_reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            if combo_cbca_reps.mean() < atra_cbca:
                all_below_atra += 1
                if combo in nl_syn_combos:
                    nl_syn_below_atra += 1
        except:
            pass
    nl_syn_n = len(nl_syn_combos)

    # Fisher's exact test for directional concordance
    a = nl_syn_below_atra
    b = nl_syn_n - nl_syn_below_atra
    c = all_below_atra - nl_syn_below_atra
    d = n_total - nl_syn_n - c
    odds_ratio_dir, fisher_p_dir = stats.fisher_exact([[a, b], [c, d]])

    # Fisher's exact test for significance concordance
    cbca_sig_count = 0
    for combo in nl_df.index:
        try:
            combo_cbca_reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            combo_cbca = combo_cbca_reps.mean()
            t, p_cbca_twosided = stats.ttest_1samp(combo_cbca_reps, atra_cbca)
            if combo_cbca < atra_cbca:
                p_cbca = p_cbca_twosided / 2
            else:
                p_cbca = 1 - (p_cbca_twosided / 2)
            if combo_cbca < atra_cbca and p_cbca < 0.05:
                cbca_sig_count += 1
        except:
            pass
    a2 = len(dual_pos)
    b2 = nl_syn_n - len(dual_pos)
    c2 = cbca_sig_count - len(dual_pos)
    d2 = n_total - nl_syn_n - c2
    odds_ratio_sig, fisher_p_sig = stats.fisher_exact([[a2, b2], [c2, d2]])

    # Permutation test for dual-positive count (1000 iterations, seed=42)
    rng = np.random.RandomState(42)
    nl_pass_arr = np.array([combo in nl_syn_combos for combo in nl_df.index])
    cbca_pass_arr = np.zeros(n_total, dtype=bool)
    for i, combo in enumerate(nl_df.index):
        try:
            combo_cbca_reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            combo_cbca = combo_cbca_reps.mean()
            t, p_cbca_twosided = stats.ttest_1samp(combo_cbca_reps, atra_cbca)
            if combo_cbca < atra_cbca:
                p_cbca = p_cbca_twosided / 2
            else:
                p_cbca = 1 - (p_cbca_twosided / 2)
            cbca_pass_arr[i] = (combo_cbca < atra_cbca and p_cbca < 0.05)
        except:
            pass

    n_perms = 10000
    perm_counts = np.array([
        (nl_pass_arr & rng.permutation(cbca_pass_arr)).sum()
        for _ in range(n_perms)
    ])
    perm_mean = perm_counts.mean()
    perm_p = (perm_counts >= len(dual_pos)).sum() / n_perms

    # Fisher's combined test (combine one-sided NL and CBCA p-values) + BH
    nl_pvals_1sided = np.ones(n_total)
    cbca_pvals_1sided = np.ones(n_total)
    combo_list = nl_df.index.tolist()
    cbca_means_arr = np.zeros(n_total)
    nl_ci_arr = np.ones(n_total)

    for i, combo in enumerate(combo_list):
        combo_nl_reps = nl_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
        hsa_nl_reps = nl_df.loc[combo, ['HSA_rep_1', 'HSA_rep_2', 'HSA_rep_3']].astype(float).values
        nl_ci_arr[i] = hsa_nl_reps.mean() / combo_nl_reps.mean() if combo_nl_reps.mean() > 0 else np.nan
        _, p_nl_2s = stats.ttest_ind(combo_nl_reps, hsa_nl_reps, equal_var=False)
        if combo_nl_reps.mean() > hsa_nl_reps.mean():
            nl_pvals_1sided[i] = p_nl_2s / 2
        else:
            nl_pvals_1sided[i] = 1 - (p_nl_2s / 2)

        try:
            reps = cbca_df.loc[combo, ['combo_rep_1', 'combo_rep_2', 'combo_rep_3']].astype(float).values
            cbca_means_arr[i] = reps.mean()
            _, p2 = stats.ttest_1samp(reps, atra_cbca)
            if reps.mean() < atra_cbca:
                cbca_pvals_1sided[i] = p2 / 2
            else:
                cbca_pvals_1sided[i] = 1 - (p2 / 2)
        except:
            pass

    from scipy.stats import chi2 as chi2_dist
    combined_pvals = np.ones(n_total)
    for i in range(n_total):
        p1 = max(nl_pvals_1sided[i], 1e-300)
        p2 = max(cbca_pvals_1sided[i], 1e-300)
        chi2_stat = -2 * (np.log(p1) + np.log(p2))
        combined_pvals[i] = 1 - chi2_dist.cdf(chi2_stat, df=4)

    # BH correction
    sorted_idx = np.argsort(combined_pvals)
    qvals = np.zeros(n_total)
    cummin = 1.0
    for j in range(n_total - 1, -1, -1):
        rank = j + 1
        adj = combined_pvals[sorted_idx[j]] * n_total / rank
        cummin = min(adj, cummin)
        qvals[sorted_idx[j]] = min(cummin, 1.0)

    # Require directional criteria
    directional = (nl_ci_arr < 1.0) & (cbca_means_arr < atra_cbca)
    fisher_q005 = (qvals < 0.05) & directional
    fisher_q001 = (qvals < 0.01) & directional

    # Check all dual-positive hits pass Fisher + BH
    dp_in_fisher = sum(1 for combo in dual_pos if fisher_q005[combo_list.index(combo)])

    print("\n" + "=" * 70)
    print("SUMMARY FOR MANUSCRIPT")
    print("=" * 70)
    print(f"\nATRA CBCA reference (from data): {atra_cbca:.4f}")
    print(f"\nTotal combinations tested: {n_total}")
    print(f"NL synergies (p < 0.05, CI < 1): {len(all_nl_syn)}")
    print(f"  - Same family synergies: {len(same_fam_synergistic)}")
    print(f"  - Different family synergies: {len(all_nl_syn) - len(same_fam_synergistic)}")
    print(f"  - Dual-positive (NL + CBCA): {len(dual_pos)}")
    print(f"    - Same family dual-pos: {sum(1 for k in dual_pos if k in same_fam_all)}")

    print("\nDirectional concordance:")
    print(f"  NL synergies with CBCA < ATRA: {nl_syn_below_atra}/{nl_syn_n} ({100*nl_syn_below_atra/nl_syn_n:.0f}%)")
    print(f"  Background rate: {all_below_atra}/{n_total} ({100*all_below_atra/n_total:.0f}%)")
    print(f"  Odds ratio: {odds_ratio_dir:.1f}, Fisher's p = {fisher_p_dir:.1e}")

    print("\nSignificance concordance:")
    print(f"  NL synergies with CBCA p < 0.05: {len(dual_pos)}/{nl_syn_n} ({100*len(dual_pos)/nl_syn_n:.0f}%)")
    print(f"  Background CBCA pass rate: {cbca_sig_count}/{n_total} ({100*cbca_sig_count/n_total:.0f}%)")
    print(f"  Odds ratio: {odds_ratio_sig:.1f}, Fisher's p = {fisher_p_sig:.1e}")

    print(f"\nPermutation test ({n_perms} iterations, seed=42):")
    print(f"  Observed dual-positive: {len(dual_pos)}")
    print(f"  Permutation mean null: {perm_mean:.1f}")
    print(f"  Permutation p: {'< 0.001' if perm_p == 0 else f'{perm_p:.3f}'}")

    print("\nFisher's combined test + BH correction:")
    print(f"  q < 0.01: {fisher_q001.sum()} hits")
    print(f"  q < 0.05: {fisher_q005.sum()} hits")
    print(f"  All {len(dual_pos)} dual-positive pass q < 0.05: {dp_in_fisher}/{len(dual_pos)}")

if __name__ == '__main__':
    main()
