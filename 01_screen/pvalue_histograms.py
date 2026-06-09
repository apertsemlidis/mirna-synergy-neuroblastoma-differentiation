#!/usr/bin/env python3
"""
Generate p-value distribution histograms for NL and CBCA.
Two-panel figure for supplementary materials — companion to generate_qqplots.py.

Same data and p-value computation as the QQ-plot script (NL: two-sided Welch's
t-test vs HSA; CBCA: one-sided one-sample t-test vs the ATRA reference). The
histogram shows the signal-vs-noise structure (a spike near 0 over a roughly
uniform tail) without drawing a null confidence band — appropriate here because
the n=3-replicate t-tests are not exactly uniform under the null and the 946
combinations share component miRNAs (non-independence), so a Beta-derived QQ
band would be mis-calibrated.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import os

os.chdir(Path(__file__).resolve().parent)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Helvetica Neue",
    "Helvetica",
    "Arial",
    "DejaVu Sans",
]
plt.rcParams["font.size"] = 10


def hist_plot(ax, pvals, title, n_bins=20, highlight_color="#FF8C00"):
    """Histogram of p-values with the uniform-null expectation marked. The
    first bin (p < 0.05 for n_bins=20) is highlighted."""
    pvals = np.array(pvals)
    pvals = pvals[~np.isnan(pvals)]
    n = len(pvals)

    counts, _, patches = ax.hist(
        pvals,
        bins=n_bins,
        range=(0, 1),
        color="gray",
        alpha=0.55,
        edgecolor="white",
        linewidth=0.6,
        zorder=2,
    )
    # Highlight the first bin (p < 1/n_bins, i.e. p < 0.05 at 20 bins)
    patches[0].set_facecolor(highlight_color)
    patches[0].set_alpha(0.85)

    expected = n / n_bins
    ax.axhline(
        expected,
        color="red",
        ls="--",
        linewidth=1.4,
        alpha=0.8,
        label=f"Uniform null ({expected:.0f}/bin)",
        zorder=3,
    )

    n_sig = int(counts[0])
    ax.scatter(
        [],
        [],
        color=highlight_color,
        marker="s",
        s=45,
        label=f"p < 0.05 (n={n_sig}, {100 * n_sig / n:.0f}%)",
    )

    ax.set_xlabel("p-value", fontsize=16)
    ax.set_ylabel("Number of combinations", fontsize=16)
    ax.set_title(title, fontsize=18)
    ax.tick_params(axis="both", labelsize=14)
    ax.legend(fontsize=11, frameon=True, fancybox=False, edgecolor="black")
    ax.set_xlim(0, 1)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def main():
    print("Generating p-value histograms...")

    # Load data (identical to generate_qqplots.py)
    nl_df = pd.read_csv("nl_hsa_scores.csv", index_col=0)
    cbca_df = pd.read_csv("cbca_scores.csv", index_col=0)

    # ATRA reference
    atra_cols = [str(i) for i in range(96)]
    atra_cbca = cbca_df.iloc[0][atra_cols].astype(float).values.mean()

    # NL p-values (two-sided Welch's t-test vs HSA)
    nl_pvals = []
    for combo in nl_df.index:
        combo_reps = (
            nl_df.loc[combo, ["combo_rep_1", "combo_rep_2", "combo_rep_3"]]
            .astype(float)
            .values
        )
        hsa_reps = (
            nl_df.loc[combo, ["HSA_rep_1", "HSA_rep_2", "HSA_rep_3"]]
            .astype(float)
            .values
        )
        _, p = stats.ttest_ind(combo_reps, hsa_reps, equal_var=False)
        nl_pvals.append(p)

    # CBCA p-values (one-sided one-sample t-test vs ATRA)
    cbca_pvals = []
    for combo in cbca_df.index:
        try:
            combo_reps = (
                cbca_df.loc[combo, ["combo_rep_1", "combo_rep_2", "combo_rep_3"]]
                .astype(float)
                .values
            )
            combo_mean = combo_reps.mean()
            _, p_two = stats.ttest_1samp(combo_reps, atra_cbca)
            if combo_mean < atra_cbca:
                p_one = p_two / 2
            else:
                p_one = 1 - (p_two / 2)
            cbca_pvals.append(p_one)
        except Exception:
            pass

    print(f"NL p-values: {len(nl_pvals)} (sig: {sum(1 for p in nl_pvals if p < 0.05)})")
    print(
        f"CBCA p-values: {len(cbca_pvals)} (sig: {sum(1 for p in cbca_pvals if p < 0.05)})"
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    hist_plot(ax1, nl_pvals, "Neurite Length (two-sided)")
    hist_plot(ax2, cbca_pvals, "Cell Body Cluster Area (one-sided)")
    fig.tight_layout(w_pad=3)

    outdir = Path("QC_images")
    outdir.mkdir(exist_ok=True)
    fig.savefig(outdir / "pvalue_histograms_NL_CBCA.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "pvalue_histograms_NL_CBCA.svg", bbox_inches="tight")
    plt.close()

    print(f"\n✓ Created: {outdir / 'pvalue_histograms_NL_CBCA.png'}")
    print(f"✓ Created: {outdir / 'pvalue_histograms_NL_CBCA.svg'}")


if __name__ == "__main__":
    main()
