#!/usr/bin/env python3
"""
Additional file 3 (v15) — MYCN-stratified Kaplan-Meier survival analysis
for four synergistic miRNA pairs.

Layout: 4 rows x 2 columns. Each row is one pair; the two columns within
the row split patients by MYCN amplification status (Non-Amp | Amp).
Per-stratum log-rank p-values are Bonferroni-corrected (× 2 for the two
strata tested per pair).

Re-render path: each panel drawn natively (not stitched from
per-pair PNGs). Per-pair plotting logic mirrors
survival/km_mycn_stratified/km_mycn_stratified_all_v14.py.

Pairs (panel ordering matches Figure 6 v15):
  A: miR-124-3p + miR-34b-5p
  B: miR-137-3p + miR-450b-5p
  C: miR-19b-3p + miR-34b-5p
  D: miR-124-3p + miR-363-3p

Output: figures/Additional file 3 v15.{png,pdf}
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "04_survival" / "data"
OUT_DIR = ROOT / "figures"

PAIRS = [
    ("A", "124+34b", ["hsa-miR-124-3p", "hsa-miR-34b-5p"]),
    ("B", "137+450b", ["hsa-miR-137-3p", "hsa-miR-450b-5p"]),
    ("C", "19b+34b", ["hsa-miR-19b-3p", "hsa-miR-34b-5p"]),
    ("D", "124+363", ["hsa-miR-124-3p", "hsa-miR-363-3p"]),
]

COLORS = {"not_both_high": "#0072B2", "both_high": "#D55E00"}
N_STRATA = 2

DETECTION_CUTOFF = {
    "hsa-miR-200a-3p": 1.0,
    "hsa-miR-211-5p": 1.0,
    "hsa-miR-429": 1.0,
    "hsa-miR-449a": 1.0,
    "hsa-miR-449b-5p": 1.0,
}


def short_name(mirna: str) -> str:
    return mirna.replace("hsa-", "")


def classify_high(series: pd.Series, mirna: str) -> pd.Series:
    if mirna in DETECTION_CUTOFF:
        return (series > DETECTION_CUTOFF[mirna]).astype(int)
    return (series >= series.median()).astype(int)


def load_data() -> pd.DataFrame:
    expr = pd.read_csv(DATA_DIR / "miRNA_expression_data.csv")
    surv = pd.read_csv(DATA_DIR / "survival_data.csv")
    df = surv.merge(expr, on="patient_id")
    df["MYCN_amp"] = df["mycn_amplified_4.0"].astype(int)
    return df


def plot_one_stratum(ax, subset: pd.DataFrame, mycn_label: str, is_left: bool) -> None:
    """Draw one MYCN stratum subpanel."""
    grp0 = subset[subset["both_high"] == 0]
    grp1 = subset[subset["both_high"] == 1]
    n0, n1 = len(grp0), len(grp1)

    if n1 == 0:
        kmf0 = KaplanMeierFitter()
        kmf0.fit(grp0["survival_time"], grp0["event"], label=f"Not both high (n={n0})")
        kmf0.plot_survival_function(
            ax=ax, color=COLORS["not_both_high"], linewidth=1.8, ci_show=False
        )
        ax.text(
            0.97,
            0.97,
            "Both high: n=0",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="gray",
                alpha=0.85,
                linewidth=0.5,
            ),
        )
    elif n0 == 0:
        kmf1 = KaplanMeierFitter()
        kmf1.fit(grp1["survival_time"], grp1["event"], label=f"Both high (n={n1})")
        kmf1.plot_survival_function(
            ax=ax, color=COLORS["both_high"], linewidth=1.8, ci_show=False
        )
    else:
        kmf0 = KaplanMeierFitter()
        kmf1 = KaplanMeierFitter()
        kmf0.fit(grp0["survival_time"], grp0["event"], label=f"Not both high (n={n0})")
        kmf1.fit(grp1["survival_time"], grp1["event"], label=f"Both high (n={n1})")
        kmf0.plot_survival_function(
            ax=ax, color=COLORS["not_both_high"], linewidth=1.8, ci_show=False
        )
        kmf1.plot_survival_function(
            ax=ax, color=COLORS["both_high"], linewidth=1.8, ci_show=False
        )

        lr = logrank_test(
            grp0["survival_time"], grp1["survival_time"], grp0["event"], grp1["event"]
        )
        pval_adj = min(lr.p_value * N_STRATA, 1.0)
        pstr = f"p = {pval_adj:.3f}" if pval_adj >= 0.001 else "p < 0.001"
        pstr += f"\n(Bonferroni × {N_STRATA})"
        ax.text(
            0.97,
            0.97,
            pstr,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="gray",
                alpha=0.85,
                linewidth=0.5,
            ),
        )

    ax.set_title(mycn_label, fontsize=9, fontweight="bold")
    ax.set_xlabel("Time (days)", fontsize=8)
    if is_left:
        ax.set_ylabel("Survival Probability", fontsize=8)
    else:
        ax.set_yticklabels([])
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(left=0)
    ax.legend(loc="lower right", fontsize=6.5, frameon=True, framealpha=0.9)
    ax.tick_params(axis="both", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    df = load_data()

    fig, axes = plt.subplots(4, 2, figsize=(9.5, 13), sharey=True)
    for row, (letter, pair_short, mirna_pair) in enumerate(PAIRS):
        work = df.copy()
        for m in mirna_pair:
            work[f"{m}_high"] = classify_high(work[m], m)
        work["both_high"] = (
            (work[f"{mirna_pair[0]}_high"] == 1) & (work[f"{mirna_pair[1]}_high"] == 1)
        ).astype(int)
        title = f"{short_name(mirna_pair[0])} + {short_name(mirna_pair[1])}"

        ax_left = axes[row, 0]
        ax_right = axes[row, 1]

        plot_one_stratum(
            ax_left, work[work["MYCN_amp"] == 0], "MYCN Non-Amp", is_left=True
        )
        plot_one_stratum(
            ax_right, work[work["MYCN_amp"] == 1], "MYCN Amp", is_left=False
        )

        # Row-level pair title — sits above the two MYCN subpanels.
        # Anchored above ax_left so its left edge aligns with the row.
        ax_left.text(
            -0.18,
            1.30,
            f"({letter})",
            transform=ax_left.transAxes,
            fontsize=14,
            fontweight="bold",
            va="bottom",
            ha="left",
        )
        ax_left.text(
            1.06,
            1.30,
            title,
            transform=ax_left.transAxes,
            fontsize=11,
            fontweight="bold",
            va="bottom",
            ha="center",
        )

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.65, wspace=0.06, top=0.96)

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Additional file 3 v15.png"
    out_pdf = OUT_DIR / "Additional file 3 v15.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
