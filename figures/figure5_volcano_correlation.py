#!/usr/bin/env python3
"""
Figure 5 composite (v15) — Dual-phenotype filtering identifies neurite
length synergies. Two-panel native re-render (no PNG stitching):

  (A) Neurite length volcano. x = log2(HSA NL / Combination NL),
      y = -log10(p-value). Dual-positive (NL + CBCA) highlighted in
      orange; NL-synergy-only in red; same-family pairs in blue.
  (B) Independence of NL synergy and CBCA improvement. x = NL synergy
      (combo - HSA), y = CBCA improvement (ATRA - combo). Dual-positive
      hits highlighted in upper-right quadrant.

Computes both panels' underlying statistics from the screen CSVs in
data/screen/, then plots onto a shared figure with descriptor labels.

Output: figures/Figure 5 v15.{png,pdf}
"""

from math import log10
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "01_screen"
OUT_DIR = ROOT / "figures"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Helvetica Neue",
    "Helvetica",
    "Arial",
    "DejaVu Sans",
]

REP_COLS = ["combo_rep_1", "combo_rep_2", "combo_rep_3"]
HSA_COLS = ["HSA_rep_1", "HSA_rep_2", "HSA_rep_3"]


def _same_family(combo: str, fams: pd.DataFrame) -> bool:
    one, two = combo.split(" + ")
    if one == "hsa-miR-137":
        one = "hsa-miR-137-3p"
    if two == "hsa-miR-137":
        two = "hsa-miR-137-3p"
    try:
        f1 = fams.loc[fams["MiRBase ID"] == one, "miR family"].values[0]
        f2 = fams.loc[fams["MiRBase ID"] == two, "miR family"].values[0]
        return f1 == f2
    except IndexError:
        return False


def compute_panel_data(
    nl_df: pd.DataFrame, cbca_df: pd.DataFrame, fams: pd.DataFrame, atra_cbca: float
) -> pd.DataFrame:
    rows = []
    for combo in nl_df.index:
        combo_nl = nl_df.loc[combo, REP_COLS].astype(float).values
        hsa_nl = nl_df.loc[combo, HSA_COLS].astype(float).values
        combo_nl_mean = combo_nl.mean()
        hsa_nl_mean = hsa_nl.mean()
        _, p_nl = stats.ttest_ind(combo_nl, hsa_nl, equal_var=False)
        ci = hsa_nl_mean / combo_nl_mean if combo_nl_mean > 0 else np.nan
        nl_synergy = combo_nl_mean - hsa_nl_mean

        try:
            combo_cbca = cbca_df.loc[combo, REP_COLS].astype(float).values
            combo_cbca_mean = combo_cbca.mean()
            _, p_cbca_two = stats.ttest_1samp(combo_cbca, atra_cbca)
            if combo_cbca_mean < atra_cbca:
                p_cbca = p_cbca_two / 2
            else:
                p_cbca = 1 - (p_cbca_two / 2)
            cbca_imp = atra_cbca - combo_cbca_mean
            cbca_pass = (combo_cbca_mean < atra_cbca) and (p_cbca < 0.05)
        except KeyError:
            combo_cbca_mean = np.nan
            cbca_imp = np.nan
            cbca_pass = False

        dual_pos = (
            (ci is not np.nan)
            and (not np.isnan(ci))
            and ci < 1.0
            and p_nl < 0.05
            and cbca_pass
        )
        log2ci = (
            np.log2(ci) if (ci is not None and not np.isnan(ci) and ci > 0) else np.nan
        )

        rows.append(
            {
                "combo": combo,
                "same_family": _same_family(combo, fams),
                "log2ci": log2ci,
                "ci": ci,
                "neglog10p_nl": -log10(p_nl) if p_nl > 0 else 10.0,
                "nl_synergy": nl_synergy,
                "cbca_improvement": cbca_imp,
                "dual_positive": dual_pos,
            }
        )
    return pd.DataFrame(rows).set_index("combo")


def plot_volcano(ax: plt.Axes, df: pd.DataFrame) -> None:
    thresh = -log10(0.05)
    nl_syn_mask = (df["ci"] < 1.0) & (df["neglog10p_nl"] > thresh)
    dual = df["dual_positive"]
    nl_only = nl_syn_mask & (~dual)
    same_fam_all = df["same_family"]
    nonsig = (~dual) & (~nl_syn_mask)

    ax.scatter(
        df.loc[nonsig, "log2ci"],
        df.loc[nonsig, "neglog10p_nl"],
        color="gray",
        alpha=0.5,
        s=50,
        zorder=1,
    )
    ax.scatter(
        df.loc[same_fam_all, "log2ci"],
        df.loc[same_fam_all, "neglog10p_nl"],
        color="blue",
        s=50,
        zorder=2,
        label="Same families",
    )
    ax.scatter(
        df.loc[nl_only, "log2ci"],
        df.loc[nl_only, "neglog10p_nl"],
        color="red",
        s=50,
        zorder=3,
        label="NL synergy only",
    )
    ax.scatter(
        df.loc[dual, "log2ci"],
        df.loc[dual, "neglog10p_nl"],
        color="#FF8C00",
        s=50,
        zorder=4,
        label="Dual-positive (NL + CBCA)",
    )

    ax.axhline(y=thresh, color="black", ls="--", linewidth=1)
    ax.axvline(x=0, color="black", ls="--", linewidth=1)
    ax.set_xlabel("log₂(HSA NL / Combination NL)", fontsize=18)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=18)
    ax.tick_params(axis="both", labelsize=13)

    valid = df["log2ci"].dropna()
    xmax = max(abs(valid.min()), abs(valid.max()))
    ax.set_xlim(-xmax * 1.1, xmax * 1.1)
    ax.set_ylim(-0.2, df["neglog10p_nl"].max() + 0.5)
    ax.legend(
        frameon=True,
        loc="upper right",
        fontsize=11,
        fancybox=False,
        edgecolor="black",
        framealpha=1,
    )

    label_positions = {
        "miR-137 + miR-17-5p": (-3.2, 3.05),
        "miR-20a-5p + miR-449b-5p": (-3.2, 2.72),
        "miR-106a-5p + miR-449b-5p": (-3.2, 2.45),
        "miR-211-5p + miR-449b-5p": (-1.1, 2.25),
    }
    top = df[nl_syn_mask | dual].sort_values("neglog10p_nl", ascending=False).head(8)
    for combo, row in top.iterrows():
        if row["neglog10p_nl"] > 2.5:
            label = combo.replace("hsa-", "")
            lx, ly = label_positions.get(
                label, (row["log2ci"] + 0.12, row["neglog10p_nl"])
            )
            ax.annotate(
                label,
                xy=(row["log2ci"], row["neglog10p_nl"]),
                xytext=(lx, ly),
                fontsize=9,
                ha="left",
                va="center",
                zorder=6,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.5, alpha=0.7),
            )
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def plot_correlation(ax: plt.Axes, df: pd.DataFrame) -> None:
    valid = df.dropna(subset=["nl_synergy", "cbca_improvement"])
    dual = valid["dual_positive"]
    x = valid["nl_synergy"].values
    y = valid["cbca_improvement"].values
    r, _ = stats.pearsonr(x, y)

    ax.scatter(
        x[~dual],
        y[~dual],
        color="gray",
        alpha=0.4,
        s=50,
        zorder=1,
        label="Non-significant",
    )
    ax.scatter(
        x[dual],
        y[dual],
        color="#FF8C00",
        s=80,
        zorder=3,
        edgecolors="black",
        linewidths=0.5,
        label=f"Dual-positive (n={int(dual.sum())})",
    )

    ax.axhline(y=0, color="black", linestyle="--", linewidth=1, alpha=0.7, zorder=2)
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1, alpha=0.7, zorder=2)
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, slope * xs + intercept, "b-", alpha=0.3, linewidth=2, zorder=2)

    ax.set_xlabel("NL Synergy (combo − HSA)", fontsize=18)
    ax.set_ylabel("CBCA Improvement (ATRA − combo)", fontsize=18)
    ax.tick_params(axis="both", labelsize=13)
    ax.legend(
        frameon=True,
        loc="lower right",
        fontsize=11,
        fancybox=False,
        edgecolor="black",
        framealpha=1,
    )
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.grid(True, alpha=0.3, zorder=0)
    ax.text(
        0.02,
        0.98,
        f"Pearson r = {r:.2f}",
        transform=ax.transAxes,
        fontsize=11,
        va="top",
        ha="left",
    )


def main() -> None:
    nl_df = pd.read_csv(SCREEN / "nl_hsa_scores.csv", index_col=0)
    cbca_df = pd.read_csv(SCREEN / "cbca_scores.csv", index_col=0)
    fams = pd.read_csv(SCREEN / "mirna_family_info.csv")

    atra_cols = [str(i) for i in range(96)]
    atra_cbca = cbca_df.iloc[0][atra_cols].astype(float).values.mean()

    panel_data = compute_panel_data(nl_df, cbca_df, fams, atra_cbca)
    print(
        f"Computed {len(panel_data)} combinations; "
        f"dual-positive = {int(panel_data['dual_positive'].sum())}; "
        f"ATRA CBCA = {atra_cbca:.4f}"
    )

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
    plot_volcano(axL, panel_data)
    plot_correlation(axR, panel_data)

    for ax, letter in [(axL, "A"), (axR, "B")]:
        ax.text(
            -0.10,
            1.02,
            f"({letter})",
            transform=ax.transAxes,
            fontsize=20,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    OUT_DIR.mkdir(exist_ok=True)
    out_png = OUT_DIR / "Figure 5 v15.png"
    out_pdf = OUT_DIR / "Figure 5 v15.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
