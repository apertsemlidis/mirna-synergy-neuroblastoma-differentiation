#!/usr/bin/env python3
"""
Step 3 of the screen-processing pipeline.

Reads complete.parquet and reproduces the two superhit tables consumed
downstream by the Figure 3 heatmaps (the "superhits outlined" rectangles)
and by the volcano plot scripts.

Pipeline:
  1. For every combo, compute mean NL, mean CBCA, mean NBP over the
     screen window (t = 96..126h, averaged over time and replicates).
  2. Keep combos where combo_CBCA < ATRA_CBCA (cytostatic) and
     combo_NL > max(one_NL, two_NL) (synergy vs HSA).
  3. Add a Welch t-test of combo NL replicates vs the winning
     single-agent's NL replicates → nl_p_value_over_HSA.
  4. Filter at p < 0.05 → HSAhits_p05.csv.
  5. Add Welch t-tests vs ATRA for CBCA and NL → cbca_pvalue_rel_ATRA,
     nl_pvalue_rel_ATRA.
  6. Filter cbca_pvalue_rel_ATRA < 0.05 → superhits.csv.

Outputs:
    data/screen/HSAhits_p05.csv
    data/screen/superhits.csv

Run `--verify` to byte-compare the freshly written CSVs against the
existing ones. Falls back to numeric diff if the CSVs differ only in
pandas-version float-repr.
"""

import argparse
import filecmp
import shutil
from pathlib import Path

import pandas as pd
import scipy.stats as stats

from screen_helpers import load_complete

ROOT = Path(__file__).resolve().parents[2]
SCREEN = ROOT / "data" / "screen"
PARQUET = SCREEN / "complete.parquet"

TIME_SLICE = slice(96, 126)
ATRA = "ATRA (25 uM)"

# Column order matches the existing HSAhits_p05.csv produced by the
# original pipeline (combo_nbp sits between "nl distance from HSA" and "one").
HSAHITS_P05_COLS = [
    "combo_cbca",
    "combo_nl",
    "nl distance from HSA",
    "combo_nbp",
    "one",
    "one_nl",
    "two",
    "two_nl",
    "nl_p_value_over_HSA",
]

# Existing superhits.csv uses a different name for the
# HSA-pvalue column (compressed "pvalue", no underscore around "value").
# This is an artefact of the original pipeline renaming between save points;
# matched so the file is byte-identical.
HSAHITS_CYTO_COLS = [
    "combo_cbca",
    "combo_nl",
    "nl distance from HSA",
    "combo_nbp",
    "one",
    "one_nl",
    "two",
    "two_nl",
    "nl_pvalue_over_HSA",
    "cbca_pvalue_rel_ATRA",
    "nl_pvalue_rel_ATRA",
]


def metric_mean(complete, condition, metric, time=TIME_SLICE):
    return complete[metric].loc[condition].loc[:, time].mean().mean()


def metric_per_replicate(complete, condition, metric, time=TIME_SLICE):
    return complete[metric].loc[condition].loc[:, time].mean(axis=1)


def winner_single_agent(complete, combo):
    one, two = combo.split(" + ")
    o = metric_mean(complete, one, "neurite length")
    t = metric_mean(complete, two, "neurite length")
    return one if o >= t else two


def compute(complete):
    """Run the full superhits pipeline. Returns (hits_p05_df, hits_cyto_df)."""
    cytostatic = metric_mean(complete, ATRA, "cell body cluster area")

    nl_index = complete["neurite length"].index.unique()
    combos = sorted([c for c in nl_index if "+" in c])

    rows = {}
    for combo in combos:
        one, two = combo.split(" + ")
        one_nl = metric_mean(complete, one, "neurite length")
        two_nl = metric_mean(complete, two, "neurite length")
        combo_nl = metric_mean(complete, combo, "neurite length")
        combo_cbca = metric_mean(complete, combo, "cell body cluster area")
        combo_nbp = metric_mean(complete, combo, "neurite branch points")
        diff_nl = combo_nl - max(one_nl, two_nl)
        if combo_cbca < cytostatic and diff_nl > 0:
            rows[combo] = {
                "combo_cbca": combo_cbca,
                "combo_nl": combo_nl,
                "nl distance from HSA": diff_nl,
                "combo_nbp": combo_nbp,
                "one": one,
                "one_nl": one_nl,
                "two": two,
                "two_nl": two_nl,
            }

    superhits = pd.DataFrame.from_dict(rows, orient="index")

    # NL t-test vs the winning single agent.
    for combo in superhits.index:
        win = winner_single_agent(complete, combo)
        combo_reps = metric_per_replicate(complete, combo, "neurite length")
        win_reps = metric_per_replicate(complete, win, "neurite length")
        superhits.loc[combo, "nl_p_value_over_HSA"] = stats.ttest_ind(
            combo_reps, win_reps, equal_var=False
        ).pvalue

    hits_p05 = superhits[superhits["nl_p_value_over_HSA"] < 0.05].copy()
    hits_p05 = hits_p05.sort_values("nl distance from HSA", ascending=False)
    hits_p05 = hits_p05[HSAHITS_P05_COLS]

    # Build the cytostatic-gated table from the p<0.05 subset, with
    # additional ATRA-relative t-tests and the renamed HSA column.
    hits_cyto = hits_p05.rename(columns={"nl_p_value_over_HSA": "nl_pvalue_over_HSA"})
    atra_cbca = metric_per_replicate(complete, ATRA, "cell body cluster area")
    atra_nl = metric_per_replicate(complete, ATRA, "neurite length")
    for combo in hits_cyto.index:
        combo_cbca_reps = metric_per_replicate(
            complete, combo, "cell body cluster area"
        )
        combo_nl_reps = metric_per_replicate(complete, combo, "neurite length")
        hits_cyto.loc[combo, "cbca_pvalue_rel_ATRA"] = stats.ttest_ind(
            combo_cbca_reps, atra_cbca, equal_var=False
        ).pvalue
        hits_cyto.loc[combo, "nl_pvalue_rel_ATRA"] = stats.ttest_ind(
            combo_nl_reps, atra_nl, equal_var=False
        ).pvalue

    hits_cyto = hits_cyto[hits_cyto["cbca_pvalue_rel_ATRA"] < 0.05]
    hits_cyto = hits_cyto[HSAHITS_CYTO_COLS]

    return hits_p05, hits_cyto


def verify_pair(ref_path, fresh_path):
    if not ref_path.exists():
        print(f"  [SKIP] no reference: {ref_path.relative_to(ROOT)}")
        return True
    if filecmp.cmp(ref_path, fresh_path, shallow=False):
        print(f"  [OK ] byte-identical: {fresh_path.name}")
        return True
    a = pd.read_csv(ref_path, index_col=0)
    b = pd.read_csv(fresh_path, index_col=0)
    if list(a.columns) != list(b.columns):
        print(f"  [FAIL] column mismatch: {fresh_path.name}")
        print(f"         ref:   {list(a.columns)}")
        print(f"         fresh: {list(b.columns)}")
        return False
    common_idx = a.index.intersection(b.index)
    extra_ref = sorted(set(a.index) - set(b.index))
    extra_new = sorted(set(b.index) - set(a.index))
    if extra_ref or extra_new:
        print(f"  [FAIL] index mismatch: {fresh_path.name}")
        if extra_ref:
            print(f"         only in ref: {extra_ref}")
        if extra_new:
            print(f"         only in fresh: {extra_new}")
        return False
    a = a.loc[common_idx].select_dtypes(include="number")
    b = b.loc[common_idx, a.columns]
    max_diff = float((a - b).abs().max().max())
    ok = max_diff < 1e-6
    tag = "OK " if ok else "FAIL"
    print(f"  [{tag}] {fresh_path.name} — max numeric |new - old| = {max_diff:.2e}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--parquet",
        type=Path,
        default=PARQUET,
        help=f"Path to complete.parquet (default: {PARQUET.relative_to(ROOT)})",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=SCREEN,
        help="Where the HSAhits_*.csv files go (default: data/screen/)",
    )
    ap.add_argument(
        "--verify",
        action="store_true",
        help="After writing, diff against the existing CSVs and assert match.",
    )
    args = ap.parse_args()

    if not args.parquet.exists():
        raise SystemExit(
            f"missing {args.parquet} — run the plate-loading step first"
        )

    out_p05 = args.out_dir / "HSAhits_p05.csv"
    out_cyto = args.out_dir / "superhits.csv"

    references = []
    if args.verify and args.out_dir == SCREEN:
        snap = Path("/tmp") / "screen_step3_ref"
        snap.mkdir(parents=True, exist_ok=True)
        for src, name in [(out_p05, "p05"), (out_cyto, "cyto")]:
            if src.exists():
                dst = snap / f"HSAhits_{name}.csv"
                shutil.copy2(src, dst)
                references.append(dst)
            else:
                references.append(src)

    print(f"Computing superhits from {args.parquet.relative_to(ROOT)} ...")
    complete = load_complete(args.parquet)
    hits_p05, hits_cyto = compute(complete)
    hits_p05.to_csv(out_p05)
    hits_cyto.to_csv(out_cyto)
    print(f"  wrote {out_p05.relative_to(ROOT)} ({len(hits_p05)} hits)")
    print(f"  wrote {out_cyto.relative_to(ROOT)} ({len(hits_cyto)} hits)")

    if args.verify:
        print("\nVerifying ...")
        ok1 = verify_pair(references[0], out_p05)
        ok2 = verify_pair(references[1], out_cyto)
        if not (ok1 and ok2):
            raise SystemExit("VERIFY FAILED — split would alter downstream data.")
        print("Step 3 verified.")


if __name__ == "__main__":
    main()
