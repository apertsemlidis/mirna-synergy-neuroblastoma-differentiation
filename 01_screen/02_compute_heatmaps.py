#!/usr/bin/env python3
"""
Step 2 of the screen_analysis_v14.py split
(monolith now at scripts/archive/screen_outdated/screen_analysis_v14.py).

Reads data/screen/complete.parquet (produced by load_screen.py) and
writes the 44x44 HSA and ABS heatmap tables for the canonical screen
slice — the same files the monolith produced via lines 542-556 of the
archived screen_analysis_v14.py:

    for m in ['neurite length', 'cell body cluster area']:
        for t in [slice(96,126)]:
            d = HSA_df(m, t, complete)
            d.to_csv(f'HSA_dfs/{initials(m)}_{t}.csv')
            d = abs_heatmap_df(m, t, 'ATRA (25 uM)', complete)
            d.to_csv(f'ABS_dfs/{initials(m)}_{t}.csv')

Outputs:
    data/screen/HSA_dfs/{nl,cbca}_slice(96, 126, None).csv
    data/screen/ABS_dfs/{nl,cbca}_slice(96, 126, None).csv

Run `--verify` to byte-compare the freshly written CSVs against the
existing ones (the ones currently consumed by
screen_heatmaps_composite_v17 and the volcano stats pipeline). If the
parquet path matches the monolith path, the diff should be zero.

Functions are imported from load_screen.py so step 1 stays the source of
truth for `complete` and the HSA/ABS math.
"""

import argparse
import filecmp
import shutil
from pathlib import Path

import pandas as pd

from screen_helpers import abs_heatmap_df, hsa_df, load_complete

ROOT = Path(__file__).resolve().parents[2]
SCREEN = ROOT / "data" / "screen"
PARQUET = SCREEN / "complete.parquet"

# (file_initials, metric_name) — matches the screen_analysis_v14 convention.
METRICS = [
    ("nl", "neurite length"),
    ("cbca", "cell body cluster area"),
]
TIME_SLICE = slice(96, 126)


def slice_tag(s: slice) -> str:
    """Reproduce the monolith's filename suffix: str(slice(96,126)) →
    'slice(96, 126, None)'."""
    return str(s)


def write_tables(parquet_path: Path, out_dir: Path):
    complete = load_complete(parquet_path)
    hsa_dir = out_dir / "HSA_dfs"
    abs_dir = out_dir / "ABS_dfs"
    hsa_dir.mkdir(parents=True, exist_ok=True)
    abs_dir.mkdir(parents=True, exist_ok=True)

    written = []
    tag = slice_tag(TIME_SLICE)
    for ini, metric in METRICS:
        hsa = hsa_df(complete, metric, TIME_SLICE)
        abs_ = abs_heatmap_df(complete, metric, TIME_SLICE)

        hsa_path = hsa_dir / f"{ini}_{tag}.csv"
        abs_path = abs_dir / f"{ini}_{tag}.csv"
        hsa.to_csv(hsa_path)
        abs_.to_csv(abs_path)
        written.extend([hsa_path, abs_path])
        print(f"  wrote {hsa_path.relative_to(ROOT)}")
        print(f"  wrote {abs_path.relative_to(ROOT)}")
    return written


def verify_against(reference_paths, fresh_paths):
    """Compare each freshly written file against an existing reference.
    Reports byte-identity and falls back to numeric diff."""
    all_ok = True
    for ref, fresh in zip(reference_paths, fresh_paths):
        if not ref.exists():
            print(f"  [SKIP] no reference: {ref.relative_to(ROOT)}")
            continue
        if filecmp.cmp(ref, fresh, shallow=False):
            print(f"  [OK ] byte-identical: {fresh.name}")
            continue
        # Bytes differ — check numerically.
        a = pd.read_csv(ref, index_col=0).astype(float)
        b = pd.read_csv(fresh, index_col=0).astype(float)
        b = b.loc[a.index, a.columns]
        max_diff = float((a - b).abs().max().max())
        ok = max_diff < 1e-4
        all_ok &= ok
        tag = "OK " if ok else "FAIL"
        print(f"  [{tag}] {fresh.name} — max numeric |new - old| = {max_diff:.2e}")
    return all_ok


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
        help="Where {HSA,ABS}_dfs/ should be written (default: data/screen/)",
    )
    ap.add_argument(
        "--verify",
        action="store_true",
        help="After writing, diff against the pre-existing CSVs and assert match.",
    )
    args = ap.parse_args()

    if not args.parquet.exists():
        raise SystemExit(
            f"missing {args.parquet} — run scripts/screen/load_screen.py first"
        )

    # Snapshot the existing CSVs before we overwrite, so --verify has
    # something to compare against even if we're writing in-place.
    tag = slice_tag(TIME_SLICE)
    reference_paths = []
    if args.verify and args.out_dir == SCREEN:
        snap_dir = Path("/tmp") / "screen_step2_ref"
        snap_dir.mkdir(parents=True, exist_ok=True)
        for ini, _ in METRICS:
            for kind in ("HSA", "ABS"):
                src = SCREEN / f"{kind}_dfs" / f"{ini}_{tag}.csv"
                if src.exists():
                    dst = snap_dir / f"{kind}_{ini}_{tag}.csv"
                    shutil.copy2(src, dst)
                    reference_paths.append(dst)
                else:
                    reference_paths.append(src)  # will trigger SKIP

    print(f"Writing HSA/ABS tables from {args.parquet.relative_to(ROOT)} ...")
    written = write_tables(args.parquet, args.out_dir)

    if args.verify:
        if not reference_paths:
            # Out-of-tree dir: compare new files against the in-tree canonical ones.
            reference_paths = [
                SCREEN / "HSA_dfs" / f"nl_{tag}.csv",
                SCREEN / "ABS_dfs" / f"nl_{tag}.csv",
                SCREEN / "HSA_dfs" / f"cbca_{tag}.csv",
                SCREEN / "ABS_dfs" / f"cbca_{tag}.csv",
            ]
        print("\nVerifying ...")
        if not verify_against(reference_paths, written):
            raise SystemExit("VERIFY FAILED — split would alter downstream data.")
        print("Step 2 verified.")


if __name__ == "__main__":
    main()
