#!/usr/bin/env python3
"""
Step 1 of the screen-processing pipeline.

Loads raw Incucyte plate files, applies well-mapping, cleaning, and
two-point normalization, and saves the resulting
`complete` dict as a single multi-index parquet at:

    data/screen/complete.parquet

Also saves the list of control conditions detected during loading at:

    data/screen/controls.csv

Downstream steps read these
artefacts instead of re-loading 152 MB of raw plates.

Run `--verify` to additionally recompute the 44x44 HSA and ABS heatmap
tables from the saved parquet and assert byte-equivalence (within float
tolerance) against the existing {HSA,ABS}_dfs/*.csv files. This confirms
the data path is preserved.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

from screen_helpers import abs_heatmap_df, hsa_df, load_complete

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "screen_raw" / "Incucyte_raw_data_v2"
ARCHIVE_DIR = ROOT / "scripts" / "archive" / "screen_outdated"
OUT_DIR = ROOT / "data" / "screen"

ECHO_CSV = ARCHIVE_DIR / "44pairs_echo_instructions.csv"
NONVAR_JSON = ARCHIVE_DIR / "nonvarwells.txt"
PLATE4_JSON = ARCHIVE_DIR / "plate4nonvarwells.txt"

PLATE_DIRS = [f"plate {i}-{j}" for i in range(1, 5) for j in range(1, 4)]
ALL_MEASURES = [
    "cell body cluster area",
    "cell body clusters",
    "neurite branch points per body cluster area",
    "neurite branch points per body cluster",
    "neurite branch points",
    "neurite length per body cluster area",
    "neurite length per body cluster",
    "neurite length",
]
NORMALIZE_TIME = 120


# --- helpers kept inline so this file has no relative-path import
#     dependency. ---


def mapper(df, plate, instructions, nonvarwells, plate4_fixed):
    clip = instructions.loc[instructions["Destination Plate Barcode"] == int(plate)]
    cols = {}
    for column in df.columns:
        mirs = clip.loc[clip["Destination well"] == column]["Sample Group"]
        if not mirs.empty:
            cols[column] = " + ".join(sorted(set(mirs.values)))
    df = df.rename(columns=cols).rename(columns=nonvarwells)
    if int(plate) == 4:
        df = df.rename(columns=plate4_fixed)
    return df


def normalize_cellmetric(replicates, neg, pos, time):
    for key, df in replicates.items():
        if key[2] == "cell body cluster area":
            p = df.loc[pos].median()[time]
            n = df.loc[neg].median()[time]
            replicates[key] = df.apply(lambda x: (p - x) / (p - n))
    return replicates


def normalize_neuritemetric(replicates, neg, pos, time):
    for key, df in replicates.items():
        if key[2] not in ("cell body cluster area", "cell body clusters"):
            p = df.loc[pos].median()[time]
            n = df.loc[neg].median()[time]
            replicates[key] = df.apply(lambda x: (x - n) / (p - n))
    return replicates


def load_replicates(measures):
    instructions = pd.read_csv(ECHO_CSV)
    with open(NONVAR_JSON) as f:
        nonvarwells = json.load(f)
    with open(PLATE4_JSON) as f:
        plate4_fixed = json.load(f)

    replicates = {}
    for folder in PLATE_DIRS:
        for measure in measures:
            txt = RAW_DIR / folder / f"{measure}.txt"
            df = pd.read_table(txt, header=1)
            df = mapper(df, folder[6], instructions, nonvarwells, plate4_fixed)
            df = df.drop("Date Time", axis=1)

            df = df.T
            df.columns = df.loc["Elapsed"].map(round)
            df = df.drop("Elapsed", axis=0)
            df.index = df.index.str.replace(r"\.[0-9]+", "", regex=True)
            for t in (138, 144):
                if t in df.columns:
                    df = df.drop(t, axis=1)
            df = df.replace(r"\s+", 0, regex=True).dropna().astype("float64")

            key = (folder[6], folder[8], measure)
            replicates[key] = df
    return replicates


def build_complete(replicates, measures):
    """Concat the 12 (plate, rep) frames for each metric. Mirrors the
    original pipeline's `complete[m] = pd.concat([...])` step."""
    complete = {}
    for measure in measures:
        frames = [
            replicates[(str(i), str(j), measure)]
            for i in range(1, 5)
            for j in range(1, 4)
        ]
        complete[measure] = pd.concat(frames)
    return complete


def detect_controls(replicates):
    """Controls have 8 replicate wells per plate; combos have 1.
    Identify them by counting index value occurrences in plate 1 rep 1."""
    rep11 = replicates[("1", "1", "cell body cluster area")]
    counts = rep11.index.value_counts()
    return counts[counts == 8].index.tolist()


def save_complete(complete, out_path):
    """Pack the per-metric dict into a single long DataFrame and write
    parquet. Layout:
        index columns: [metric, condition]
        data columns:  time points (one column per hour)
        values:        float64 measurement
    """
    frames = []
    for metric, df in complete.items():
        frame = df.copy()
        # Cast columns (time points) to str so they round-trip cleanly.
        frame.columns = [str(c) for c in frame.columns]
        frame.insert(0, "condition", frame.index)
        frame.insert(0, "metric", metric)
        frames.append(frame.reset_index(drop=True))
    combined = pd.concat(frames, ignore_index=True)
    combined.to_parquet(out_path, index=False)


def verify(parquet_path):
    """Recompute HSA + ABS tables from parquet and diff against the
    existing data/screen/{HSA,ABS}_dfs/*.csv produced by the original pipeline."""
    complete = load_complete(parquet_path)
    time_slice = slice(96, 126)
    pairs = {"nl": "neurite length", "cbca": "cell body cluster area"}

    all_ok = True
    for ini, metric in pairs.items():
        for kind, fn in (("HSA", hsa_df), ("ABS", abs_heatmap_df)):
            new = fn(complete, metric, time_slice)
            csv = OUT_DIR / f"{kind}_dfs" / f"{ini}_slice(96, 126, None).csv"
            old = pd.read_csv(csv, index_col=0)
            # Align on the ordered index of the new table.
            old = old.loc[new.index, new.columns]
            diff = (new.astype(float) - old.astype(float)).abs()
            max_diff = float(diff.max().max())
            ok = max_diff < 1e-4
            all_ok &= ok
            tag = "OK " if ok else "FAIL"
            print(f"  [{tag}] {kind:3s} {ini:4s}  max |new - old| = {max_diff:.2e}")
    return all_ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--measures",
        nargs="*",
        default=ALL_MEASURES,
        help="Subset of metrics to load (default: all 8).",
    )
    ap.add_argument(
        "--verify",
        action="store_true",
        help="After saving, recompute HSA/ABS tables and diff against existing CSVs.",
    )
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_parquet = OUT_DIR / "complete.parquet"
    out_controls = OUT_DIR / "controls.csv"

    print(
        f"Loading {len(args.measures)} metric(s) × {len(PLATE_DIRS)} plate/rep frames from "
        f"{RAW_DIR.relative_to(ROOT)} ..."
    )
    replicates = load_replicates(args.measures)
    print(f"  loaded {len(replicates)} (plate, rep, metric) frames")

    replicates = normalize_cellmetric(
        replicates, "miRNA mimic pool (10 nM)", "siPLK1 (10 nM)", NORMALIZE_TIME
    )
    replicates = normalize_neuritemetric(
        replicates, "miRNA mimic pool (10 nM)", "ATRA (25 uM)", NORMALIZE_TIME
    )
    print(f"  normalized at t={NORMALIZE_TIME}h")

    complete = build_complete(replicates, args.measures)
    save_complete(complete, out_parquet)
    print(
        f"  saved {out_parquet.relative_to(ROOT)} "
        f"({out_parquet.stat().st_size / 1024:.1f} KiB)"
    )

    controls = detect_controls(replicates)
    pd.Series(controls, name="condition").to_csv(out_controls, index=False)
    print(f"  saved {out_controls.relative_to(ROOT)} ({len(controls)} controls)")

    if args.verify:
        print("\nVerifying against existing HSA/ABS CSVs ...")
        if not verify(out_parquet):
            raise SystemExit("VERIFY FAILED — split would alter downstream data.")
        print("All HSA/ABS tables match within tolerance. Step 1 verified.")


if __name__ == "__main__":
    main()
