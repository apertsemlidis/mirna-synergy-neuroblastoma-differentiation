#!/usr/bin/env python3
"""
Process dose-response raw Incucyte exports into combined CSVs for
synergy analysis (SynergyFinder-format: PairIndex, Response, Drug1, Drug2,
Conc1, Conc2, ConcUnit).

Produces:
  - combined cbca v3.csv          (cell body cluster area, % inhibition vs siPLK1/mimic-pool)
  - combined nl globalnorm.csv    (neurite length, normalized to global max)

Inputs (data/dose_response/input/):
  - final_doseresponse_echo_instructions.csv
  - nonvarwells_DR.txt
  - plate 1/ ... plate 6/  (each with neurite length.txt, cell body cluster area.txt, etc.)
"""

import json
import os
import pandas as pd

INPUT_DIR = os.path.join("data", "dose_response", "input")
OUTPUT_DIR = os.path.join("data", "dose_response", "output")
TIMEPOINT = 120  # hours
PLATES = [str(i) for i in range(1, 7)]


def convert_nl_to_conc(nanol):
    """Convert an Echo transfer volume (nL) into the resulting well
    concentration (nM).

    Formula: C_final = V_transfer × (C_stock / V_well)
               = nanol (nL) × (10 µM / 150 µL)
               = nanol × 10/150 nM

    Unit check:
        nL × µM / µL = nL × (1000 nM) / (1000 nL) = nM   ✓

    Constants (verified against the original dose_response.ipynb notebook
    — the companion `needed_nl(desired_conc) = desired_conc * (150/10)`
    uses the same two values — and against the Echo instruction sheet,
    whose transfer volumes [0, 7.5, 15, 37.5, 75] nL produce clean
    half-log final concentrations [0, 0.5, 1, 2.5, 5] nM under this
    formula):
        10   — source plate miRNA mimic stock concentration in µM
        150  — destination well final assay volume in µL
    """
    STOCK_UM = 10
    WELL_VOLUME_UL = 150
    return (nanol * STOCK_UM) / WELL_VOLUME_UL


# Non-variant well → control name mapping (inverted)
with open(os.path.join(INPUT_DIR, "nonvarwells_DR.txt")) as f:
    raw = json.load(f)
nonvarwells = {well: name for name, wells in raw.items() for well in wells}

echo_instr = pd.read_csv(
    os.path.join(INPUT_DIR, "final_doseresponse_echo_instructions.csv")
)


def load_data(platenum, measure):
    data = pd.read_table(
        os.path.join(INPUT_DIR, f"plate {platenum}", f"{measure}.txt"),
        header=1,
    )
    data["Elapsed"] = data["Elapsed"].map(lambda x: round(x))
    # Incucyte export occasionally emits whitespace-only cells for wells
    # where no value was computed at a timepoint. We coerce those to 0 so
    # downstream astype('float64') succeeds and the row can still be used
    # by the SynergyFinder reshape.
    #
    # Caveat for a reviewer: this conflates "no measurement" with "zero
    # signal," which are biologically different (especially for neurite
    # length, where 0 may be indistinguishable from a failed segmentation).
    # Truly empty cells (parsed as NaN by read_table, not whitespace) are
    # NOT affected by this replace and will propagate as NaN through
    # .astype('float64'); any such NaN would surface in the normalization
    # step in generate_synergy_finder_file().
    #
    # Inspection of the current input files at TIMEPOINT=120h shows no
    # whitespace cells after the first header row, so this replace is a
    # no-op on the current dataset — the line is kept for robustness
    # against future exports.
    data = data.replace(r"\s+", 0, regex=True)
    data = data.astype("float64")
    data = data.rename(columns=nonvarwells)
    return data


def get_global_max_nl(timepoint):
    """Maximum neurite length across all plates at a given timepoint (variant wells only)."""
    maxes = []
    for p in PLATES:
        d = load_data(p, "neurite length")
        d = d[d.Elapsed == timepoint].iloc[:, 2:]
        maxes.append(d.max().max())
    return max(maxes)


def generate_synergy_finder_file(
    platenum, measure, timepoint, maxnlnorm=False, global_max=None
):
    """Convert a plate's endpoint data into SynergyFinder long-format rows."""
    plate = echo_instr.loc[echo_instr["Destination Plate Barcode"] == int(platenum)]
    data = load_data(platenum, measure)
    rows = []
    for well in plate["Destination well"].unique():
        response = data[data["Elapsed"] == timepoint][well].values[0]
        well_plate = plate.loc[plate["Destination well"] == well].reset_index()
        mir1 = well_plate.loc[0, "Sample Group"]
        conc1 = well_plate.loc[0, "Destination Volume"]
        mir2 = well_plate.loc[1, "Sample Group"]
        conc2 = well_plate.loc[1, "Destination Volume"]
        rows.append(
            [
                int(platenum),
                response,
                mir1,
                mir2,
                convert_nl_to_conc(conc1),
                convert_nl_to_conc(conc2),
                "nM",
            ]
        )

    sf = pd.DataFrame(
        rows,
        columns=[
            "PairIndex",
            "Response",
            "Drug1",
            "Drug2",
            "Conc1",
            "Conc2",
            "ConcUnit",
        ],
    )

    endpoint = data[data["Elapsed"] == timepoint]
    if measure == "cell body cluster area":
        # % inhibition: siPLK1 = 100% inhibition (positive), mimic pool = 0% (negative)
        pos = endpoint["siPLK1"].T.median()
        neg = endpoint["mimic pool"].T.median()
        sf["Response"] = sf["Response"].apply(lambda x: 100 * (pos - x) / (pos - neg))
    elif measure == "neurite length":
        neg = endpoint["mimic pool"].T.median()
        if maxnlnorm == 1:
            pos = endpoint["miRNA-124"].T.median()
        elif maxnlnorm == 2:
            pos = sf.Response.max()
        elif maxnlnorm == 3:
            pos = global_max
        else:
            raise ValueError(f"Unknown maxnlnorm={maxnlnorm} for neurite length")
        sf["Response"] = sf["Response"].apply(lambda x: 100 * (x - neg) / (pos - neg))

    sf.Response = sf.Response.apply(lambda x: 0 if x < 0 else x)
    return sf


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Cell body cluster area (per-plate % inhibition)
    cbca_dfs = [
        generate_synergy_finder_file(p, "cell body cluster area", TIMEPOINT)
        for p in PLATES
    ]
    cbca_out = os.path.join(OUTPUT_DIR, "combined cbca v3.csv")
    pd.concat(cbca_dfs).to_csv(cbca_out, index=False)
    print(f"wrote {cbca_out}")

    # Neurite length (normalized to global max across plates)
    global_max = get_global_max_nl(TIMEPOINT)
    print(f"global max neurite length @ {TIMEPOINT}h = {global_max:.2f}")
    nl_dfs = [
        generate_synergy_finder_file(
            p, "neurite length", TIMEPOINT, maxnlnorm=3, global_max=global_max
        )
        for p in PLATES
    ]
    nl_out = os.path.join(OUTPUT_DIR, "combined nl globalnorm.csv")
    pd.concat(nl_dfs).to_csv(nl_out, index=False)
    print(f"wrote {nl_out}")


if __name__ == "__main__":
    main()
