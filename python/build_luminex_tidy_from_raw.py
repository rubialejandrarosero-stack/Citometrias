#!/usr/bin/env python3
"""Rebuild output/tidy/luminex_tidy.csv from the two raw Luminex xPONENT export CSVs
(26-plex + 4-plex panels), replacing whatever the previous pipeline produced, so all
downstream Luminex figures (PCA, heatmap, key-cytokines panel) can be regenerated from
a verified source of truth.

Plate layout / sample naming (confirmed identical in both raw files, well positions match):
  CI1/CI2/CI3   = Control Inflammatory, timepoint 1/2/3 (24/48/96h), 4 technical well replicates each
  CN1/CN2/CN3   = Control Basal (non-inflammatory), same pattern
  I.D{1-4}.{t}  = Inflammatory, Resting PBMC (no asterisk), donor D1-4, time t in {24,48,96}
  I*.D{1-4}.{t} = Inflammatory, Activated PBMC (asterisk), donor, time
  N.D{1-4}.{t}  = Basal, Resting PBMC, donor, time
  N*.D{1-4}.{t} = Basal, Activated PBMC, donor, time
Excluded (not part of the experimental design): CMG, CM24, MG.1, MG.2, M24.1, M24.2
Uses the "Result" DataType block (standard-curve-interpolated concentration, pg/mL) -- NOT
"Median" (raw bead fluorescence intensity), which the first version of this script incorrectly
used. Verified: for the 78 real (non-Standard/Background) sample rows, "Result" is 100% clean
numeric data (no Invalid/OOR/NaN) -- those flags only appear on Standard curve-fitting rows.

Output schema matches the existing output/tidy/luminex_tidy.csv:
  cytokine, environment, pbmc, time, donor, value
  pbmc in {None, Resting, Activated}; donor is empty for Control (None), 1-4 otherwise.
"""
import csv
import re
import pandas as pd

RAW_FILES = [
    "/mnt/c/Users/57319/Downloads/RUBI 02.08.2024_26CK.csv",
    "/mnt/c/Users/57319/Downloads/RUBY 090224_4CK.csv",
]
EXCLUDE = {"CMG", "CM24", "MG.1", "MG.2", "M24.1", "M24.2"}
TIME_OF_CODE = {"1": 24, "2": 48, "3": 96}


def parse_sample(sample):
    """Return (environment, pbmc, time, donor) or None if not part of the design."""
    if sample in EXCLUDE or sample.startswith("Standard") or sample.startswith("Background"):
        return None

    m = re.match(r"^C([IN])(\d)$", sample)
    if m:
        env = "Inflammatory" if m.group(1) == "I" else "Basal"
        return env, "None", TIME_OF_CODE[m.group(2)], None

    m = re.match(r"^([IN])(\*?)\.D(\d)\.(\d+)$", sample)
    if m:
        env = "Inflammatory" if m.group(1) == "I" else "Basal"
        pbmc = "Activated" if m.group(2) == "*" else "Resting"
        donor = int(m.group(3))
        time = int(m.group(4))
        return env, pbmc, time, donor

    return None  # unrecognised sample name -> skip, don't guess


def extract_result_block(path):
    with open(path, newline="", encoding="latin-1") as f:
        rows = list(csv.reader(f))
    start = next(i for i, r in enumerate(rows) if r[:1] == ["DataType:"] and r[1:2] == ["Result"])
    header = rows[start + 1]
    analytes = header[2:-1]  # drop Location, Sample, Total Events
    records = []
    for row in rows[start + 2:]:
        if not row or not row[0].strip():
            break  # blank line ends the block
        if row[0] == "Location":
            break
        sample = row[1]
        parsed = parse_sample(sample)
        if parsed is None:
            continue
        environment, pbmc, time, donor = parsed
        for analyte, val in zip(analytes, row[2:2 + len(analytes)]):
            if val == "" or val is None:
                continue
            v = val.strip()
            # right/left-censored values from xPONENT ("> 960" = saturated above the standard
            # curve's top point, "< 5" = below its bottom point) -> use the limit itself as a
            # conservative point estimate, per standard Luminex practice. Genuinely missing/
            # unusable points ("Invalid(...)", "NaN") are dropped, not guessed.
            if v.startswith(">") or v.startswith("<"):
                v = v[1:].strip()
            try:
                fval = float(v)
            except ValueError:
                continue  # "Invalid(...)", "NaN", etc. -> skip, genuinely unusable
            records.append({"cytokine": analyte, "environment": environment, "pbmc": pbmc,
                            "time": time, "donor": donor, "sample": sample, "value": fval})
    return pd.DataFrame(records)


def main():
    frames = [extract_result_block(p) for p in RAW_FILES]
    df = pd.concat(frames, ignore_index=True)

    # Control (pbmc == None) has 4 technical well replicates per cytokine/env/time -> average
    # them to a single value, matching the previous tidy file's convention (4 repeated rows were
    # actually the 4 replicate values, not 4 donors) -- verify replicate count then collapse.
    ctrl = df[df.pbmc == "None"]
    rep_counts = ctrl.groupby(["cytokine", "environment", "time"]).size()
    print("Control replicate counts (should be 4 each):")
    print(rep_counts.describe())

    pbmc_df = df[df.pbmc != "None"].drop(columns=["sample"])
    ctrl_df = (ctrl.groupby(["cytokine", "environment", "pbmc", "time"], as_index=False)
               .agg(value=("value", "mean")))
    ctrl_df["donor"] = pd.NA

    out = pd.concat([pbmc_df, ctrl_df[pbmc_df.columns]], ignore_index=True)
    out = out.sort_values(["cytokine", "environment", "pbmc", "time", "donor"]).reset_index(drop=True)

    print("\nTotal analytes:", out.cytokine.nunique())
    print("Analytes:", sorted(out.cytokine.unique()))
    print("\nRows:", len(out))
    print(out.head(10))

    out.to_csv("output/tidy/luminex_tidy.csv", index=False)
    print("\nsaved output/tidy/luminex_tidy.csv")


if __name__ == "__main__":
    main()
