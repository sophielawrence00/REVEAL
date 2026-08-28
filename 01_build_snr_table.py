"""
Averages replicate SNR columns per lipid/condition from the SNR Excel sheet,
expands abbreviated condition names to match REVEAL's naming, and writes
a tidy CSV: mz_key, protein, condition, snr_mean
"""
import os
import pandas as pd
import re

# Paths are resolved relative to this script, so it can be run from any directory.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "Example Data")   # <-- point this at your data folder

SNR_XLSX = os.path.join(DATA_DIR, "snr_data.xlsx")
OUT_CSV  = os.path.join(DATA_DIR, "snr_long.csv")

df = pd.read_excel(SNR_XLSX, sheet_name=0)
df = df.rename(columns={df.columns[0]: "mz_key"})

# Parse each replicate column name into (protein, condition, rep)
# Expects full condition names already, e.g. mGlyR_Glutamate_1, mGlyR_1 (apo),
# AmtB_BrainLipids_1. Any protein name is accepted, not a fixed list.
col_pattern = re.compile(r"^(?P<protein>[A-Za-z0-9]+)_(?:(?P<cond>[A-Za-z0-9]+)_)?(?P<rep>\d+)$")

groups = {}  # (protein, condition) -> list of column names
for col in df.columns[1:]:
    m = col_pattern.match(col)
    if not m:
        print(f"  WARNING: could not parse column '{col}', skipping")
        continue
    protein = m.group("protein")
    cond    = m.group("cond") or "apo"
    groups.setdefault((protein, cond), []).append(col)

rows = []
for (protein, cond), cols in groups.items():
    mean_snr = df[cols].mean(axis=1, skipna=True)
    for mz_key, val in zip(df["mz_key"], mean_snr):
        if pd.isna(val):
            continue
        rows.append({"mz_key": mz_key, "protein": protein, "condition": cond, "snr_mean": val})

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)
print(f"Wrote {len(out)} rows to {OUT_CSV}")
print(out.head(15))
