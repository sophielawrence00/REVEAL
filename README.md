# REVEAL Lipid Binder Classification Pipeline

Three scripts that together identify which lipids bind a given membrane
protein, by correlating each lipid's reconstructed MS2 precursor spectrum following application of REVEAL against the MS1 envelope and cross-checking the result against signal-to-noise data. All spectral input data (per-lipid MS2 reconstructions and MS1 envelopes) is stored as individual tab-separated text files, one per sample/lipid.

## Pipeline overview

Run the scripts in this order — each one consumes the output of the last:

```
1. build\_snr\_table.py
        │  snr\_data.xlsx  →  snr\_long.csv
        ▼
2. REVEAL\_pearson\_vs\_MS1\_permutateAllSamples.py
        │  raw lipid .txt files + MS1 .txt files + Master.csv
        │  →  {protein}\_pearson\_results.csv, heatmaps, r-distributions
        ▼
3. make\_scatter.py
           snr\_long.csv + {protein}\_pearson\_results.csv
           →  r-vs-SNR scatter plots, binder count charts, confident\_binders.xlsx
```

\---

## 1\. `build\_snr\_table.py`

Averages replicate SNR (signal-to-noise ratio) columns per lipid/condition
from an Excel sheet and writes a tidy long-format CSV.

**What "SNR" means here:** for each lipid's reconstructed precursor
spectrum, the signal-to-noise ratio is calculated by taking the area under
the curve in the m/z region where the expected MS1 signal should fall
("signal" — the m/z space where intensity is seen for the protein's charge
state distribution), and comparing it to the area under the curve in an
m/z region where no MS1 signal is expected ("noise" — m/z space with no
protein charge state distribution present). All intensities (in the
per-lipid text files, the MS1 files, and the SNR calculation) are raw,
unnormalized values. This per-replicate SNR value is calculated manually
(not by any script in this repo) and provided as input in `snr\_data.xlsx`;
`build\_snr\_table.py` then averages the replicate SNR values for each
lipid/condition to produce one mean SNR per lipid/condition, which is used
downstream as a data-quality filter alongside the Pearson r classification.

**Input:** `snr\_data.xlsx` — first column is the lipid `mz\_key`; remaining
columns are named `{protein}\_{condition}\_{rep}` (e.g. `mGlyR\_Glutamate\_1`)
or `{protein}\_{rep}` for the apo condition (e.g. `mGlyR\_1`).

**Output:** `snr\_long.csv` with columns `mz\_key, protein, condition, snr\_mean`.

**Config (edit at top of file):**

* `SNR\_XLSX` — path to the input Excel file
* `OUT\_CSV` — path for the output CSV

\---

## 2\. `REVEAL\_pearson\_vs\_MS1\_permutateAllSamples.py`

The core analysis. Auto-discovers protein/condition/replicate experiments
from subfolder names, averages each lipid's replicate reconstructed MS2
precursor spectra, correlates that reconstructed spectrum against the MS1
envelope (Pearson r), derives a significance threshold via permutation
testing, and produces per-protein heatmaps and results tables.

**Expected folder structure** (auto-detected under `ROOT\_DIR`):

```
ROOT\_DIR/
├── mGlyR\_1/                    ← protein=mGlyR, condition=apo, rep=1
├── mGlyR\_2/
├── mGlyR\_3/
├── mGlyR\_Glycine\_1/            ← protein=mGlyR, condition=Glycine, rep=1
├── mGluR2\_1/
├── mGluR2\_Glut\_1/              ← condition folder name differs from MS1 name (mapped)
├── AmtB\_BrainLipids\_1/
├── mGlyR\_MS1.txt                ← MS1 envelope for mGlyR apo
├── mGlyR\_Glycine\_MS1.txt        ← MS1 envelope for mGlyR Glycine
├── mGluR2\_Glutamate\_MS1.txt     ← MS1 envelope for mGluR2 Glut (name expanded)
└── Master.csv                   ← mz\_key → lipid name(s) lookup
```

Each replicate subfolder should contain one `.txt` file per lipid, holding
that lipid's reconstructed MS2 precursor spectrum as tab-separated
`mz, intensity` pairs; the filename stem (minus `.txt`) becomes the
lipid's `mz\_key`. The MS1 envelope for each protein/condition is likewise
stored as its own individual `.txt` file (tab-separated `mz, intensity`).

**What it does:**

1. Loads `Master.csv` to map `mz\_key` → lipid name(s).
2. Discovers all `{protein}\_{condition}\_{rep}` subfolders.
3. Resolves the matching MS1 file for each protein/condition (applying
`CONDITION\_TO\_MS1` to expand abbreviated names, e.g. `Glut` → `Glutamate`).
4. Averages each lipid's replicate reconstructed MS2 precursor spectra into one matrix per condition.
5. Computes Pearson r of each lipid's reconstructed spectrum against the MS1 envelope.
6. Builds a null distribution via permutation (shuffling the real MS1
vector) and takes a percentile (default: 99th, i.e. p < 0.01) as the
binder threshold — pooled across all conditions for one protein.
7. Outputs, per protein, into `OUT\_DIR`:

   * `{protein}\_pearson\_heatmap.pdf` — r-value heatmap with MS1 intensity strip
   * `{protein}\_binder\_counts.pdf` — bar chart of binders per condition
   * `{protein}\_{condition}\_r\_distribution.pdf` — histogram with threshold line
   * `{protein}\_pearson\_results.csv` — full results table (wide format, one row per lipid)

**Config (edit at top of file):**

* `ROOT\_DIR` — root folder containing the experiment subfolders and MS1 files
* `MASTER\_CSV`, `OUT\_DIR` — derived from `ROOT\_DIR` by default
* `CONDITION\_TO\_MS1` — mapping for condition names that differ between folder names and MS1 filenames
* `PROTEINS\_TO\_PROCESS` — list of proteins to process, or `None` for all discovered
* `N\_REPLICATES`, `N\_PERMUTATIONS`, `PERM\_THRESHOLD` — permutation test parameters
* `THRESHOLD\_MODE` — `"permutation"` or `"hard"`; `HARD\_THRESHOLD` used if `"hard"`
* `MS1\_MODE` — `"own"` (each condition uses its own MS1) or `"apo"` (all conditions share the apo MS1)

\---

## 3\. `make\_scatter.py`

Merges the Pearson r results from script 2 with the SNR table from
script 1, and classifies "confident binders" as lipids with
**r > threshold AND SNR > 3**.

**Input:**

* `{protein}\_pearson\_results.csv` files (from script 2's `OUT\_DIR`)
* `snr\_long.csv` (from script 1)

**Output**, written to `RESULTS\_DIR`:

* `{protein}\_r\_vs\_snr\_scatter.pdf` — scatter of r vs SNR, coloured by condition, with threshold lines
* `{protein}\_binder\_counts\_dual.pdf` — bar chart of confident binder counts per condition
* `mGlyR\_vs\_mGluR2\_apo\_scatter.pdf` / `mGlyR\_vs\_mGluR2\_apo\_binder\_counts.pdf` — cross-protein apo comparison (only if both proteins have apo data)
* `confident\_binders.xlsx` — one styled sheet per protein listing all confident binders

**Config (edit at top of file):**

* `ROOT\_DIR` — must match the `ROOT\_DIR` used in script 2
* `RESULTS\_DIR`, `SNR\_CSV` — derived from `ROOT\_DIR` by default
* `SNR\_CUTOFF` — minimum SNR for a lipid to count as a confident binder (default `3.0`)

\---

## Requirements

```
pandas
numpy
scipy
matplotlib
openpyxl
```

Install with:

```bash
pip install pandas numpy scipy matplotlib openpyxl
```

## Usage

1. Edit the path variables (`SNR\_XLSX`, `ROOT\_DIR`, etc.) at the top of each script to match your data locations.
2. Run in order:

```bash
python build\_snr\_table.py
python REVEAL\_pearson\_vs\_MS1\_permutateAllSamples.py
python make\_scatter.py
```

## Notes

* All plots are saved as high-resolution (300 dpi) PDFs.
* `REVEAL\_pearson\_vs\_MS1\_permutateAllSamples.py` uses a fixed random seed (`seed=42`) for the permutation test, so results are reproducible between runs on the same data.
* Lipids with an ambiguous name-to-mz mapping (more than one candidate name in `Master.csv`) are flagged with a trailing `\*` in heatmap labels.

