# REVEAL lipid binding correlation classification

Three python scripts that help analyse the data from the REVEAL wrokflow. These scripts correlate each lipid's reconstructed MS2 precursor spectrum (following application of REVEAL) against the MS1 envelope and cross-check the result against signal-to-noise data. All spectral input data (per-lipid MS2 reconstructions and MS1 envelopes) is unnornmalised and stored as individual tab-separated text files, one per sample/lipid.

## Pipeline overview
Out of the box all three scripts point at the bundled `Example Data/` folder, so they can be
run as-is to reproduce the AmtB example. Paths are resolved relative to each script, so the
scripts can be run from any working directory. To use your own data, edit the `DATA_DIR` /
`ROOT_DIR` variable at the top of each script.

Run the scripts in filename order:

```
01_build_snr_table.py
        │  snr_data.xlsx  →  snr_long.csv
        ▼
02_pearson_vs_ms1.py
        │  raw lipid .txt files + MS1 .txt files + Master.csv
        ▼  →  {protein}_pearson_results.csv, heatmaps, r-distributions
03_classify_binders.py
           snr_long.csv + {protein}_pearson_results.csv
           →  r-vs-SNR scatter plots, binder count charts, confident_binders.xlsx
```
## `01_build_snr_table.py`

Averages replicate signal-to-noise ratio (SNR) columns per lipid/condition
from an Excel sheet and writes a tidy long-format CSV. For each lipid's reconstructed precursor
spectrum, the SNR is calculated by taking the unnormalised area under
the curve in the m/z region where the expected MS1 signal should fall, and comparing it to the area under the curve in an
m/z region where no MS1 signal is expected. This per-replicate SNR value is calculated manually
(not by any script in this repo) and provided as input in `snr_data.xlsx`.
`01_build_snr_table.py` then averages the replicate SNR values for each
lipid/condition to produce one mean SNR per lipid/condition, which is used
downstream as a data-quality filter alongside the Pearson r classification.

**Input:** `snr_data.xlsx` - first column is the lipid `mz_key`; remaining
columns are named `{protein}_{condition}_{rep}` (e.g. `mGlyR_Glutamate_1`)
or `{protein}_{rep}` for the apo condition (e.g. `mGlyR_1`).

**Output:** `snr_long.csv` with columns `mz_key, protein, condition, snr_mean`.

## `02_pearson_vs_ms1.py`

Auto-discovers protein/condition/replicate experiments
from subfolder names, averages each lipid's replicate reconstructed MS2
precursor spectra, correlates that reconstructed spectrum against the MS1
envelope (Pearson r), derives a significance threshold via permutation
testing, and produces per-protein heatmaps and results tables.

**Inputs:** All data must be converted into a `.txt`. Place all reconstructed MS2 precursor spectra as tab-separated `mz, intensity` pairs into subfolders for each condition/replicate. The `.txt` filename is the detected mz and this stem (minus `.txt`) becomes the
lipid's `mz_key`. 

The MS1 envelope for each protein/condition is likewise
stored as its own individual `.txt` file (tab-separated `mz, intensity`). 
All intensitites are raw intensity.

`Master.csv` is a `.csv` file where the first column is the lipid mz (same as the `.txt` file names), the second column is the assigned lipid and the third column is any further possibly lipid identities.

**Example folder structure** (auto-detected under `ROOT_DIR`):

```
ROOT_DIR/
├── mGlyR_1/                    ← protein=mGlyR, condition=apo, rep=1
├── mGlyR_2/
├── mGlyR_3/
├── mGlyR_Glycine_1/            ← protein=mGlyR, condition=Glycine, rep=1
├── mGluR2_1/
├── mGluR2_Glut_1/
├── mGlyR_MS1.txt                ← MS1 envelope for mGlyR apo
├── mGlyR_Glycine_MS1.txt        
├── mGluR2_Glutamate_MS1.txt     
└── Master.csv                   ← mz_key → lipid name(s) lookup
```

**Outputs:** (per protein, into `OUT_DIR`)

   * `{protein}_pearson_heatmap.pdf` - r-value heatmap with MS1 intensity strip
   * `{protein}_binder_counts.pdf` - bar chart of binders per condition
   * `{protein}_{condition}_r_distribution.pdf` - histogram with threshold line
   * `{protein}_pearson_results.csv` - full results table (wide format, one row per lipid)

## `03_classify_binders.py`

Merges the Pearson r results from script 02 with the SNR table from
script 01, and classifies "confident binders" as lipids with
**r > threshold AND SNR > 3**. 

**Inputs:**
* `{protein}_pearson_results.csv` files (from script 02)
* `snr_long.csv` (from script 01)

**Output**, written to `RESULTS_DIR`:
* `{protein}_r_vs_snr_scatter.pdf` - scatter plot of r vs SNR, coloured by condition, with threshold lines
* `{protein}_binder_counts_dual.pdf` — bar chart of number of "confident binders" per condition. Separate bar charts for each protein.
* `mGlyR_vs_mGluR2_apo_scatter.pdf` / `mGlyR_vs_mGluR2_apo_binder_counts.pdf` - cross-protein apo comparison (only if both proteins have apo data)
* `confident_binders.xlsx` - an excel sheet per protein listing all confident binders

## Notes

* All plots are saved as high-resolution (300 dpi) PDFs.
* `02_pearson_vs_ms1.py` uses a fixed random seed (`seed=42`) for the permutation test, so results are reproducible between runs on the same data.
* Lipids with an ambiguous name-to-mz mapping (more than one candidate name in `Master.csv`) are flagged with a trailing `*` in heatmap labels.

