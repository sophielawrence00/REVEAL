"""
REVEAL — Lipid Binder Classification by Pearson Correlation vs MS1 Envelope
=============================================================================
Automatically discovers all protein/condition/replicate combinations from
subfolder names, maps each to its corresponding MS1 file, averages replicates,
computes Pearson r vs MS1, applies a permutation null threshold, and outputs
heatmaps + results tables.

FOLDER STRUCTURE (auto-detected):
  ROOT/
  ├── mGlyR_1/                    ← protein=mGlyR, condition=apo, rep=1
  ├── mGlyR_2/
  ├── mGlyR_3/
  ├── mGlyR_Glycine_1/            ← protein=mGlyR, condition=Glycine, rep=1
  ├── mGluR2_1/
  ├── mGluR2_Glut_1/              ← condition folder name != MS1 name (mapped below)
  ├── AmtB_BrainLipids_1/
  ├── mGlyR_MS1.txt               ← MS1 for mGlyR apo
  ├── mGlyR_Glycine_MS1.txt       ← MS1 for mGlyR Glycine
  ├── mGluR2_Glutamate_MS1.txt    ← MS1 for mGluR2 Glut (name expanded)
  └── Master.csv
"""

import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import pearsonr

# ═══════════════════════════════════════════════════════════════════════════
#  USER CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════

# Paths are resolved relative to this script, so it can be run from any directory.
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.join(BASE_DIR, "Example Data")   # <-- point this at your data folder

MASTER_CSV = os.path.join(ROOT_DIR, "Master.csv")
OUT_DIR    = os.path.join(ROOT_DIR, "REVEAL_correlation_output")

# Map subfolder condition string → MS1 filename condition string
# Only needed where they differ. Exact match on the condition part of the folder name.
# e.g. mGluR2_Glut_1 → condition="Glut" → MS1 uses "Glutamate"
CONDITION_TO_MS1 = {
    "Glut"      : "Glutamate",
    # add more here if needed, e.g. "Ala": "Alanine" if they differ
}

# Proteins to process (set to None to process all discovered proteins)
PROTEINS_TO_PROCESS = ["mGlyR", "mGluR2", "AmtB"]   # or None for all

N_REPLICATES   = 3
N_PERMUTATIONS = 1000
PERM_THRESHOLD = 99     # percentile → p < 0.01

THRESHOLD_MODE = "permutation"   # "hard" or "permutation"
HARD_THRESHOLD = 0.5      # used when THRESHOLD_MODE == "hard"

MS1_MODE = "own"   # "own" = each condition uses its own MS1
                   # "apo" = all conditions use the apo MS1

# Aesthetics
CMAP_HEATMAP = "RdBu_r"
COL_BINDER   = "#886dda"

# ═══════════════════════════════════════════════════════════════════════════
#  1. LOAD MASTER.CSV
# ═══════════════════════════════════════════════════════════════════════════

def load_master(path):
    master_primary, master_all = {}, {}
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            parts = [p.strip() for p in line.strip().split(",")]
            if not parts or not parts[0]:
                continue
            mz_key = parts[0]
            names  = [p for p in parts[1:] if p and p.lower() != "n/a"]
            master_primary[mz_key] = names[0] if names else "unknown"
            master_all[mz_key]     = names    if names else ["unknown"]
    return master_primary, master_all

master_primary, master_all = load_master(MASTER_CSV)
print(f"Master loaded: {len(master_primary)} lipid entries.")

# ═══════════════════════════════════════════════════════════════════════════
#  2. AUTO-DISCOVER EXPERIMENTS FROM SUBFOLDER NAMES
# ═══════════════════════════════════════════════════════════════════════════

def discover_experiments(root_dir, proteins_filter=None):
    """
    Scan subfolders matching {protein}_{condition}_{rep} or {protein}_{rep}.
    Apo condition is inferred when there is no condition string (just protein_N).
    Returns dict: {protein: {condition: [rep_folder, ...]}}
    """
    experiments = {}
    pattern = re.compile(r'^(.+?)_(\d+)$')   # greedy protein+condition, then _N

    for entry in sorted(os.scandir(root_dir), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        m = pattern.match(entry.name)
        if not m:
            continue

        prefix = m.group(1)   # everything before the trailing _N
        rep    = int(m.group(2))

        # Split prefix into protein + condition
        # Known proteins — split on first token
        # Strategy: try each protein filter first, else use first underscore token
        protein, condition = None, None
        if proteins_filter:
            for p in proteins_filter:
                if prefix == p:
                    protein, condition = p, "apo"
                    break
                elif prefix.startswith(p + "_"):
                    protein   = p
                    condition = prefix[len(p)+1:]
                    break
        if protein is None:
            # Generic: first token is protein, rest is condition (or apo)
            parts     = prefix.split("_", 1)
            protein   = parts[0]
            condition = parts[1] if len(parts) > 1 else "apo"

        if proteins_filter and protein not in proteins_filter:
            continue

        experiments.setdefault(protein, {}).setdefault(condition, [])
        experiments[protein][condition].append(entry.path)

    return experiments

experiments = discover_experiments(ROOT_DIR, PROTEINS_TO_PROCESS)
print("\nDiscovered experiments:")
for prot, conds in experiments.items():
    for cond, folders in conds.items():
        print(f"  {prot} / {cond}: {len(folders)} replicate(s)")

# ═══════════════════════════════════════════════════════════════════════════
#  3. RESOLVE MS1 FILE FOR EACH PROTEIN+CONDITION
# ═══════════════════════════════════════════════════════════════════════════

def resolve_ms1(root_dir, protein, condition, cond_map=CONDITION_TO_MS1):
    """
    Try to find the MS1 file. Naming convention:
      {protein}_MS1.txt          for apo
      {protein}_{condition}_MS1.txt  for others
    Applies cond_map to expand shortened condition names.
    """
    ms1_cond = cond_map.get(condition, condition)   # expand e.g. Glut→Glutamate

    if condition == "apo":
        candidates = [
            os.path.join(root_dir, f"{protein}_MS1.txt"),
        ]
    else:
        candidates = [
            os.path.join(root_dir, f"{protein}_{ms1_cond}_MS1.txt"),
            os.path.join(root_dir, f"{protein}_{condition}_MS1.txt"),
            os.path.join(root_dir, f"{protein}_MS1.txt"),   # fallback
        ]

    for c in candidates:
        if os.path.isfile(c):
            return c
    return None

# ═══════════════════════════════════════════════════════════════════════════
#  4. LOAD HELPERS
# ═══════════════════════════════════════════════════════════════════════════

def load_lipid_file(path):
    # Read robustly: take only first two columns, skip malformed lines
    rows = []
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    mz  = float(parts[0])
                    val = float(parts[1])
                    rows.append((mz, val))
                except ValueError:
                    continue   # skip header / non-numeric lines
    df = pd.DataFrame(rows, columns=["mz", "intensity"]).sort_values("mz")
    df = df.groupby("mz")["intensity"].sum()   # collapse duplicate m/z rows
    return df

def load_ms1(path):
    rows = []
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
    return pd.DataFrame(rows, columns=["mz", "intensity"]).sort_values("mz")

def bin_ms1_to_grid(ms1_df, grid_index, bin_width=10):
    ms1_binned = pd.Series(0.0, index=grid_index)
    for centre in grid_index:
        mask = ((ms1_df["mz"] >= centre - bin_width / 2) &
                (ms1_df["mz"] <  centre + bin_width / 2))
        ms1_binned[centre] = ms1_df.loc[mask, "intensity"].sum()
    return ms1_binned

def build_profile_matrix(rep_folders):
    """Average lipid intensity profiles across replicate folders."""
    rep_matrices = []
    for folder in rep_folders:
        lipid_files = sorted(glob.glob(os.path.join(folder, "*.txt")))
        if not lipid_files:
            print(f"    WARNING: no .txt files in {folder}")
            continue
        rep_data = {}
        for fp in lipid_files:
            mz_key = os.path.splitext(os.path.basename(fp))[0]
            rep_data[mz_key] = load_lipid_file(fp)
        rep_matrices.append(pd.DataFrame(rep_data).fillna(0).sort_index())

    if not rep_matrices:
        return None

    all_idx = rep_matrices[0].index
    for rm in rep_matrices[1:]:
        all_idx = all_idx.union(rm.index)
    aligned = [rm.reindex(all_idx, fill_value=0) for rm in rep_matrices]
    return pd.concat(aligned).groupby(level=0).mean()

# ═══════════════════════════════════════════════════════════════════════════
#  5. PEARSON r + PERMUTATION NULL
# ═══════════════════════════════════════════════════════════════════════════

def pearson_vs_ms1(profile_matrix, ms1_binned):
    ms1_vec = ms1_binned.reindex(profile_matrix.index, fill_value=0).values
    r_vals, p_vals = {}, {}
    for col in profile_matrix.columns:
        vec = profile_matrix[col].values
        if vec.std() == 0 or ms1_vec.std() == 0:
            r_vals[col], p_vals[col] = 0.0, 1.0
        else:
            r, p = pearsonr(vec, ms1_vec)
            r_vals[col], p_vals[col] = r, p
    return pd.Series(r_vals), pd.Series(p_vals)

def permutation_null_real_data(matrix, ms1_binned, n_perm=N_PERMUTATIONS, seed=42):
    """
    Builds a null distribution of |r| using REAL data only:
    shuffles the real, binned MS1 vector (breaking its true correspondence
    to position), then correlates it against every REAL lipid profile in
    `matrix`. This tests the actual question of interest — "how extreme is
    a lipid's r-value compared to what any real lipid profile would show
    against a randomly-mismatched real MS1 envelope?" — rather than
    comparing against synthetic Gaussian noise.

    Returns a flat list of |r| values (not yet reduced to a threshold),
    so multiple conditions' null distributions can be pooled together
    before taking a percentile.
    """
    rng = np.random.default_rng(seed)
    ms1_vals = ms1_binned.reindex(matrix.index, fill_value=0).values
    null_rs = []
    for _ in range(n_perm):
        shuffled_ms1 = rng.permutation(ms1_vals)
        if shuffled_ms1.std() == 0:
            continue
        for col in matrix.columns:
            vec = matrix[col].values
            if vec.std() == 0:
                continue
            r, _ = pearsonr(vec, shuffled_ms1)
            null_rs.append(abs(r))
    return null_rs

def detect_gap_threshold(r_series, min_r=0.1, gap_quantile=0.95):
    """
    Find the natural break in the r-value distribution.
    Strategy: among lipids with r > min_r, find the largest gap between
    consecutive sorted r values. The threshold is set at the lower edge
    of that gap. Falls back to permutation threshold if no clear gap exists.

    min_r      : ignore near-zero r values (noise floor) when searching
    gap_quantile: only consider gaps in the top gap_quantile of gap sizes
    """
    vals = np.sort(r_series[r_series > min_r].values)
    if len(vals) < 3:
        return None   # not enough points to detect a gap
    gaps = np.diff(vals)
    if gaps.max() < 0.05:
        return None   # no meaningful gap exists
    # Find largest gap
    gap_idx   = np.argmax(gaps)
    threshold = vals[gap_idx]          # lower edge of the gap
    return float(threshold)


def plot_r_distribution(r_series, protein, condition, threshold,
                        threshold_label, out_dir):
    """Histogram of r values with threshold line marked."""
    fig, ax = plt.subplots(figsize=(5, 3.5))
    vals = r_series.dropna().values
    ax.hist(vals, bins=20, color="#886dda", edgecolor="white",
            linewidth=0.5, alpha=0.85)
    ax.axvline(threshold, color="#e05252", linewidth=1.8,
               linestyle="--", label=f"{threshold_label}  r = {threshold:.3f}")
    ax.set_xlabel("Pearson r vs MS1", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title(f"{protein} / {condition}", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, framealpha=0.7)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    safe_cond = condition.replace("/", "_")
    out_path  = os.path.join(out_dir,
                    f"{protein}_{safe_cond}_r_distribution.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
#  6. PLOTTING
# ═══════════════════════════════════════════════════════════════════════════

def plot_heatmap(results_df, protein, threshold, ms1_norm_lookup,
                out_dir, condition_thresholds=None):
    results_sorted = results_df.loc[
        sorted(results_df.index, key=lambda k: float(k))
    ]

    def label(mz_key):
        name   = master_primary.get(mz_key, mz_key)
        ambig  = master_all.get(mz_key, [])
        suffix = " *" if len(ambig) > 1 else ""
        return f"{name}{suffix}\n({mz_key})"

    row_labels = [label(k) for k in results_sorted.index]
    n_rows     = len(results_sorted)
    n_cols     = len(results_sorted.columns)

    fig_h = max(6, n_rows * 0.40)
    fig_w = max(6, n_cols * 1.8 + 3.5)

    fig, (ax_heat, ax_ms1) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [n_cols, 0.55], "wspace": 0.04}
    )

    norm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=1.0)
    im   = ax_heat.imshow(results_sorted.values, cmap=CMAP_HEATMAP,
                          norm=norm, aspect="auto")

    ax_heat.set_xticks(range(n_cols))
    ax_heat.set_xticklabels(
        [c.replace("_", " ") for c in results_sorted.columns],
        fontsize=10, fontweight="bold"
    )
    ax_heat.set_yticks(range(n_rows))
    ax_heat.set_yticklabels(row_labels, fontsize=7.5)
    ax_heat.set_title(
        f"{protein}  —  Pearson r vs MS1 Envelope\n"
        f"● = binder  (r > {threshold:.3f}, permutation threshold, p<0.01)",
        fontsize=10, fontweight="bold", pad=10
    )

    for i, mz_key in enumerate(results_sorted.index):
        for j, cond in enumerate(results_sorted.columns):
            r_val = results_sorted.loc[mz_key, cond]
            if pd.isna(r_val):
                continue
            is_binder = r_val > threshold
            color     = "white" if abs(r_val) > 0.45 else "black"
            txt       = f"{r_val:.2f}" + (" ●" if is_binder else "")
            ax_heat.text(j, i, txt, ha="center", va="center",
                         fontsize=6, color=color,
                         fontweight="bold" if is_binder else "normal")

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.015)
    cbar.set_label("Pearson r", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # MS1 intensity strip
    ms1_vals = np.array([ms1_norm_lookup.get(float(k), 0)
                         for k in results_sorted.index])
    ax_ms1.barh(range(n_rows), ms1_vals, color=COL_BINDER,
                alpha=0.75, edgecolor="none")
    ax_ms1.set_xlim(0, 1.1)
    ax_ms1.set_ylim(-0.5, n_rows - 0.5)
    ax_ms1.invert_yaxis()
    ax_ms1.set_yticks([])
    ax_ms1.set_xlabel("MS1\n(norm.)", fontsize=8)
    ax_ms1.set_xticks([0, 0.5, 1])
    ax_ms1.tick_params(axis="x", labelsize=7)
    ax_ms1.spines[["top", "right", "left"]].set_visible(False)

    fig.tight_layout()
    out_path = os.path.join(out_dir, f"{protein}_pearson_heatmap.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_binder_summary(results_df, protein, threshold, out_dir):
    counts = {c: (results_df[c] > threshold).sum()
              for c in results_df.columns}
    fig, ax = plt.subplots(figsize=(max(4, len(counts) * 1.2), 3.5))
    bars = ax.bar(list(counts.keys()), list(counts.values()),
                  color=COL_BINDER, edgecolor="white", linewidth=0.5)
    ax.bar_label(bars, fontsize=9)
    ax.set_ylabel("Classified binders", fontsize=10)
    ax.set_title(f"{protein}  —  Binders per condition\n"
                 f"(r > {threshold:.3f})", fontsize=10, fontweight="bold")
    ax.set_ylim(0, max(counts.values()) * 1.35 if counts else 5)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    out_path = os.path.join(out_dir, f"{protein}_binder_counts.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")

# ═══════════════════════════════════════════════════════════════════════════
#  7. MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════

os.makedirs(OUT_DIR, exist_ok=True)

for protein, conditions in experiments.items():
    print(f"\n{'═'*60}\n  {protein}\n{'═'*60}")

    all_r                = {}
    condition_thresholds = {}   # {condition: (threshold, label)}
    last_ms1             = None     # keep for MS1 strip normalisation
    condition_data       = {}   # {condition: (matrix, ms1_binned)} — collected in pass 1

    # ── Pass 1: load real data for every condition (no threshold yet) ──
    for condition, rep_folders in conditions.items():
        # Resolve MS1 according to mode
        if MS1_MODE == "apo":
            ms1_path = resolve_ms1(ROOT_DIR, protein, "apo")
            if not ms1_path:
                print(f"  WARNING: no apo MS1 found for {protein} — skipping")
                continue
            print(f"\n  {condition}  →  MS1: {os.path.basename(ms1_path)} (apo, shared)")
        else:
            ms1_path = resolve_ms1(ROOT_DIR, protein, condition)
            if not ms1_path:
                print(f"  WARNING: no MS1 found for {protein}/{condition} — skipping")
                continue
            print(f"\n  {condition}  →  MS1: {os.path.basename(ms1_path)}")

        matrix = build_profile_matrix(rep_folders)
        if matrix is None or matrix.empty:
            print(f"    SKIP — no lipid data.")
            continue

        ms1_raw    = load_ms1(ms1_path)
        ms1_binned = bin_ms1_to_grid(ms1_raw, matrix.index)
        last_ms1   = (ms1_raw, matrix.index)

        condition_data[condition] = (matrix, ms1_binned)

    if not condition_data:
        print(f"  No results for {protein}.")
        continue

    # ── Build one pooled, real-data permutation threshold from ALL conditions ──
    if THRESHOLD_MODE == "hard":
        threshold = HARD_THRESHOLD
        print(f"\n  Hard threshold: r > {threshold}")
    else:
        pooled_null = []
        for condition, (matrix, ms1_binned) in condition_data.items():
            pooled_null.extend(permutation_null_real_data(matrix, ms1_binned))
        threshold = float(np.percentile(pooled_null, PERM_THRESHOLD))
        print(f"\n  Permutation threshold (p<0.01), pooled across "
              f"{len(condition_data)} condition(s), {len(pooled_null)} null values: "
              f"r > {threshold:.4f}")

    # ── Pass 2: compute r-values per condition using the shared threshold ──
    for condition, (matrix, ms1_binned) in condition_data.items():
        r_series, _ = pearson_vs_ms1(matrix, ms1_binned)

        # Use the single, consistent per-protein threshold everywhere (matches heatmap & bar chart)
        cond_threshold = threshold
        thresh_label   = "permutation" if THRESHOLD_MODE == "permutation" else "hard"
        gap_thresh = detect_gap_threshold(r_series)
        if gap_thresh is not None:
            print(f"    (Gap-detected threshold would be r > {gap_thresh:.4f} — not used; "
                  f"using consistent {thresh_label} threshold r > {cond_threshold:.4f} instead)")
        else:
            print(f"    No clear gap found — using {thresh_label} threshold: "
                  f"r > {cond_threshold:.4f}")

        # Store gap threshold per condition for heatmap annotation
        condition_thresholds[condition] = (cond_threshold, thresh_label)

        # Plot r distribution
        plot_r_distribution(r_series, protein, condition,
                            cond_threshold, thresh_label, OUT_DIR)

        all_r[condition] = r_series

        binders = r_series[r_series > cond_threshold].sort_values(ascending=False)
        print(f"    Lipids tested: {len(r_series)}  |  Binders: {len(binders)}")
        for mz_key, r_val in binders.items():
            print(f"      {mz_key:8s}  {master_primary.get(mz_key, mz_key):40s}"
                  f"  r = {r_val:.4f}")

    if not all_r:
        print(f"  No results for {protein}.")
        continue

    results_df = pd.DataFrame(all_r).fillna(np.nan)

    # MS1 normalisation strip — use apo MS1 if available, else last loaded
    ms1_norm_lookup = {}
    if last_ms1:
        ms1_raw_strip, grid = last_ms1
        binned = bin_ms1_to_grid(ms1_raw_strip,
                                  results_df.index.astype(float))
        mx = binned.max()
        if mx > 0:
            ms1_norm_lookup = (binned / mx).to_dict()

    plot_heatmap(results_df, protein, threshold, ms1_norm_lookup,
                 OUT_DIR, condition_thresholds)
    plot_binder_summary(results_df, protein, threshold, OUT_DIR)

    # Export CSV
    out_df = results_df.copy()
    out_df.index.name = "mz_key"
    out_df.insert(0, "lipid_primary",
                  [master_primary.get(k, k) for k in out_df.index])
    out_df.insert(1, "lipid_candidates",
                  [" | ".join(master_all.get(k, [k])) for k in out_df.index])
    out_df["classified_binder_any"] = out_df[list(all_r.keys())].gt(threshold).any(axis=1)
    out_df["threshold_r"]           = threshold
    csv_path = os.path.join(OUT_DIR, f"{protein}_pearson_results.csv")
    out_df.to_csv(csv_path)
    print(f"\n  CSV: {csv_path}")

print(f"\n{'═'*60}\n  Done. Output: {OUT_DIR}\n{'═'*60}")
