#!/usr/bin/env python

################################################################################
# SuSiE Fine‑mapping Template (Python + R via rpy2)
# Author: Reagan M. Mogire
#
# USAGE (run directly – no wrapper required):
#   python sussie.py <lead_snp_position> <chromosome> <sample_size>
#
#   Example:
#     python sussie.py 87383538 7 30392
#
#   The script expects:
#     • A bgzipped/tab‑delimited GWAS summary file at
#         "path/to/sorted_<PHENOTYPE>_data.tbl.gz"
#     • A reference FASTA for allele harmonisation (edit REF_FASTA variable)
#     • A PLINK reference panel for LD (edit PLINK_PANEL variable)
#     • The PLINK executable in $PATH or specify PLINK_BIN.
#
# All hard‑coded paths have been replaced with generic placeholders so you can
# adapt them to your own file system. The analysis logic and workflow remain
# unchanged from the original script.
################################################################################

import gwaslab as gl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import subprocess
import sys

# RPy2 interface
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri
numpy2ri.activate()

# ── Parse CLI arguments ──────────────────────────────────────────────────────
if len(sys.argv) != 4:
    print("Usage: python sussie.py <position> <chromosome> <sample_size>")
    sys.exit(1)

lead_pos   = int(sys.argv[1])      # Lead SNP position (base‑1)
chromosome = sys.argv[2]           # Chromosome as string (e.g. "7")
N_sample   = int(sys.argv[3])      # Total sample size for the GWAS

# ── User‑editable file paths ─────────────────────────────────────────────────
PHENOTYPE      = "AST"                                    # edit if needed
SUMSTATS_FILE  = f"path/to/sorted_{PHENOTYPE}_data.tbl.gz"  # bgzipped summary file
REF_FASTA      = "path/to/reference/genome.fasta"          # reference FASTA for harmonise()
PLINK_PANEL    = f"path/to/plink/reference_panel/chr{chromosome}"  # .bed/.bim/.fam prefix
PLINK_BIN      = "plink"                                   # path to PLINK executable

# ── Output settings ──────────────────────────────────────────────────────────
OUT_PREFIX = f"sig_{chromosome}_{lead_pos}"
OUT_DIR    = f"susie_output_chr{chromosome}_{lead_pos}"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Step 1: Load GWAS summary statistics ─────────────────────────────────────
df = pd.read_csv(SUMSTATS_FILE, sep="\t", dtype={"P.value": str})
df["P.value"] = df["P.value"].astype(float)

sumstats = gl.Sumstats(
    df,
    snpid="SNP", chrom="CHR", pos="POS",
    ea="Allele1", nea="Allele2",
    beta="BETA", se="SE", p="P.value",
    build="38")

sumstats.data["N"] = N_sample
sumstats.basic_check()

# ── Step 2: Define ±250 kb region around lead SNP ────────────────────────────
WINDOW = 250_000
start_pos = lead_pos - WINDOW
end_pos   = lead_pos + WINDOW
locus = sumstats.filter_value(f'CHR=={chromosome} & POS>{start_pos} & POS<{end_pos}')
if locus.data.empty:
    print("No SNPs found in the specified region.")
    sys.exit(1)

locus.harmonize(basic_check=False, ref_seq=REF_FASTA)

# Export locus data and SNP list
locus_tsv   = os.path.join(OUT_DIR, f"{OUT_PREFIX}.tsv")
locus_snps  = os.path.join(OUT_DIR, f"{OUT_PREFIX}.snplist.all")
locus.data.to_csv(locus_tsv, sep="\t", index=False)
locus.data["SNPID"].to_csv(locus_snps, sep="\t", index=False, header=False)

# ── Step 3: Compute LD matrices with PLINK ───────────────────────────────────
for rtype, suffix in [("r", ""), ("r2", "_r2")]:
    cmd = [
        PLINK_BIN,
        "--bfile", PLINK_PANEL,
        "--keep-allele-order",
        f"--{rtype}", "square",
        "--write-snplist",
        "--extract", locus_snps,
        "--out", os.path.join(OUT_DIR, f"{OUT_PREFIX}{suffix}")
    ]
    subprocess.run(cmd, check=True)

# ── Step 4: Filter summary stats to LD SNP list ──────────────────────────────
keep_snps = pd.read_csv(os.path.join(OUT_DIR, f"{OUT_PREFIX}.snplist"), header=None, names=["SNPID"])
locus_df  = pd.read_csv(locus_tsv, sep="\t")
filter_df = locus_df[locus_df["SNPID"].isin(keep_snps["SNPID"])]
filter_df = filter_df.set_index("SNPID").reindex(keep_snps["SNPID"]).reset_index()
if filter_df.empty:
    print("Filtered dataset is empty. Aborting.")
    sys.exit(1)

# ── Step 5: Read LD matrices ────────────────────────────────────────────────
LD   = pd.read_csv(os.path.join(OUT_DIR, f"{OUT_PREFIX}.ld"),   sep="\t", header=None).values
LD2  = pd.read_csv(os.path.join(OUT_DIR, f"{OUT_PREFIX}_r2.ld"), sep="\t", header=None).values

# ── Optional: visualise LD
import seaborn as sns; import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(18, 7))
sns.heatmap(LD,  cmap="Spectral", ax=ax[0]); ax[0].set_title("LD r matrix")
sns.heatmap(LD2, cmap="Spectral", ax=ax[1]); ax[1].set_title("LD r² matrix")
fig.savefig(os.path.join(OUT_DIR, f"ld_heatmaps_{chromosome}_{lead_pos}.png"), dpi=300)

# ── Step 6: SuSiE fine‑mapping via rpy2 ──────────────────────────────────────
susieR = importr("susieR")
fit = susieR.susie_rss(
    bhat = filter_df["BETA"].values.reshape((-1, 1)),
    shat = filter_df["SE"].values.reshape((-1, 1)),
    R    = LD,
    L    = 10,
    n    = N_sample
)

# ── Step 7: Annotate credible sets & plot ───────────────────────────────────
cs    = susieR.susie_get_cs(fit, coverage=0.95, min_abs_corr=0.5, Xcorr=LD)[0]
pip   = susieR.susie_get_pip(fit)
filter_df["cs"]  = 0
for i, idx in enumerate(cs):
    filter_df.loc[np.array(idx) - 1, "cs"] = i + 1
filter_df["pip"] = np.array(pip)
filter_df["MLOG10P"] = -np.log10(filter_df["P"])

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(15, 7), gridspec_kw={"height_ratios": [4, 1]})
axes[0].scatter(filter_df["POS"], filter_df["MLOG10P"], c=LD[filter_df["P"].idxmin()]**2, cmap="viridis")
axes[0].set_ylabel("-log10(P)")
axes[1].scatter(filter_df["POS"], filter_df["pip"], c=LD[filter_df["P"].idxmin()]**2, cmap="viridis")
axes[1].set_xlabel("Position (bp)"); axes[1].set_ylabel("PIP")
axes[0].set_xlim(start_pos, end_pos)
plt.tight_layout()
plot_file = os.path.join(OUT_DIR, f"credible_set_plot_{chromosome}_{lead_pos}.pdf")
plt.savefig(plot_file, dpi=300)
print(f"✓ Fine‑mapping plot saved to {plot_file}")


