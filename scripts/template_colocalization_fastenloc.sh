#!/bin/bash
#SBATCH --job-name=fastenloc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

################################################################################
# FastENLOC Colocalization Template
# Author: Reagan M. Mogire
#
# USAGE:
#   bash template_colocalization_fastenloc.sh ALT
#
#   ALT (or any phenotype label) replaces $PHENOTYPE in paths below.
#
# Description:
#   1. Converts GWAS summary stats to FastENLOC format (locus_id SNPID Beta SE).
#   2. Extracts total variant count from a user-provided log file (optional).
#   3. Runs FastENLOC colocalization against a single GTEx tissue (default Liver).
#
# Inputs (edit for your environment):
#   - sorted_${PHENOTYPE}_data.tbl.gz : bgzipped GWAS summary table
#   - GTEx eQTL VCF or annotation (.vcf[.gz])
#   - FastENLOC executable path
#
# Output:
#   - Colocalization results with chosen tissue prefix
################################################################################

PHENOTYPE=$1
if [[ -z "$PHENOTYPE" ]]; then
  echo "Usage: bash $(basename $0) <PHENOTYPE>"; exit 1
fi

# ─── User paths (edit) ───────────────────────────────────────────────────────
GWAS_FILE="path/to/sorted_${PHENOTYPE}_data.tbl.gz"    # bgzipped GWAS table
LOG_FILE="path/to/${PHENOTYPE}_meta_analysis.log"      # optional: log with SNP count line
EQTL_VCF="path/to/GTEx_Liver.vcf.gz"                  # liver eQTL VCF
FASTENLOC_BIN="/path/to/fastenloc"                   # FastENLOC executable
TISSUE="Liver"
OUT_PREFIX="fastenloc_${PHENOTYPE}_${TISSUE}"

# ─── Create output dir ───────────────────────────────────────────────────────
mkdir -p results
cd results

# ─── Step 1: Convert GWAS to FastENLOC format ───────────────────────────────
GWAS_OUT="${PHENOTYPE}_gwas_beta_se.tsv"

echo "Converting ${GWAS_FILE} → ${GWAS_OUT} …"
zcat "$GWAS_FILE" | awk 'NR>1 {print $2":"$3, $1, $6, $7}' > "$GWAS_OUT"

# ─── Step 2: Get total variant count ─────────────────────────────────────────
if [[ -f "$LOG_FILE" ]]; then
  TOTAL_VAR=$(awk -F': ' "/Number of SNPs analyzed:/ {print \$2}" "$LOG_FILE")
else
  TOTAL_VAR=$(wc -l < "$GWAS_OUT")
fi

echo "Total variants = $TOTAL_VAR"

# ─── Step 3: Run FastENLOC ───────────────────────────────────────────────────
"$FASTENLOC_BIN" \
  -gs "$GWAS_OUT" \
  -total_variants "$TOTAL_VAR" \
  -eqtl "$EQTL_VCF" \
  -tissue "$TISSUE" \
  -thread 8 \
  -prefix "$OUT_PREFIX"

echo "✓ FastENLOC completed: ${OUT_PREFIX}*"
