#!/bin/bash
#SBATCH --job-name=meta_metasoft
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

################################################################################
# METASOFT Meta‑analysis Template
# Author: Reagan M. Mogire
#
# USAGE (run directly – no wrapper required):
#   bash template_meta_analysis_metasoft.sh ALT
#
#   Where "ALT" (or any phenotype label) substitutes $PHENOTYPE in file names.
#
# Description:
#   1. Decompress and sort per‑cohort summary files by SNP ID.
#   2. Merge them with GNU join, padding missing values with "NA".
#   3. Run Han‑Eskin random‑effects meta‑analysis in METASOFT.
#
# Expectation:
#   Each input file is bgzipped/tab‑delimited and named:
#       path/to/cohort?/sorted_${PHENOTYPE}_data.tbl.gz
#   with columns: <SNP> <A1> <A2> <BETA> <SE> …
#
# Edit the COHORT* variables, METASOFT_JAR, and PVALUE_TABLE to match
# your environment. The analysis logic remains identical to the original.
################################################################################

# ── Input argument ───────────────────────────────────────────────────────────
PHENOTYPE=$1
if [[ -z "$PHENOTYPE" ]]; then
  echo "Usage: bash $(basename $0) <PHENOTYPE>"; exit 1
fi

# ── Generic cohort file paths (edit) ─────────────────────────────────────────
COHORT1="path/to/cohort1/sorted_${PHENOTYPE}_data.tbl.gz"
COHORT2="path/to/cohort2/sorted_${PHENOTYPE}_data.tbl.gz"
COHORT3="path/to/cohort3/sorted_${PHENOTYPE}_data.tbl.gz"

# ── Software locations (edit) ───────────────────────────────────────────────
METASOFT_JAR="/path/to/Metasoft.jar"
PVALUE_TABLE="/path/to/HanEskinPvalueTable.txt"

# ── Output prefix ───────────────────────────────────────────────────────────
OUT_PREFIX="results/${PHENOTYPE}_meta"
mkdir -p results logs

# ── Temporary filenames ─────────────────────────────────────────────────────
TMP1="${PHENOTYPE}_sorted_cohort1.txt"
TMP2="${PHENOTYPE}_sorted_cohort2.txt"
TMP3="${PHENOTYPE}_sorted_cohort3.txt"
MERGED12="${PHENOTYPE}_merged_12.txt"
FINAL_MERGE="${PHENOTYPE}_final_merged.txt"

# ── Step 1: Extract, strip header, sort ──────────────────────────────────────
zcat $COHORT1 | tail -n +2 | sort -k1,1 > $TMP1
zcat $COHORT2 | tail -n +2 | sort -k1,1 > $TMP2
zcat $COHORT3 | tail -n +2 | sort -k1,1 > $TMP3

# ── Step 2: Merge cohort1 + cohort2 ──────────────────────────────────────────
join -a 1 -a 2 -e "NA" -t $'\t' -1 1 -2 1 \
  -o '0,1.6,1.7,2.6,2.7' \
  $TMP1 $TMP2 > $MERGED12

# ── Step 3: Merge with cohort3 ───────────────────────────────────────────────
join -a 1 -a 2 -e "NA" -t $'\t' -1 1 -2 1 \
  -o '0,1.2,1.3,1.4,1.5,2.6,2.7' \
  $MERGED12 $TMP3 > $FINAL_MERGE

# ── Step 4: Run METASOFT ─────────────────────────────────────────────────────
java -Xmx4g -jar $METASOFT_JAR \
  -input $FINAL_MERGE \
  -pvalue_table $PVALUE_TABLE \
  -output ${OUT_PREFIX}_results.txt \
  -log ${OUT_PREFIX}.log \
  -mvalue \
  -mvalue_method sam

# ── Step 5: Clean‑up ─────────────────────────────────────────────────────────
rm -f $TMP1 $TMP2 $TMP3 $MERGED12 $FINAL_MERGE

echo "✓ METASOFT meta‑analysis complete for ${PHENOTYPE}"
