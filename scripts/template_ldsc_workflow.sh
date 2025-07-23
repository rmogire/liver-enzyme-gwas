#!/bin/bash
#SBATCH --job-name=ldsc_workflow
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

################################################################################
# LDSC Workflow Template
# Author:  Reagan M. Mogire
#
# USAGE (one trait-pair per run):
#   bash template_ldsc_workflow.sh  ALP  ALT  32840  28989
#
#   where
#     TRAIT1 = ALP         (label used in filenames)
#     TRAIT2 = ALT         (second trait label)
#     N1     = 32840       (total sample size for trait1)
#     N2     = 28989       (total sample size for trait2)
#
# INPUT EXPECTATION
# ─────────────────
#   • Raw GWAS summary files (bgzipped)             :  raw/<TRAIT>.tsv.gz
#       must contain at least  (CHR  BP  A1  A2  BETA  SE  P)
#   • 1000G/UKBB LD-score reference (AFR/EUR/…)     :  ref/UKBB.AFR.*
#   • ldsc.py and munge_sumstats.py in $PATH
#
# OUTPUT
# ──────
#   munged/ <TRAIT>.sumstats.gz        (formatted for LDSC)
#   results/<TRAIT>_h2.log             (heritability)
#   results/<T1>_<T2>_rg.log           (genetic correlation)
#
################################################################################

# ── Command-line arguments ────────────────────────────────────────────────────
TRAIT1=$1   # e.g. ALP
TRAIT2=$2   # e.g. ALT
N1=$3       # sample size for trait1
N2=$4       # sample size for trait2
[[ -z "$TRAIT1" || -z "$TRAIT2" || -z "$N1" || -z "$N2" ]] && {
  echo "Usage: bash $(basename $0) TRAIT1 TRAIT2 N1 N2" ; exit 1 ; }

# ── Directories (edit these) ─────────────────────────────────────────────────
RAW_DIR="raw"                # bgzipped raw GWAS tables
MUNGED_DIR="munged"          # munged .sumstats.gz
REF_DIR="ref/UKBB.AFR"       # basename for *.l2.ldscore.gz and *.weights.gz
OUT_DIR="results"
mkdir -p $MUNGED_DIR $OUT_DIR  logs

# ── Helper: Reformat function ────────────────────────────────────────────────
reformat() {
  local TRAIT=$1
  echo "→ Re-formatting $TRAIT"
  zcat ${RAW_DIR}/${TRAIT}.tsv.gz | \
  awk 'BEGIN{OFS="\t"}
       NR==1 {print "SNP","A1","A2","BETA","SE","P"; next}
       {print $1":"$2":"$3":"$4, $3, $4, $5, $6, $7}' \
  > ${MUNGED_DIR}/${TRAIT}_munge_input.txt
}

# ── Step 0: Reformat raw files for both traits ───────────────────────────────
for T in $TRAIT1 $TRAIT2; do reformat $T; done

# ── Step 1: Munge (format for LDSC) ──────────────────────────────────────────
module load ldsc   # or ensure ldsc env is activated

munge_sumstats.py --sumstats ${MUNGED_DIR}/${TRAIT1}_munge_input.txt \
                  --out      ${MUNGED_DIR}/${TRAIT1} \
                  --N $N1

munge_sumstats.py --sumstats ${MUNGED_DIR}/${TRAIT2}_munge_input.txt \
                  --out      ${MUNGED_DIR}/${TRAIT2} \
                  --N $N2

# ── Step 2: SNP-based heritability for each trait ────────────────────────────
ldsc.py --h2 ${MUNGED_DIR}/${TRAIT1}.sumstats.gz \
        --ref-ld-chr ${REF_DIR} \
        --w-ld-chr   ${REF_DIR} \
        --out ${OUT_DIR}/${TRAIT1}_h2

ldsc.py --h2 ${MUNGED_DIR}/${TRAIT2}.sumstats.gz \
        --ref-ld-chr ${REF_DIR} \
        --w-ld-chr   ${REF_DIR} \
        --out ${OUT_DIR}/${TRAIT2}_h2

# ── Step 3: Pair-wise genetic correlation ────────────────────────────────────
ldsc.py --rg  ${MUNGED_DIR}/${TRAIT1}.sumstats.gz,\
${MUNGED_DIR}/${TRAIT2}.sumstats.gz \
        --ref-ld-chr ${REF_DIR} \
        --w-ld-chr   ${REF_DIR} \
        --out ${OUT_DIR}/${TRAIT1}_vs_${TRAIT2}_rg

echo "✓ LDSC workflow complete for $TRAIT1 and $TRAIT2"
