#!/bin/bash
#SBATCH --job-name=run_gwas_saige
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=24:00:00

################################################################################
# Template: SAIGE single-variant association test
# Author: Reagan M. Mogire
# Description: This template demonstrates how to run a GWAS using SAIGE
#              for a given phenotype, ancestry group, and cohort.
#
# Requirements:
#   - SAIGE (v0.44+)
#   - Genotype file in PLINK or VCF+GRM format
#   - Phenotype file with columns: FID, IID, phenotype
#   - Covariate file with appropriate fixed-effect variables
#
# This script performs:
#   1. Null model fitting
#   2. Single-variant association testing
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# Configurable user parameters (EDIT AS NEEDED)
# ─────────────────────────────────────────────────────────────────────────────

PHENO="ALT"                          # Phenotype name
COHORT="AADM"                        # Cohort name or ancestry label
PHENO_FILE="data/${COHORT}_pheno.txt"  # Tab-delimited file with FID, IID, PHENO
COVAR_FILE="data/${COHORT}_covar.txt"  # Tab-delimited file with FID, IID, covariates
GRM_DIR="genotypes/${COHORT}/grm"      # SAIGE-formatted GRM files
VCF_FILE="genotypes/${COHORT}/chr22.dose.vcf.gz"  # Imputed genotype file
SAMPLE_FILE="${COHORT}.sample"       # Sample list used in VCF/GRM
OUTDIR="results/${PHENO}_${COHORT}"
THREADS=6

# Covariates for model (comma-separated list)
COVARS="age,sex,PC1,PC2,PC3,PC4"

# Create output directory
mkdir -p ${OUTDIR}

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Fit null GLMM model
# ─────────────────────────────────────────────────────────────────────────────

Rscript step1_fitNULLGLMM.R \
  --plinkFile=${GRM_DIR}/${COHORT} \
  --phenoFile=${PHENO_FILE} \
  --phenoCol=${PHENO} \
  --covarColList=${COVARS} \
  --sampleIDColinphenoFile=IID \
  --traitType=quantitative \
  --invNormalize=TRUE \
  --outputPrefix=${OUTDIR}/null_${PHENO}_${COHORT} \
  --nThreads=${THREADS} \
  --LOCO=FALSE

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Single-variant association testing
# ─────────────────────────────────────────────────────────────────────────────

Rscript step2_SPAtests.R \
  --vcfFile=${VCF_FILE} \
  --vcfFileIndex=${VCF_FILE}.tbi \
  --vcfField=DS \
  --chrom=22 \
  --minMAC=5 \
  --sampleFile=${SAMPLE_FILE} \
  --GMMATmodelFile=${OUTDIR}/null_${PHENO}_${COHORT}.rda \
  --varianceRatioFile=${OUTDIR}/null_${PHENO}_${COHORT}.varianceRatio.txt \
  --SAIGEOutputFile=${OUTDIR}/saige_assoc_${PHENO}_${COHORT}.txt \
  --numLinesOutput=2 \
  --IsOutputAFinCaseCtrl=TRUE \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth=0.05 \
  --is_fastTest=TRUE \
  --nThreads=${THREADS}

echo "SAIGE GWAS completed for ${PHENO} in ${COHORT}"
