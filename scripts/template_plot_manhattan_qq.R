#!/usr/bin/env Rscript
###############################################################################
# template_plot_manhattan_qq.R
#
# USAGE — run for ONE phenotype file:
#   Rscript template_plot_manhattan_qq.R <PHENOTYPE_LABEL>
#
#   expects a file named:
#       ./gwas/<PHENOTYPE_LABEL>_imputed_for_metasoft.out
#   containing at least the columns:
#       RSID  PVALUE_FE  BETA_FE  STD_FE
#
# OUTPUTS (in working dir):
#   <PHENOTYPE>_manhattan.png
#   <PHENOTYPE>_qq.png
#
# You can loop over many phenotypes in Bash/SLURM, but the script itself
# handles exactly one file per invocation (matching your original workflow).
###############################################################################

## ---- Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(qqman)
})

options(scipen = 999)   # Disable scientific notation in plots

## ---- CLI -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript template_plot_manhattan_qq.R <PHENOTYPE>")
}
pheno <- args[1]

## ---- Load data -------------------------------------------------------------
infile <- file.path("gwas", paste0(pheno, "_imputed_for_metasoft.out"))
if (!file.exists(infile)) stop("File not found: ", infile)

dat <- read.table(infile, header = TRUE, fill = TRUE) |>
  filter(!is.na(PVALUE_FE))

## ---- Tidy ------------------------------------------------------------------
dat <- dat |>
  separate(RSID,
           into = c("CHR", "BP", "Allele1", "Allele2"),
           sep  = "[:_]",
           remove = FALSE,
           convert = TRUE) |>
  mutate(
    CHR     = as.numeric(gsub("^chr", "", CHR, ignore.case = TRUE)),
    POS     = as.numeric(BP),
    P.value = as.numeric(PVALUE_FE),
    SNP     = paste("chr", CHR, POS, Allele1, Allele2, sep = ":")
  )

# Replace zeros with smallest non-zero P
min_p <- min(dat$P.value[dat$P.value > 0], na.rm = TRUE)
dat$P.value[dat$P.value == 0] <- min_p

dat <- dat |>
  filter(is.finite(P.value)) |>
  rename(BETA = BETA_FE, SE = STD_FE) |>
  select(SNP, CHR, POS, Allele1, Allele2, BETA, SE, P.value) |>
  na.omit()

## ---- Manhattan plot --------------------------------------------------------
png(paste0(pheno, "_manhattan.png"), width = 1500, height = 300)
manhattan(
  dat,
  chr = "CHR", bp = "POS", p = "P.value", snp = "SNP",
  col = c("blue2", "steelblue3"),
  ylim = c(0, -log10(1e-20)),
  cex.axis = 1.2, cex.main = 1.5, cex.lab = 1.5
)
dev.off()

## ---- QQ plot ---------------------------------------------------------------
png(paste0(pheno, "_qq.png"), width = 800, height = 800)
qq(dat$P.value, main = paste("QQ Plot –", pheno), ylim = c(0, 50))
dev.off()

cat("✓ Manhattan and QQ plots generated for", pheno, "\n")
