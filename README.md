# Liver-Enzyme GWAS (African-Ancestry) – Template Scripts  
_Companion code for Mogire RM *et al.*_

This repository contains **fully documented templates** that reproduce
every analysis step we ran for the manuscript  
**“Multi-cohort genome-wide association analyses reveal loci underlying circulating liver enzyme levels in African-ancestry populations.”**

> ⚠️ **No data are included.**  
> Each script has placeholder paths (`path/to/...`).  
> Replace those with your own file locations before running.

---

## Directory map

| Path                     | Purpose                                   |
|--------------------------|-------------------------------------------|
| `scripts/`               | Bash / Python / R templates for each step |
| `configs/`               | Optional config files (e.g. FUMA)         |
| `results/ logs/`         | Created at run-time; ignored by Git       |

---

## Template overview

| Stage                              | Script / file |
|------------------------------------|---------------|
| Single-cohort GWAS (SAIGE)         | `scripts/template_run_gwas_saige.sh` |
| Meta-analysis (METASOFT)           | `scripts/template_meta_analysis_metasoft.sh` |
| Fine-mapping (SuSiE + PLINK)       | `scripts/template_finemap_susie.py` |
| eQTL colocalisation (FastENLOC)    | `scripts/template_colocalization_fastenloc.sh` |
| SNP-heritability & rg (LDSC)       | `scripts/template_ldsc_rg.sh` or `template_ldsc_workflow.sh` |
| Manhattan / QQ plots (R)           | `scripts/template_plot_manhattan_qq.R` |

Each file starts with a **usage block** showing the exact command to run after you have edited the paths.

---

## Software / environment

A minimal conda specification is provided in [`environment.yml`](./environment.yml) with:

* Python ≥ 3.10 (pandas, numpy, seaborn, matplotlib, gwaslab, rpy2)  
* R ≥ 4.2 (susieR, data.table, qqman)  
* SAIGE, METASOFT, PLINK 2, LDSC  
* Java 8 (for METASOFT)

> **FastENLOC** must be compiled separately; see its GitHub page.

Create and activate the environment:

```bash
mamba env create -f environment.yml
conda activate liver-gwas
