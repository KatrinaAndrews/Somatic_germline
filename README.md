# Somatic mutation data for germline variant classification in developmental disorder genes

This repository contains all analysis scripts associated with the manuscript:

> **[Using somatic data to aid germline clinical variant
interpretation in developmental disorders]**
> [Katrina A Andrews, Matthew DC Neville, Iñigo Martincorena, Raheleh Rahbari, Helen
Firth, Sarah Lindsay, Marc Tischkowitz, and Matthew E Hurles]
> Unpublished

---

## Overview

This project investigates whether somatic mutation data from cancer (COSMIC), sperm Nanoseq, and buccal Nanoseq can be used as evidence to support germline variant classification in dominant developmental disorder (DD) genes. We show that genes driving somatic clonal selection disproportionately accumulate germline gain-of-function variants, and that codon-level somatic counts provide likelihood ratio evidence applicable to ACMG/AMP variant classification.

Key analyses:
- Benchmarking of somatic features and computational predictors for discriminating pathogenic germline variants from benign population variants, across gene groups (RASopathy genes, altered-function DD genes, loss-of-function DD genes)
- Logistic regression models combining REVEL with somatic codon-level counts
- Likelihood ratios for COSMIC missense codon count thresholds, mapped to ACMG/AMP evidence levels
- Benchmarking against deep mutational scanning (MAVE) data for PTPN11, PTEN, DDX3X, and HRAS
- Application to automatic ACMG reclassification of RASopathy VUS and conflicting variants

---

## Repository structure

```
Scripts/
├── Data_processing_and_annotation/   # Variant annotation and data preparation
├── Analysis/
│   ├── Gene_group_models/            # Gene-group-level benchmarking models
│   ├── Gene_group_hold_out_gene_models/  # Per-gene hold-out models
│   ├── Other_visualisations/         # Upset plot gene overlap figure
│   ├── Figures_for_somatic_germline_paper.Rmd   # Generates all manuscript figures
│   ├── Tables_for_somatic_germline_paper.Rmd    # Generates all supplementary tables
│   ├── Numbers_for_somatic_germline_paper.Rmd   # Computes all inline manuscript numbers
│   └── gene_panel_functions_for_manuscript_figures.R  # Shared plotting functions
├── MAVE_benchmarking/                # MAVE vs somatic data benchmarking (PTPN11, PTEN, DDX3X, HRAS)
└── ACMG_and_LRs/                     # Likelihood ratio calculations and ACMG classification
```

---

## Data availability

Publicly available data used in this analysis can be accessed as follows:

| Dataset | Source |
|---|---|
| ClinVar pathogenic/benign/VUS missense variants | ClinVar (NCBI) |
| De novo mutations (31k cohort) | Kaplanis et al. 2020, *Nature* (doi:10.1038/s41586-020-2832-5) |
| UK Biobank population variants | UK Biobank |
| COSMIC somatic mutations | COSMIC v98 (cancer.sanger.ac.uk) |
| TwinsUK sperm Nanoseq | Neville et al. 2025, *Nature* (doi:10.1038/s41586-025-09448-3) |
| Buccal Nanoseq | Lawson et al. 2025, *Nature* (doi:10.1038/s41586-025-09584-w) |
| DDG2P gene list | GENE2PHENOTYPE (ebi.ac.uk/gene2phenotype) |
| RASopathy VCEP gene list | ClinGen |
| PTPN11 MAVE scores | Jiang et al. 2025, *Nature Communications* (doi:10.1038/s41467-025-60641-3) |
| PTEN MAVE scores (activity) | Mighell et al. 2018, *Am. J. Hum. Genet.* (doi:10.1016/j.ajhg.2018.03.018) |
| PTEN MAVE scores (abundance) | Matreyek et al. 2018, *Nature Genetics* (doi:10.1038/s41588-018-0122-z), obtained via the integrated dataset of Matreyek et al. 2021, *Genome Medicine* (doi:10.1186/s13073-021-00984-x) |
| DDX3X MAVE scores | Radford et al. 2023, *Nature Communications* (doi:10.1038/s41467-023-43041-4) |
| HRAS MAVE scores | Hidalgo et al. 2022, *eLife* (doi:10.7554/eLife.76595) |

---

## Dependencies

All analyses were performed in **R 4.5**. The following packages are required:

```r
install.packages(c(
  "tidyverse", "caret", "patchwork", "ggrepel", "ggtext",
  "ggnewscale", "cowplot", "janitor", "UpSetR", "DT",
  "knitr", "rmarkdown", "ragg", "parallel"
))
```

---

## Execution order

Scripts must be run in the following order.

**1. Data processing** *(cluster; skip if using pre-processed data)*
- `Scripts/Data_processing_and_annotation/Annotating_with_clean_NumberSubmissions.Rmd`
- `Scripts/Data_processing_and_annotation/Define_cancer_DD_gene_mechanism_sets.Rmd`
- `Scripts/Data_processing_and_annotation/Define_somatic_selection_gene_groups.Rmd`
- `Scripts/Data_processing_and_annotation/Data_prep_gene_groups.Rmd`

**2. Benchmarking models**
- `Scripts/Analysis/Gene_group_models/Gene_group_single_feature_AUC_AUPRC_F1_MCC.Rmd`
- `Scripts/Analysis/Gene_group_models/Gene_group_LR_models.Rmd`
- `Scripts/Analysis/Gene_group_hold_out_gene_models/Gene_group_hold_out_single_feature_AUC_AUPRC_F1_MCC.Rmd`
- `Scripts/Analysis/Gene_group_hold_out_gene_models/Gene_group_hold_out_LR_models_parallel.R`

**3. MAVE benchmarking**
- `Scripts/MAVE_benchmarking/PTPN11_variant_benchmarking.Rmd`
- `Scripts/MAVE_benchmarking/PTEN_variant_benchmarking.Rmd`
- `Scripts/MAVE_benchmarking/DDX3X_variant_benchmarking.Rmd`
- `Scripts/MAVE_benchmarking/HRAS_variant_benchmarking.Rmd`

**4. Likelihood ratios and ACMG classification**
- `Scripts/ACMG_and_LRs/Likelihood_ratio_threshold_calcs.Rmd`
- `Scripts/ACMG_and_LRs/Likelihood_ratio_threshold_plots.Rmd`
- `Scripts/ACMG_and_LRs/Prep_data_for_Auto_ACMG.Rmd`
- `Scripts/ACMG_and_LRs/Auto_ACMG_RASopathy_results_analysis.Rmd`

**5. Manuscript outputs**
- `Scripts/Analysis/Figures_for_somatic_germline_paper.Rmd` — all figures
- `Scripts/Analysis/Tables_for_somatic_germline_paper.Rmd` — all supplementary tables
- `Scripts/Analysis/Numbers_for_somatic_germline_paper.Rmd` — all inline numbers

---

## Contact

[Katrina Andrews] — [ka13@sanger.ac.uk]

Wellcome Sanger Institute
