#!/usr/bin/env Rscript
# Helper: re-run single-feature AUC section of Gene_group_LR_models.Rmd
# Run this AFTER Gene_group_LR_models.Rmd has finished its LR model training,
# to overwrite the buggy all-NA AUC_table_gene_groups_single_features CSV.

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
})

date <- "240426"
sets <- c("set1")
input_folder <- "../../../Data/ML_model_input/Gene_group_models/"

bases <- c("Clinvar_vs_Biobank_long", "Clinvar_vs_Biobank_unique",
           "DNM_vs_Biobank_long",     "DNM_vs_Biobank_unique")

for (base in bases) {
  for (set in sets) {
    obj_name  <- paste0(base, "_", set, "_", date)
    file_path <- file.path(input_folder, paste0(base, "_", set, "_", date, ".csv.gz"))
    assign(obj_name, read_csv(file_path, show_col_types = FALSE))
  }
}

dataframes <- paste0(as.vector(outer(bases, sets, paste, sep = "_")), "_", date)

gene_groups <- c(
  "altfunc_DD_cancer_gene", "altfunc_minus_RASopathy_gene",
  "lof_DD_cancer_gene", "other_DD_cancer_gene", "RASopathy_gene", "all"
)

for (df in dataframes) {
  tmp     <- get(df)
  tmp$all <- TRUE
  assign(df, tmp)
}

AUC_cols_all <- c(
  "n_COSMIC", "n_COSMIC_genome", "n_COSMIC_targeted",
  "n_COSMIC_codon_count", "n_COSMIC_genome_codon_count", "n_COSMIC_targeted_codon_count",
  "hotcodon_cancer",
  "n_sperm", "n_sperm_exomes", "n_sperm_codon_count", "n_sperm_exomes_codon_count",
  "n_buccal", "n_buccal_codon_count",
  "AlphaMissense", "REVEL", "ESM1b",
  "SIFT", "Polyphen2_HVAR", "MutationAssessor", "MetaSVM", "M-CAP",
  "PrimateAI", "VARITY_R", "CADD_raw", "popEVE", "ESM1v"
)

add_AUC_single_features <- function(dat, predictor_col_name) {
  dat <- dat %>%
    mutate(Pathogenicity_category = case_when(
      Pathogenicity_category == "Pathogenic" ~ 1,
      Pathogenicity_category == "Benign"     ~ 0
    ))
  if (length(unique(dat$Pathogenicity_category)) < 2)
    return(list(predictor_col_name, NA, NA, NA, NA, NA))
  r <- roc(response  = dat$Pathogenicity_category,
           predictor = as.numeric(dat[[predictor_col_name]]),
           levels    = c(0, 1), direction = "auto", quiet = TRUE, na.rm = TRUE)
  auc_val  <- round(auc(r)[1], 4)
  ci       <- ci.auc(r)
  dir_sym  <- r$direction
  dir_comp <- if (identical(dir_sym, "<")) "case>control" else "control>case"
  dir_read <- if (identical(dir_sym, "<")) "controls < cases" else "controls > cases"
  list(predictor_col_name, auc_val, round(ci[1], 4), round(ci[3], 4), dir_comp, dir_read)
}

make_AUC_df <- function(my_data, cols) {
  bind_rows(lapply(cols, \(col) {
    r <- add_AUC_single_features(my_data, col)
    data.frame(
      Predictor         = as.character(r[[1]]),
      AUC_full          = as.numeric(r[[2]]),
      AUC_full_lower    = as.numeric(r[[3]]),
      AUC_full_upper    = as.numeric(r[[4]]),
      Case_vs_Control   = as.character(r[[5]]),
      Direction_verbose = as.character(r[[6]]),
      stringsAsFactors  = FALSE
    )
  }))
}

AUC_results_list <- list()

for (df_name in dataframes) {
  df <- get(df_name)
  for (gene_group in gene_groups) {
    df_filtered <- if (gene_group == "all") df else df %>% filter(.data[[gene_group]] == TRUE)
    AUC_df <- make_AUC_df(df_filtered, AUC_cols_all) %>%
      mutate(data_name = df_name, gene_group = gene_group) %>%
      mutate(
        data_name = str_replace(data_name, "^Clinvar_vs_Biobank_long",  "clinvar_long_uncontaminated"),
        data_name = str_replace(data_name, "^Clinvar_vs_Biobank_unique", "clinvar_unique"),
        data_name = str_replace(data_name, "^DNM_vs_Biobank_long",      "DNM_long_uncontaminated"),
        data_name = str_replace(data_name, "^DNM_vs_Biobank_unique",    "DNM_unique"),
        data_name = str_remove(data_name, "_\\d{6}$")
      ) %>%
      mutate(gene_group = case_when(
        gene_group == "altfunc_DD_cancer_gene"           ~ "altfunc",
        gene_group == "altfunc_minus_RASopathy_gene"     ~ "altfunc_minus_RAS",
        gene_group == "lof_DD_cancer_gene"               ~ "lof",
        gene_group == "other_DD_cancer_gene"             ~ "other",
        gene_group == "RASopathy_gene"                   ~ "RASopathy",
        gene_group == "all"                              ~ "all"
      ))
    AUC_results_list <- append(AUC_results_list, list(AUC_df))
  }
}

AUC_results <- bind_rows(AUC_results_list) %>%
  unique() %>%
  extract(data_name,
          into  = c("data_name", "gene_set"),
          regex = "^(.*)_(set\\d+)(?:_.*)?$",
          remove = FALSE)

write_csv(AUC_results,
          paste0("../../../Data/AUC_tables/Gene_group_models/AUC_table_gene_groups_single_features_", date, ".csv"))

cat("Done. Wrote", nrow(AUC_results), "rows to AUC_table_gene_groups_single_features_", date, ".csv\n")
