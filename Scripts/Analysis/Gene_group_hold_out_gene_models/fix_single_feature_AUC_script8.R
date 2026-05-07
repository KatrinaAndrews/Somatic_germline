#!/usr/bin/env Rscript
# Helper: re-run single-feature AUC section of Gene_group_hold_out_LR_models.R
# Run this AFTER Gene_group_hold_out_LR_models.R has finished its LR model training,
# to overwrite the buggy all-NA Single_feature_AUCs_individual_genes CSV.

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
})

date <- "240426"
input_folder <- "../../../Data/ML_model_input/Gene_group_models/"

clinvar_mutation_data_unique <- read_csv(file.path(input_folder, paste0("Clinvar_vs_Biobank_unique_set1_", date, ".csv.gz")), show_col_types = FALSE)
clinvar_mutation_data_long   <- read_csv(file.path(input_folder, paste0("Clinvar_vs_Biobank_long_set1_",   date, ".csv.gz")), show_col_types = FALSE)
DNM_mutation_data_unique     <- read_csv(file.path(input_folder, paste0("DNM_vs_Biobank_unique_set1_",     date, ".csv.gz")), show_col_types = FALSE)
DNM_mutation_data_long       <- read_csv(file.path(input_folder, paste0("DNM_vs_Biobank_long_set1_",       date, ".csv.gz")), show_col_types = FALSE)

data_sets <- c("clinvar_mutation_data_unique", "clinvar_mutation_data_long",
               "DNM_mutation_data_unique",     "DNM_mutation_data_long")

for (dataset in data_sets) {
  dat <- get(dataset) %>%
    mutate(
      n_COSMIC_codon_count_log = log1p(n_COSMIC_codon_count),
      n_sperm_codon_count_log  = log1p(n_sperm_codon_count),
      n_buccal_codon_count_log = log1p(n_buccal_codon_count)
    )
  assign(dataset, dat, envir = .GlobalEnv)
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
  reverse_dir <- c("popEVE", "ESM1b", "ESM1v", "SIFT")
  chosen_dir  <- if (predictor_col_name %in% reverse_dir) ">" else "<"
  dat <- dat %>%
    mutate(Pathogenicity_category = case_when(
      Pathogenicity_category == "Pathogenic" ~ 1,
      Pathogenicity_category == "Benign"     ~ 0
    ))
  if (length(unique(dat$Pathogenicity_category)) < 2)
    return(list(predictor_col_name, NA, NA, NA, NA))
  r <- roc(response  = dat$Pathogenicity_category,
           predictor = as.numeric(dat[[predictor_col_name]]),
           levels    = c(0, 1), direction = chosen_dir, quiet = TRUE, na.rm = TRUE)
  ci      <- ci.auc(r)
  dir_read <- if (identical(r$direction, "<")) "case>control" else "control>case"
  list(predictor_col_name, round(auc(r)[1], 4), round(ci[1], 4), round(ci[3], 4), dir_read)
}

make_AUC_df <- function(my_data, cols) {
  bind_rows(lapply(cols, \(col) {
    r <- add_AUC_single_features(my_data, col)
    data.frame(
      Predictor      = as.character(r[[1]]),
      AUC_full       = as.numeric(r[[2]]),
      AUC_full_lower = as.numeric(r[[3]]),
      AUC_full_upper = as.numeric(r[[4]]),
      Direction      = as.character(r[[5]]),
      stringsAsFactors = FALSE
    )
  }))
}

single_feature_AUC_results <- data.frame()

for (df_name in data_sets) {
  df      <- get(df_name)
  nm_parts <- str_split(df_name, "_")[[1]]
  dataset <- paste0(nm_parts[1], "_", nm_parts[4])

  genes <- df %>%
    group_by(gene) %>%
    count(Pathogenicity_category) %>%
    pivot_wider(values_from = n, names_from = Pathogenicity_category) %>%
    filter(coalesce(Benign, 0) > 9, coalesce(Pathogenic, 0) > 9) %>%
    pull(gene)

  for (g in genes) {
    gene_df <- filter(df, gene == g)
    pat_df  <- filter(gene_df, Pathogenicity_category == "Pathogenic")
    if (nrow(pat_df) == 0) next

    missing_pct <- pat_df %>%
      summarise(across(all_of(AUC_cols_all), ~ mean(is.na(.)) * 100)) %>%
      as.list() %>% unlist()

    AUC_cols_subset <- names(missing_pct)[!is.na(missing_pct) & missing_pct < 20]

    AUC_cols_subset2 <- gene_df %>%
      select(Pathogenicity_category, any_of(AUC_cols_subset)) %>%
      pivot_longer(-Pathogenicity_category, names_to = "predictor",
                   values_to = "score", values_drop_na = TRUE) %>%
      count(predictor, Pathogenicity_category, name = "n") %>%
      pivot_wider(names_from = Pathogenicity_category, values_from = n, values_fill = 0) %>%
      filter(coalesce(Benign, 0) >= 10, coalesce(Pathogenic, 0) >= 10) %>%
      pull(predictor)

    if (length(AUC_cols_subset2) > 0) {
      AUC_df <- make_AUC_df(gene_df, AUC_cols_subset2) %>%
        mutate(dataset = dataset, gene = g)
      single_feature_AUC_results <- bind_rows(single_feature_AUC_results, AUC_df)
    }
  }
}

write_csv(single_feature_AUC_results,
          paste0("../../../Data/AUC_tables/Gene_group_hold_out_models/Single_feature_AUCs_individual_genes_", date, ".csv"))

cat("Done. Wrote", nrow(single_feature_AUC_results), "rows to Single_feature_AUCs_individual_genes_", date, ".csv\n")
