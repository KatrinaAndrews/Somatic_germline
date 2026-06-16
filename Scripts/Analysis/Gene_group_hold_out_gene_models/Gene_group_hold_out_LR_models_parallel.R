#!/usr/bin/env Rscript
# Gene-group hold-out LR models (parallelised)
#
# Single-feature AUC/AUPRC/F1/MCC is computed separately in:
#   Gene_group_hold_out_single_feature_AUC_AUPRC_F1_MCC.Rmd
# Run that script first to produce Single_feature_AUCs_individual_genes_<date>.csv.
#
# Parallelisation strategy: mclapply (fork-based) over independent units.
#   - Hold-out LR models: parallelised over (dataset × gene_group × gene × CEP) combinations.
#   allowParallel = FALSE in trainControl prevents nested parallelism inside caret.
#
# Outputs:
#   Data/AUC_tables/Gene_group_hold_out_models/AUC_table_gene_holdout_LR_models_<date>.csv
#   Data/AUC_tables/Gene_group_hold_out_models/Variant_counts_per_gene_<date>.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(caret)
  library(parallel)
})

date    <- "240426"
n_cores <- max(1L, detectCores() - 1L)
message("Using ", n_cores, " cores")

dir.create("../../../Data/AUC_tables/Gene_group_hold_out_models",
           recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Model training parameters
# ---------------------------------------------------------------------------
set.seed(42)
seeds <- vector(mode = "list", length = 101)
for (i in 1:100) seeds[[i]] <- sample.int(1000, 1)
seeds[[101]] <- sample.int(1000, 1)

twoClassSummary_pathogenic <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) lev <- levels(data$obs)
  roc_obj <- pROC::roc(response  = data$obs,
                       predictor = data[, "Pathogenic"],
                       levels    = c("Benign", "Pathogenic"), quiet = TRUE)
  c(ROC  = as.numeric(pROC::auc(roc_obj)),
    Sens = caret::sensitivity(data$pred, data$obs, positive = "Pathogenic"),
    Spec = caret::specificity(data$pred, data$obs, negative = "Benign"))
}

# allowParallel = FALSE: we parallelise at the outer loop, not inside caret
trctrl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 10, sampling = "down",
  summaryFunction = twoClassSummary_pathogenic, classProbs = TRUE,
  seeds = seeds, allowParallel = FALSE
)

# ---------------------------------------------------------------------------
# Helper: train one gene hold-out model and return AUC rows
# ---------------------------------------------------------------------------
train_gene_holdout <- function(df_minus_gene, df_gene, cep, dataset, gene_group, g) {
  model_dat <- df_minus_gene %>%
    select(Pathogenicity_category, all_of(cep),
           n_COSMIC_codon_count_log, n_sperm_codon_count_log, n_buccal_codon_count_log) %>%
    mutate(Pathogenicity_category = factor(Pathogenicity_category,
                                           levels = c("Benign", "Pathogenic"))) %>%
    drop_na()

  lr_full <- train(Pathogenicity_category ~ ., data = model_dat,
                   method = "glm", family = "binomial", trControl = trctrl, metric = "ROC")

  lr_cep  <- train(as.formula(paste0("Pathogenicity_category ~ ", cep)), data = model_dat,
                   method = "glm", family = "binomial", trControl = trctrl, metric = "ROC")

  cf <- coef(lr_full$finalModel)
  getc <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else NA_real_
  coef_tbl <- tibble(
    coef_cep    = getc(cep),
    coef_COSMIC = getc("n_COSMIC_codon_count_log"),
    coef_sperm  = getc("n_sperm_codon_count_log"),
    coef_buccal = getc("n_buccal_codon_count_log")
  )

  cf2 <- coef(lr_cep$finalModel)
  coef_cep_only <- if (cep %in% names(cf2)) unname(cf2[[cep]]) else NA_real_

  probs_tr  <- predict(lr_full, model_dat, type = "prob")[, "Pathogenic"]
  roc_tr    <- roc(model_dat$Pathogenicity_category, probs_tr, direction = "<", quiet = TRUE)
  ci_tr     <- ci.auc(roc_tr)

  auc_train_row <- bind_cols(
    tibble(dataset = dataset, gene_group = gene_group, CEP = cep, gene = g,
           test_or_train = "train",
           AUC_full       = round(auc(roc_tr), 3),
           AUC_full_lower = round(ci_tr[1], 3),
           AUC_full_upper = round(ci_tr[3], 3),
           coef_cep_only  = coef_cep_only),
    coef_tbl
  )

  test_dat <- df_gene %>%
    select(Pathogenicity_category, all_of(cep),
           n_COSMIC_codon_count_log, n_sperm_codon_count_log, n_buccal_codon_count_log) %>%
    mutate(Pathogenicity_category = factor(Pathogenicity_category,
                                           levels = c("Benign", "Pathogenic"))) %>%
    drop_na()

  probs_ts_full <- predict(lr_full, test_dat, type = "prob")[, "Pathogenic"]
  probs_ts_cep  <- predict(lr_cep,  test_dat, type = "prob")[, "Pathogenic"]
  roc_ts_full   <- roc(test_dat$Pathogenicity_category, probs_ts_full, direction = "<", quiet = TRUE)
  roc_ts_cep    <- roc(test_dat$Pathogenicity_category, probs_ts_cep,  direction = "<", quiet = TRUE)
  ci_full       <- ci.auc(roc_ts_full)
  ci_cep        <- ci.auc(roc_ts_cep)
  dl            <- roc.test(roc_ts_full, roc_ts_cep, method = "delong")

  auc_test_row <- bind_cols(
    tibble(dataset = dataset, gene_group = gene_group, CEP = cep, gene = g,
           test_or_train  = "test",
           AUC_full       = round(auc(roc_ts_full), 3),
           AUC_full_lower = round(ci_full[1], 3), AUC_full_upper = round(ci_full[3], 3),
           AUC_cep_only   = round(auc(roc_ts_cep), 3),
           AUC_cep_lower  = round(ci_cep[1], 3),  AUC_cep_upper  = round(ci_cep[3], 3),
           delong_p_value = signif(dl$p.value, 3),
           AUC_diff       = round(dl$estimate[[1]] - dl$estimate[[2]], 3),
           AUC_diff_lower = round(dl$conf.int[[1]], 3),
           AUC_diff_upper = round(dl$conf.int[[2]], 3)),
    coef_tbl
  )

  bind_rows(auc_train_row, auc_test_row)
}

# ---------------------------------------------------------------------------
# Build all valid (data_set, gene_group, gene, cep) combinations upfront,
# filtering on missingness and n criteria, then parallelise over them.
# ---------------------------------------------------------------------------
gene_groups <- c("RASopathy_gene", "altfunc_DD_cancer_gene",
                 "lof_DD_cancer_gene", "other_DD_cancer_gene")
ceps        <- c("AlphaMissense", "REVEL", "ESM1b")

lr_combos <- do.call(rbind, lapply(data_sets, function(data_set) {
  df       <- get(data_set, envir = .GlobalEnv)
  nm_parts <- str_split(data_set, "_")[[1]]
  dataset  <- paste0(nm_parts[1], "_", nm_parts[4])

  do.call(rbind, lapply(gene_groups, function(gene_group) {
    df_gg <- filter(df, !!sym(gene_group) == TRUE)

    genes <- df_gg %>%
      group_by(gene) %>%
      count(Pathogenicity_category) %>%
      pivot_wider(values_from = n, names_from = Pathogenicity_category) %>%
      filter(coalesce(Benign, 0) > 9, coalesce(Pathogenic, 0) > 9) %>%
      pull(gene)

    if (length(genes) == 0) return(NULL)

    do.call(rbind, lapply(genes, function(g) {
      do.call(rbind, lapply(ceps, function(cep) {
        gene_df <- filter(df_gg, gene == g)

        missingness <- gene_df %>%
          filter(Pathogenicity_category == "Pathogenic") %>%
          summarize(na_pct = mean(is.na(.data[[cep]])) * 100) %>%
          pull(na_pct)

        n_counts <- gene_df %>%
          filter(!is.na(.data[[cep]])) %>%
          count(Pathogenicity_category) %>%
          pivot_wider(values_from = n, names_from = Pathogenicity_category)

        n_benign <- coalesce(n_counts$Benign[1],     0L)
        n_path   <- coalesce(n_counts$Pathogenic[1], 0L)

        if (!is.na(missingness) && missingness < 20 && n_benign > 9 && n_path > 9) {
          data.frame(data_set = data_set, dataset = dataset, gene_group = gene_group,
                     gene = g, cep = cep, stringsAsFactors = FALSE)
        } else {
          NULL
        }
      }))
    }))
  }))
}))

run_one_holdout <- function(idx) {
  combo    <- lr_combos[idx, ]
  df_gg    <- filter(get(combo$data_set, envir = .GlobalEnv), !!sym(combo$gene_group) == TRUE)
  df_gene  <- filter(df_gg, gene == combo$gene)
  df_minus <- filter(df_gg, gene != combo$gene)
  tryCatch(
    train_gene_holdout(df_minus, df_gene, combo$cep, combo$dataset, combo$gene_group, combo$gene),
    error = function(e) {
      message("Error for ", combo$gene, "/", combo$cep, ": ", conditionMessage(e))
      NULL
    }
  )
}

message("Training ", nrow(lr_combos), " hold-out model combinations across ", n_cores, " cores...")
lr_results_list <- mclapply(seq_len(nrow(lr_combos)), run_one_holdout, mc.cores = n_cores)
LR_auc_results  <- bind_rows(Filter(Negate(is.null), lr_results_list))

write_csv(LR_auc_results,
          paste0("../../../Data/AUC_tables/Gene_group_hold_out_models/AUC_table_gene_holdout_LR_models_", date, ".csv"))
message("Wrote ", nrow(LR_auc_results), " LR AUC rows")

# ---------------------------------------------------------------------------
# Variant counts per gene
# ---------------------------------------------------------------------------
variant_count_DF <- do.call(rbind, lapply(data_sets, function(data_set) {
  df    <- get(data_set, envir = .GlobalEnv)
  genes <- unique(df$gene)
  do.call(rbind, lapply(genes, function(g) {
    df_gene <- filter(df, gene == g)
    data.frame(
      data_set          = data_set,
      gene              = g,
      n_path_variants   = sum(df_gene$Pathogenicity_category == "Pathogenic"),
      n_benign_variants = sum(df_gene$Pathogenicity_category == "Benign")
    )
  }))
}))

write_csv(variant_count_DF,
          paste0("../../../Data/AUC_tables/Gene_group_hold_out_models/Variant_counts_per_gene_", date, ".csv"))
