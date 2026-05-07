# ACMG_auto.R — missense-only, dominant gene ACMG classification (no somatic weighting)
#
# Applies Bayesian ACMG framework using PS1, PS2, PM1, PM2, PP2, PP3, BA1, BS1, BP4.
# VCEP rules: RASopathy VCEP allele frequency cutoffs; PP2 applied to all missense.
#
# Usage:
#   Rscript ACMG_auto.R <variants.tsv> <hotspots.tsv> <output.tsv>
#
# Arguments:
#   variants.tsv   — annotated variant TSV (see Prep_data_for_Auto_ACMG.Rmd for format)
#   hotspots.tsv   — two-column TSV with gene and codon columns
#   output.tsv     — path to write classified variants

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript ACMG_auto.R <variants.tsv> <hotspots.tsv> <output.tsv>")
variants_path <- args[1]
hotspots_path <- args[2]
output        <- args[3]

# RASopathy genes where PP2_supporting is generally applicable (all RASopathy VCEP genes)
pp2_genes <- c("PTPN11","SOS1","SOS2","RAF1","BRAF","KRAS","HRAS","NRAS",
               "MAP2K1","MAP2K2","RIT1","SHOC2")

# --- Parameters: prior & LRs (Tavtigian odds-path) ---
prior_p    <- 0.10
prior_odds <- prior_p / (1 - prior_p)
LR <- list(PS = 18.7, PM = 4.33, PP = 2.08,
           VSb = 1/350, BS = 1/18.7, BP = 1/2.08)

# Thresholds
thr <- list(
  pp3_revel = 0.70,
  bp4_revel = 0.30,
  ba1 = 0.0005, bs1_lo = 0.00025, bs1_hi = 0.0005  # allele frequency cutoffs from RASopathy VCEP
)

# --- Inputs ---
hotspots <- read_tsv(hotspots_path, show_col_types = FALSE) %>%
  mutate(codon = as.integer(codon)) %>%
  distinct(gene, codon) %>% mutate(hotspot = TRUE)

v <- read_tsv(variants_path, show_col_types = FALSE) %>%
  mutate(
    codon                        = as.integer(codon),
    gnomad_popmax_af             = as.numeric(gnomad_popmax_af),
    revel                        = as.numeric(revel),
    cadd_phred                   = as.numeric(cadd_phred),
    de_novo_confirmed            = as.logical(de_novo_confirmed),
    clinvar_exact_path_star2plus = as.logical(clinvar_exact_path_star2plus),
    gnomad_min_cov               = suppressWarnings(as.numeric(gnomad_min_cov))
  ) %>%
  left_join(hotspots, by = c("gene", "codon")) %>%
  mutate(hotspot = coalesce(hotspot, FALSE),
         cov_ok  = is.na(gnomad_min_cov) | gnomad_min_cov >= 10) %>%
  mutate(PP2_supporting = grepl("missense", consequence))

# --- Evidence flags (dominant, missense) ---
v2 <- v %>%
  mutate(
    PS1 = coalesce(clinvar_exact_path_star2plus, FALSE),
    PS2 = coalesce(de_novo_confirmed, FALSE),
    PM1 = hotspot,
    PP2 = coalesce(PP2_supporting, FALSE),
    PM2 = (is.na(gnomad_popmax_af) | gnomad_popmax_af == 0),
    PP3 = revel >= thr$pp3_revel,
    BA1 = cov_ok & (gnomad_popmax_af > thr$ba1),
    BS1 = cov_ok & (gnomad_popmax_af >= thr$bs1_lo & gnomad_popmax_af < thr$bs1_hi),
    BP4 = revel <= thr$bp4_revel
  )

# --- Bayesian combine ---
classify_from_post <- function(p) {
  case_when(
    p >= 0.99  ~ "Pathogenic",
    p >= 0.90  ~ "Likely pathogenic",
    p <= 0.001 ~ "Benign",
    p <  0.10  ~ "Likely benign",
    TRUE       ~ "VUS"
  )
}

criteria_vec <- c("PS1", "PS2", "PM1", "PM2", "PP3", "PP2", "BA1", "BS1", "BP4")

res <- v2 %>%
  rowwise() %>%
  mutate(
    log_odds = log(prior_odds) +
      as.numeric(PS1) * log(LR$PS) +
      as.numeric(PS2) * log(LR$PS) +
      as.numeric(PM1) * log(LR$PM) +
      as.numeric(PM2) * log(LR$PM) +
      as.numeric(PP3) * log(LR$PP) +
      as.numeric(PP2) * log(LR$PP) +
      as.numeric(BS1) * log(LR$BS) +
      as.numeric(BP4) * log(LR$BP) +
      as.numeric(BA1) * log(LR$VSb),
    odds_post  = exp(log_odds),
    post_prob  = odds_post / (1 + odds_post),
    acmg_class = classify_from_post(post_prob),
    criteria_met = {
      flags <- c_across(all_of(criteria_vec))
      labs  <- criteria_vec[which(!is.na(flags) & flags)]
      if (length(labs) == 0) "None" else paste(labs, collapse = ";")
    }
  ) %>%
  ungroup()

res %>%
  select(everything(), post_prob, acmg_class, criteria_met) %>%
  write_tsv(file = output)
