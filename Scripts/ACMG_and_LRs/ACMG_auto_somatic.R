# ACMG_auto_somatic.R — missense-only, dominant gene ACMG classification with somatic weighting
#
# Extends ACMG_auto.R by adding somatic evidence codes from the input variant file:
#   pathogenic_somatic_strong    → PS-level evidence
#   pathogenic_somatic_moderate  → PM-level evidence
#   pathogenic_somatic_supporting→ PP-level evidence
#   benign_somatic_supporting    → BP-level evidence
#
# Usage:
#   Rscript ACMG_auto_somatic.R <variants.tsv> <hotspots.tsv> <output.tsv>
#
# Arguments:
#   variants.tsv   — annotated variant TSV with somatic evidence columns
#                    (see Prep_data_for_Auto_ACMG.Rmd for format)
#   hotspots.tsv   — two-column TSV with gene and codon columns
#   output.tsv     — path to write classified variants

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript ACMG_auto_somatic.R <variants.tsv> <hotspots.tsv> <output.tsv>")
variants_path <- args[1]
hotspots_path <- args[2]
output        <- args[3]

# RASopathy VCEP genes where PP2 is applicable (missense z score >3.09 in gnomAD)
# BRAF: GN049, MAP2K1: GN045, PTPN11: GN043 — all others explicitly not applicable
pp2_genes <- c("BRAF", "MAP2K1", "PTPN11")

# --- Parameters: prior & LRs (Tavtigian odds-path) ---
prior_p    <- 0.10
prior_odds <- prior_p / (1 - prior_p)
LR <- list(PS = 18.7, PM = 4.33, PP = 2.08,
           VSb = 1/350, BM = 1/4.33, BS = 1/18.7, BP = 1/2.08)

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
  mutate(PP2_supporting = grepl("missense", consequence) & gene %in% pp2_genes)

# --- Evidence flags (dominant, missense) ---
v2 <- v %>%
  mutate(
    PS1 = coalesce(clinvar_exact_path_star2plus, FALSE),
    PS2 = coalesce(de_novo_confirmed, FALSE),
    PM1 = hotspot,
    PP2 = coalesce(PP2_supporting, FALSE),
    # PM2/PP3/BA1/BS1/BP4 are coalesced to FALSE: if the required gnomAD AF or
    # REVEL value is missing, the evidence code is simply not applied (the
    # posterior is left unchanged) rather than producing NA and forcing a VUS.
    # PM2 ("absent from gnomAD") requires a queried AF of 0; a missing AF does
    # not count as absent.
    PM2 = coalesce(gnomad_popmax_af == 0, FALSE),
    PP3 = coalesce(revel >= thr$pp3_revel, FALSE),
    BA1 = coalesce(cov_ok & (gnomad_popmax_af > thr$ba1), FALSE),
    BS1 = coalesce(cov_ok & (gnomad_popmax_af >= thr$bs1_lo & gnomad_popmax_af < thr$bs1_hi), FALSE),
    BP4 = coalesce(revel <= thr$bp4_revel, FALSE),
    # somatic evidence flags (coalesce to FALSE so NAs do not propagate)
    SOM_PS = coalesce(as.logical(pathogenic_somatic_strong),     FALSE),
    SOM_PM = coalesce(as.logical(pathogenic_somatic_moderate),   FALSE),
    SOM_PP = coalesce(as.logical(pathogenic_somatic_supporting), FALSE),
    SOM_BP = coalesce(as.logical(benign_somatic_supporting),     FALSE)
  ) %>%
  # Apply only the STRONGEST pathogenic somatic tier so a single variant can
  # never receive more than one pathogenic somatic evidence weight. The
  # COSMIC-count thresholds are nested (strong implies moderate implies
  # supporting), so without this a strong variant could be double-counted.
  # (SOM_BP is benign evidence at count == 0 and cannot co-occur with the
  # pathogenic tiers, so it is left untouched.)
  mutate(
    SOM_PM = SOM_PM & !SOM_PS,
    SOM_PP = SOM_PP & !SOM_PS & !SOM_PM
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

criteria_vec <- c("PS1", "PS2", "PM1", "PM2", "PP3", "PP2", "BA1", "BS1", "BP4",
                  "SOM_PS", "SOM_PM", "SOM_PP", "SOM_BP")

res <- v2 %>%
  rowwise() %>%
  mutate(
    log_odds = log(prior_odds) +
      as.numeric(PS1)    * log(LR$PS)  +
      as.numeric(PS2)    * log(LR$PM)  +  # moderate strength per VCEP
      as.numeric(PM1)    * log(LR$PM)  +
      as.numeric(PM2)    * log(LR$PP)  +  # supporting strength per RASopathy VCEP
      as.numeric(PP3)    * log(LR$PP)  +
      as.numeric(PP2)    * log(LR$PP)  +
      as.numeric(BS1)    * log(LR$BS)  +
      as.numeric(BP4)    * log(LR$BP)  +
      as.numeric(BA1)    * log(LR$VSb) +
      as.numeric(SOM_PS) * log(LR$PS)  +
      as.numeric(SOM_PM) * log(LR$PM)  +
      as.numeric(SOM_PP) * log(LR$PP)  +
      as.numeric(SOM_BP) * log(LR$BP),
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
