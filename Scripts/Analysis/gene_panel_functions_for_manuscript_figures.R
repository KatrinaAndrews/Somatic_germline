# gene_panel_functions.R (edited: compact labels + option for shared labels across multiple genes)
#
# Key features:
# - Each gene plot is a 5-row patchwork (ClinVar / De novo / COSMIC / Sperm / Buccal)
# - Panel titles above each subplot are removed (less "wordy")
# - You can choose panel labels placed on the LEFT (horizontal), or no panel labels
# - For multi-gene side-by-side layouts, use `combine_gene_panel_plots_shared_labels()`
#   to show the left labels ONCE for the whole figure (recommended).

# 1) libraries & theme
library(tidyverse)
library(patchwork)
library(grid)

theme_set(theme_bw(14))

# 2) palettes
safe_colorblind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733",
  "#332288", "#AA4499", "#44AA99", "#999933",
  "#882255", "#661100", "#6699CC", "#888888"
)
missense_colors <- c("missense" = "#882255", "non_missense" = "#6699CC")

# Build a single, consistent legend (always shows both missense / non-missense)
missense_legend_grob <- function(
  title = "",
  labels = c("Missense", "Other"),
  position = "bottom"
) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required to build the shared legend. Please install it.")
  }
  dummy <- tibble::tibble(
    x = 1:2,
    y = 1,
    missense = factor(c("missense", "non_missense"), levels = c("missense", "non_missense"))
  )

  p <- ggplot2::ggplot(dummy, ggplot2::aes(x = x, y = y, fill = missense)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(
      values = missense_colors,
      breaks = c("missense", "non_missense"),
      drop = FALSE,
      labels = labels
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = title, nrow = 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = position)

  cowplot::get_legend(p)
}

# 3) helper: nice codon axis
xaxis_codons <- function(max_codon) {
  if (is.na(max_codon) || !is.finite(max_codon) || max_codon <= 0) max_codon <- 1
  if (max_codon <= 600) {
    scale_x_continuous(limits = c(0, max_codon + 20), breaks = seq(0, max_codon + 20, by = 100))
  } else if (max_codon <= 1200) {
    scale_x_continuous(limits = c(0, max_codon + 50), breaks = seq(0, max_codon + 50, by = 250))
  } else {
    scale_x_continuous(limits = c(0, max_codon + 100), breaks = pretty(c(0, max_codon)))
  }
}

# 4) helper: left-side horizontal label as a strip next to a plot
left_labeled_panel <- function(p, label, strip_width = 0.26, label_size = 0.9) {
  strip <- wrap_elements(
    full = grid::textGrob(
      label,
      rot = 0,
      x = unit(0.02, "npc"),
      just = c("left", "centre"),
      gp = grid::gpar(fontface = "bold", cex = label_size)
    )
  )
  (strip | p) + plot_layout(widths = c(strip_width, 1))
}

# 5) helper: a 5-row strip of labels for a whole multi-gene figure
panel_label_strip <- function(labels, strip_width = 0.18, label_size = 0.9) {
  stopifnot(length(labels) == 5)
  make_one <- function(lbl) {
    wrap_elements(
      full = grid::textGrob(
        lbl,
        rot = 0,
        x = unit(0.02, "npc"),
        just = c("left", "centre"),
        gp = grid::gpar(fontface = "bold", cex = label_size)
      )
    )
  }
  strip <- make_one(labels[1]) / make_one(labels[2]) / make_one(labels[3]) / make_one(labels[4]) / make_one(labels[5])
  # Note: width is controlled when we combine with `|` via plot_layout(widths=...)
  strip
}

# 6) main function: generate a single-gene 5-panel patchwork
#    panel_labels:
#      - "left" : add left strip label to EACH panel (useful when plotting one gene)
#      - "none" : no panel labels (recommended when combining multiple genes; add shared labels once)
generate_gene_panel_plots <- function(
  gene_symbol,
  clinvar,
  COSMIC,
  DNMs,
  sperm,
  buccal,
  panel_labels = c("left", "none")
) {

  panel_labels <- match.arg(panel_labels)

  clinvar_gene <- clinvar %>% filter(gene == gene_symbol)
  COSMIC_gene  <- COSMIC  %>% filter(gene == gene_symbol)
  DNMs_gene    <- DNMs    %>% filter(gene == gene_symbol)
  Sperm_gene   <- sperm   %>% filter(gene == gene_symbol)
  buccal_gene  <- buccal  %>% filter(gene == gene_symbol)

  max_codon <- suppressWarnings(max(c(
    ifelse(nrow(clinvar_gene) > 0, max(clinvar_gene$codon, na.rm = TRUE), -Inf),
    ifelse(nrow(COSMIC_gene)  > 0, max(COSMIC_gene$codon,  na.rm = TRUE), -Inf),
    ifelse(nrow(DNMs_gene)    > 0, max(DNMs_gene$codon,    na.rm = TRUE), -Inf),
    ifelse(nrow(Sperm_gene)   > 0, max(Sperm_gene$codon,   na.rm = TRUE), -Inf),
    ifelse(nrow(buccal_gene)  > 0, max(buccal_gene$codon,  na.rm = TRUE), -Inf)
  ), na.rm = TRUE))
  if (!is.finite(max_codon) || is.na(max_codon) || max_codon <= 0) max_codon <- 1

  base_panel_theme <- theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank()  # prevent repeating labels on y-axis title
  )

  # ClinVar
  clinvar_plot <- clinvar_gene %>%
    mutate(missense = factor(missense, levels = c("missense", "non_missense"))) %>%
    ggplot(aes(x = codon, fill = missense, colour = missense, y = NumberSubmitters)) +
    geom_col(width = 1, linewidth = 0.25, na.rm = TRUE) +
    xlab("") +
    xaxis_codons(max_codon) +
    scale_fill_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE) +
    scale_colour_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE, guide = "none") +
    base_panel_theme

  # De novo
  dnm_plot <- DNMs_gene %>%
    mutate(missense = factor(missense, levels = c("missense", "non_missense"))) %>%
    ggplot(aes(x = codon, fill = missense, colour = missense, y = n_DNM)) +
    geom_col(width = 1, linewidth = 0.25, na.rm = TRUE) +
    xlab("") +
    xaxis_codons(max_codon) +
    scale_fill_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE) +
    scale_colour_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE, guide = "none") +
    base_panel_theme

  # COSMIC
  cosmic_plot <- COSMIC_gene %>%
    mutate(missense = factor(missense, levels = c("missense", "non_missense"))) %>%
    ggplot(aes(x = codon, fill = missense, colour = missense, y = n_COSMIC)) +
    geom_col(width = 1, linewidth = 0.25, na.rm = TRUE) +
    xlab("") +
    xaxis_codons(max_codon) +
    scale_fill_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE) +
    scale_colour_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE, guide = "none") +
    base_panel_theme

  # Sperm
  sperm_plot <- Sperm_gene %>%
    mutate(missense = factor(missense, levels = c("missense", "non_missense"))) %>%
    ggplot(aes(x = codon, fill = missense, colour = missense, y = n_sperm)) +
    geom_col(width = 1, linewidth = 0.25, na.rm = TRUE) +
    xlab("") +
    xaxis_codons(max_codon) +
    scale_fill_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE) +
    scale_colour_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE, guide = "none") +
    base_panel_theme

  # Buccal
  buccal_plot <- buccal_gene %>%
    mutate(missense = factor(missense, levels = c("missense", "non_missense"))) %>%
    ggplot(aes(x = codon, fill = missense, colour = missense, y = n_buccal)) +
    geom_col(width = 1, linewidth = 0.25, na.rm = TRUE) +
    xlab("amino acid codon") +
    xaxis_codons(max_codon) +
    scale_fill_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE) +
    scale_colour_manual(values = missense_colors, breaks = c("missense", "non_missense"), drop = FALSE, guide = "none") +
    base_panel_theme

  # Desired labels (used either per-panel or as a shared strip)
  labels <- c(
    "n Clinvar P/LP",
    "n DNM",
    "n COSMIC",
    "n Sperm nanoseq",
    "n Buccal nanoseq"
  )

  if (panel_labels == "left") {
    combined_plots <-
      left_labeled_panel(clinvar_plot, labels[1]) /
      left_labeled_panel(dnm_plot,     labels[2]) /
      left_labeled_panel(cosmic_plot,  labels[3]) /
      left_labeled_panel(sperm_plot,   labels[4]) /
      left_labeled_panel(buccal_plot,  labels[5])
  } else {
    combined_plots <- clinvar_plot / dnm_plot / cosmic_plot / sperm_plot / buccal_plot
  }

  gene_title <- wrap_elements(
    full = grid::textGrob(
      gene_symbol,
      x = unit(0, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontface = "bold", cex = 1.25)
    )
  )

  gene_title / combined_plots + plot_layout(heights = c(0.08, 1))
}

# 7) Combine multiple gene plots with labels shown ONCE (recommended)
# Usage:
#   p1 <- generate_gene_panel_plots("PTPN11", ..., panel_labels="none")
#   p2 <- generate_gene_panel_plots("SMAD4",  ..., panel_labels="none")
#   p3 <- generate_gene_panel_plots("PIK3CA", ..., panel_labels="none")
#   combine_gene_panel_plots_shared_labels(list(p1,p2,p3), ncol=3)
combine_gene_panel_plots_shared_labels <- function(
  gene_plots,
  ncol = length(gene_plots),
  strip_width = 0.08,
  label_size = 0.9,
  show_legend = TRUE,
  legend_height = 0.10,
  legend_title = "",
  legend_labels = c("Missense", "Other"),
  labels = c(
    "n Clinvar P/LP",
    "n DNM",
    "n COSMIC",
    "n Sperm nanoseq",
    "n Buccal nanoseq"
  )
) {
  stopifnot(length(labels) == 5)
  # Make the shared strip (5 rows)
  strip <- panel_label_strip(labels, strip_width = strip_width, label_size = label_size)

  # Arrange genes side-by-side
  genes <- wrap_plots(gene_plots, ncol = ncol)

  # Put strip on the left of the entire multi-gene figure
  top <- (strip | genes) + plot_layout(widths = c(strip_width, 1))

  if (!isTRUE(show_legend)) return(top)

  leg <- missense_legend_grob(title = legend_title, labels = legend_labels, position = "bottom")
  top / patchwork::wrap_elements(full = leg) +
    patchwork::plot_layout(heights = c(1, legend_height))
}
