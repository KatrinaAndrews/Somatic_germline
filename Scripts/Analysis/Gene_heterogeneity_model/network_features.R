# network_features.R
# -----------------------------------------------------------------------------
# Helpers for the gene-level heterogeneity model: build a STRING
# protein-protein interaction (PPI) graph and compute
# random-walk-with-restart (RWR) proximity of every gene to a seed set.
#
# Implemented with igraph.
#
# The network-proximity FEATURE used in the model is the leakage-safe CONTRAST
#     log( RWR proximity to the "successful" (positive) network /
#          RWR proximity to the "unsuccessful" (negative) network )
# computed FOLD-AWARE (seeds = training-fold genes only). The contrast cancels
# node-degree bias (well-connected genes would otherwise score highly on
# every network, regardless of seed set).
#
# STRING download From
# https://string-db.org/cgi/download (organism = Homo sapiens, 9606):
#   - 9606.protein.links.v12.0.txt.gz   (protein1 protein2 combined_score)
#   - 9606.protein.info.v12.0.txt.gz    (#string_protein_id  preferred_name ...)
# Point STRING_LINKS_FILE / STRING_INFO_FILE (below, or via args) at these.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(igraph)
})

# ---- default file locations (override by passing args to load_string_graph) ----
STRING_LINKS_FILE <- "../../../Data/STRING/9606.protein.links.v12.0.txt.gz"
STRING_INFO_FILE  <- "../../../Data/STRING/9606.protein.info.v12.0.txt.gz"

#' Load a STRING PPI graph keyed by gene symbol.
#'
#' @param threshold Minimum STRING combined_score to keep an edge (700 = high
#'   confidence, matching the reference POC).
#' @param links_file,info_file STRING download paths.
#' @return An undirected, simplified igraph whose vertex names are gene symbols.
load_string_graph <- function(threshold = 700,
                              links_file = STRING_LINKS_FILE,
                              info_file  = STRING_INFO_FILE) {
  # Fail loudly with a helpful message if the (large, not-checked-in) STRING
  # downloads haven't been staged in Data/STRING/ yet, rather than a cryptic
  # "file not found" error from read_delim further down.
  if (!file.exists(links_file))
    stop("STRING links file not found: ", links_file,
         "\nStage 9606.protein.links.v12.0.txt.gz (see header).")
  if (!file.exists(info_file))
    stop("STRING info file not found: ", info_file,
         "\nStage 9606.protein.info.v12.0.txt.gz (see header).")

  # STRING's edge list (links_file) identifies proteins by an internal
  # "ENSP"-style ID (e.g. 9606.ENSP00000288602), not by gene symbol. The info
  # file is STRING's own lookup table mapping each of those IDs to a human
  # readable gene symbol ("preferred_name"). We build a named-vector lookup
  # (id2sym) so we can translate the edge list in one vectorised step below.
  info <- read_delim(info_file, delim = "\t", show_col_types = FALSE,
                     progress = FALSE)
  # first column name starts with '#'; normalise
  names(info)[1] <- "string_protein_id"
  id2sym <- setNames(info$preferred_name, info$string_protein_id)

  # Read every protein-protein edge, keep only "high-confidence" edges
  # (STRING's own combined_score, 0-1000 scale; 700 is STRING's documented
  # high-confidence cutoff), then translate both endpoints from ENSP IDs to
  # gene symbols. A handful of edges won't translate (id2sym lookup returns
  # NA for any ID missing from the info file) or will map both endpoints to
  # the same symbol (self-loops after translation); both are dropped.
  links <- read_delim(links_file, delim = " ", show_col_types = FALSE,
                      progress = FALSE) %>%
    filter(combined_score >= threshold) %>%
    transmute(g1 = id2sym[protein1],
              g2 = id2sym[protein2]) %>%
    filter(!is.na(g1), !is.na(g2), g1 != g2)

  # Build an undirected graph from the (gene, gene) edge list. STRING scores
  # are symmetric (no directionality to "protein A interacts with protein B"),
  # so an undirected graph is the natural representation for RWR/PageRank
  # below. simplify() removes any duplicate edges and self-loops that slipped
  # through (STRING can list the same pair from multiple evidence channels).
  g <- graph_from_data_frame(links[, c("g1", "g2")], directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  g
}

#' Random-walk-with-restart proximity from a seed set (personalised PageRank).
#'
#' Restart probability `restart` = probability of teleporting back to the seeds
#' at each step (higher = more local to seeds). This maps onto igraph's
#' personalised PageRank with damping = 1 - restart and a personalisation vector
#' that is uniform over the seed genes.
#'
#' @param graph igraph from load_string_graph().
#' @param seeds Character vector of seed gene symbols (only those present in the
#'   graph are used).
#' @param restart Restart probability in (0,1); ~0.5 in the reference POC.
#' @return Named numeric vector of proximity scores over ALL graph vertices.
#'   NB: seed genes carry an inflated self-score, so only NON-seed scores are
#'   valid for evaluation (the model uses this fold-aware — seeds are training
#'   genes only, scores are read for held-out / grey-zone genes).
rwr_proximity <- function(graph, seeds, restart = 0.5) {
  # Not every seed gene symbol will necessarily be a vertex in the STRING
  # graph (STRING doesn't cover every gene, and gene-symbol spelling can
  # differ). Restrict to seeds that are actually present; if literally none
  # of the seeds are in the graph, there's nothing to walk from, so return an
  # all-NA vector rather than letting page_rank() error out.
  seeds_in <- intersect(seeds, V(graph)$name)
  if (length(seeds_in) == 0)
    return(setNames(rep(NA_real_, vcount(graph)), V(graph)$name))

  # Personalised PageRank needs a "personalisation vector": one weight per
  # vertex in the graph, specifying where a random walk restarts to. We want
  # restarts to land uniformly on the seed genes and nowhere else, so this is
  # a 0/1 indicator (1 for each seed vertex, 0 elsewhere) renormalised to sum
  # to 1 (a probability distribution over vertices).
  pvec <- as.numeric(V(graph)$name %in% seeds_in)
  pvec <- pvec / sum(pvec)

  # igraph's page_rank() implements the walk: at each step, follow a random
  # edge from the current vertex with probability `damping`, or "teleport"
  # back to a seed (drawn from pvec) with probability `1 - damping`. We
  # parameterise this in terms of `restart` (the probability of teleporting
  # back to a seed at each step) rather than igraph's own `damping`, since
  # `restart` is the more natural quantity for an RWR-proximity feature:
  # restart = 1 - damping, so a higher restart means the walk stays "more
  # local" to the seeds and diffuses less far into the rest of the network.
  pr <- page_rank(graph, algo = "prpack", directed = FALSE,
                  personalized = pvec, damping = 1 - restart)$vector

  # The result is one proximity score per vertex in the whole graph (not just
  # the seeds) -- this is exactly what we want, since we'll later read off
  # scores for query genes that are NOT seeds themselves.
  setNames(pr, V(graph)$name)
}

#' Fold-aware network-proximity CONTRAST feature for a set of query genes.
#'
#' Seeds one RWR on the positive (successful) genes and one on the negative
#' (unsuccessful) genes, then returns log10((prox_pos + eps)/(prox_neg + eps))
#' for each query gene. Call with TRAINING-fold positives/negatives as seeds and
#' the held-out (or grey-zone) genes as queries to keep it leakage-safe.
#'
#' @param graph igraph.
#' @param pos_seeds,neg_seeds Seed gene-symbol vectors (training-fold only).
#' @param query_genes Genes to score (held-out / grey-zone).
#' @param restart RWR restart probability.
#' @return Named numeric vector of the log-ratio contrast for query_genes.
#'   Genes absent from the graph return NA (median-imputed downstream).
network_contrast <- function(graph, pos_seeds, neg_seeds, query_genes,
                             restart = 0.5) {
  # Run TWO independent random walks: one seeded only on the positive
  # (high-utility) genes, one seeded only on the negative (low-utility)
  # genes. Every vertex in the graph now has two proximity scores: "how
  # close am I to the positive network" and "how close am I to the negative
  # network".
  prox_pos <- rwr_proximity(graph, pos_seeds, restart = restart)
  prox_neg <- rwr_proximity(graph, neg_seeds, restart = restart)

  # PageRank scores are always > 0 for reachable vertices (in a connected
  # graph, restart alone gives every vertex a small positive floor) but can
  # be extremely close to zero for vertices far from the seeds. Taking
  # log10() of a value that's ever exactly 0 or negative would produce -Inf
  # or NaN, so we add a small positive constant (`eps`) to both proximities
  # before taking logs. `eps` is set relative to the data itself (a tenth of
  # the smallest positive proximity actually observed on either walk) rather
  # than an arbitrary fixed number, so it stays sensible regardless of graph
  # size or restart probability.
  pos_pos <- prox_pos[prox_pos > 0]
  pos_neg <- prox_neg[prox_neg > 0]
  eps <- min(c(pos_pos, pos_neg), na.rm = TRUE) / 10

  # Pre-fill the output with NA for every query gene, then overwrite only
  # the ones that are actually vertices in the graph. Query genes absent
  # from STRING (e.g. very obscure genes) are left as NA here and are
  # median-imputed later downstream, with an explicit missingness flag.
  out <- rep(NA_real_, length(query_genes))
  names(out) <- query_genes
  in_g <- query_genes[query_genes %in% V(graph)$name]

  # The feature itself: log10(proximity-to-positive-network) minus
  # log10(proximity-to-negative-network), i.e. log10 of the RATIO of the two
  # proximities. Positive values mean "closer, in the network, to the
  # high-utility gene set than to the low-utility gene set"; negative values
  # mean the reverse. Working in log-ratio space (rather than a plain
  # difference or ratio) makes the feature symmetric around 0 and keeps it
  # on a scale comparable to the other model features.
  out[in_g] <- log10(prox_pos[in_g] + eps) - log10(prox_neg[in_g] + eps)
  out
}
