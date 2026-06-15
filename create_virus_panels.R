# =============================================================================
# PLASMA SPECIFICITY SCANNING (PSS) — STEPS 1-3: BUILD OPTIMAL VIRUS SUB-PANELS
# =============================================================================
# This script implements steps 1-3 of the PSS method:
#   1. Randomly draw virus sub-panels of varying sizes from the reference panel.
#   2. Score each sub-panel by how well it recovers the known epitope clusters
#      of the reference bnAbs (default score: homogeneity, an entropy-based
#      measure of agreement between two clusterings).
#   3. Keep the best-scoring sub-panels.
#
# INPUT
#   `mab_data`: reference bnAb (monoclonal antibody) neutralization data, with
#   one row per bnAb. Required columns:
#       - mAb     : bnAb identifier
#       - Epitope : known epitope class of the bnAb
#       - <virus columns> : neutralization measurement of the bnAb against each
#                           reference virus (one column per virus)
#
# OUTPUT (a named list)
#       - all_scores     : score of every candidate sub-panel
#       - top_panels      : the `n_top_panels` best-scoring sub-panels
#       - optimal_panels : the virus names making up each top sub-panel
#
# Runtime: ~5-10 min on a MacBook Pro M1 (8 cores) with the paper's settings.
#
# This file only DEFINES the function. Use `run_PSS.R` to format the data and
# run the full workflow.
# =============================================================================

# ---- Required packages (loaded by run_PSS.R) --------------------------------
# dplyr, clevr, foreach, doParallel, parallel
library(dplyr)
library(clevr)        # entropy-based cluster comparison (homogeneity, etc.)
library(foreach)
library(doParallel)


# =============================================================================
# MAIN FUNCTION
# =============================================================================
build_optimal_panels <- function(mab_data,
                                  min_panel_size,
                                  max_panel_size,
                                  n_panels_per_size,
                                  n_cores,
                                  n_top_panels = 100,
                                  criterion    = "homogeneity",
                                  n_clusters   = n_distinct(mab_data$Epitope)) {

  # ---------------------------------------------------------------------------
  # Helper: Spearman correlation matrix between bnAbs, NA-aware.
  #
  # `panel_data` has one row per selected virus and one column per bnAb.
  # The function returns the bnAb-by-bnAb Spearman correlation matrix, OR a
  # matrix still containing NA as soon as any bnAb pair has fewer than two
  # jointly non-missing viruses (or an undefined correlation). The caller skips
  # any panel whose matrix contains NA.
  # ---------------------------------------------------------------------------
  compute_spearman_or_na <- function(panel_data) {
    # Rank within each bnAb column, keeping NAs in place
    # (Spearman correlation = Pearson correlation on the ranks).
    ranked <- apply(panel_data, 2, function(x) rank(x, na.last = "keep"))
    n_mab  <- ncol(ranked)
    cor_matrix <- matrix(NA, n_mab, n_mab)

    for (i in 1:(n_mab - 1)) {
      cor_matrix[i, i] <- 1
      for (j in (i + 1):n_mab) {
        valid <- !is.na(ranked[, i]) & !is.na(ranked[, j])
        if (sum(valid) > 1) {
          corr <- suppressWarnings(cor(ranked[valid, i], ranked[valid, j]))
          if (is.na(corr)) {
            return(cor_matrix)   # undefined correlation -> abandon this panel
          }
          cor_matrix[i, j] <- corr
          cor_matrix[j, i] <- corr
        } else {
          return(cor_matrix)     # not enough overlap -> abandon this panel
        }
      }
    }
    cor_matrix[n_mab, n_mab] <- 1
    colnames(cor_matrix) <- colnames(panel_data)
    rownames(cor_matrix) <- colnames(panel_data)
    cor_matrix
  }

  # ---- Prepare data ---------------------------------------------------------
  # Drop the mAb and Epitope columns, then transpose so that
  # rows = viruses and columns = bnAbs. Row names hold the virus names.
  neut_by_virus <- t(as.matrix(mab_data[, -c(1, 2)]))
  colnames(neut_by_virus) <- mab_data$mAb

  # Named lookup: epitope class for each bnAb
  epitope_by_mab <- setNames(mab_data$Epitope, mab_data$mAb)

  # ---- Step 1: draw candidate sub-panels ------------------------------------
  n_sizes   <- max_panel_size - min_panel_size + 1   # number of distinct sizes
  n_viruses <- ncol(mab_data) - 2                    # viruses in reference panel

  set.seed(98)   # for reproducibility
  panels_by_size <- lapply(
    seq_len(n_sizes),
    function(s) {
      panel_size <- min_panel_size - 1 + s
      # Draw `n_panels_per_size` UNIQUE sorted virus index sets of this size
      sampled <- list()
      while (length(sampled) < n_panels_per_size) {
        sampled <- unique(c(
          sampled,
          replicate(n_panels_per_size - length(sampled),
                    sort(sample.int(n_viruses, panel_size)),
                    simplify = FALSE)
        ))
      }
      sampled
    }
  )

  n_total_panels <- n_sizes * n_panels_per_size

  # Map each flat panel id (1..n_total_panels) to its (size index, replicate index).
  # Precomputing this avoids a division/modulo inside every parallel iteration.
  panel_index <- cbind(
    size_idx = ((seq_len(n_total_panels) - 1L) %/% n_panels_per_size) + 1L,
    rep_idx  = ((seq_len(n_total_panels) - 1L) %%  n_panels_per_size) + 1L
  )

  # Split the work into chunks to reduce parallel scheduling overhead
  chunk_size <- max(100L, n_total_panels %/% (n_cores * 100L))
  chunks <- split(seq_len(n_total_panels),
                  ceiling(seq_along(seq_len(n_total_panels)) / chunk_size))

  # Resolve the scoring function (e.g. clevr::homogeneity).
  # "cv" selects a custom nearest-neighbour rule handled inline below.
  if (criterion != "cv") {
    criterion_fun <- match.fun(criterion)
  } else {
    criterion_fun <- NA
  }

  # ---- Step 2: score every candidate sub-panel (in parallel) ----------------
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  parallel::clusterExport(cl, c(
    "panels_by_size", "neut_by_virus", "compute_spearman_or_na", "criterion",
    "n_clusters", "epitope_by_mab", "criterion_fun", "panel_index"
  ), envir = environment())
  doParallel::registerDoParallel(cl)

  # Each worker returns a list of c(panel_id, score) vectors; skipped panels
  # contribute NULL. Combining lists with "c" is much faster than rbind here.
  score_list <- foreach(
    chunk = chunks,
    .packages = "clevr",
    .combine  = "c"
  ) %dopar% {
    results <- vector("list", length(chunk))

    for (i in seq_along(chunk)) {
      panel_id <- chunk[i]
      size_idx <- panel_index[panel_id, 1L]
      rep_idx  <- panel_index[panel_id, 2L]

      virus_idx  <- panels_by_size[[size_idx]][[rep_idx]]
      panel_data <- neut_by_virus[virus_idx, , drop = FALSE]

      cor_matrix <- compute_spearman_or_na(panel_data)
      if (anyNA(cor_matrix)) next   # skip panels with undefined correlations

      if (criterion == "cv") {
        # Custom rule: count bnAbs whose most-correlated neighbour shares its
        # epitope class (or whose best correlation is weak, < 0.5).
        diag(cor_matrix) <- NA
        n_match <- 0
        for (mab in colnames(cor_matrix)) {
          nearest_mab <- names(which.max(cor_matrix[mab, ]))
          nearest_cor <- cor_matrix[mab, nearest_mab]
          if (epitope_by_mab[mab] == epitope_by_mab[nearest_mab] || nearest_cor < 0.5) {
            n_match <- n_match + 1
          }
        }
        results[[i]] <- c(panel_id, n_match)
      } else {
        # Cluster bnAbs by correlation distance, then compare the resulting
        # clusters to the known epitope classes using `criterion_fun`.
        dist_mat           <- as.dist(1 - cor_matrix)
        hc                 <- hclust(dist_mat)          # complete linkage (default)
        cluster_assignment <- cutree(hc, n_clusters)
        epitope_labels     <- epitope_by_mab[names(cluster_assignment)]
        results[[i]]       <- c(panel_id, criterion_fun(epitope_labels, cluster_assignment))
      }
    }
    results
  }

  parallel::stopCluster(cl)

  # ---- Step 3: assemble scores and keep the best sub-panels -----------------
  scores_df <- as.data.frame(do.call(rbind, score_list))
  colnames(scores_df) <- c("panel_id", "score")

  scores_ordered <- scores_df %>% arrange(desc(score))

  top_panels   <- head(scores_ordered, n_top_panels)
  n_top_panels <- nrow(top_panels)   # in case fewer panels than requested exist

  # Recover the virus names that make up each top sub-panel
  optimal_panels <- vector("list", n_top_panels)
  for (i in seq_len(n_top_panels)) {
    panel_id <- top_panels$panel_id[i]
    size_idx <- (panel_id - 1) %/% n_panels_per_size + 1
    rep_idx  <- panel_id - ((size_idx - 1) * n_panels_per_size)
    optimal_panels[[i]] <- rownames(neut_by_virus)[panels_by_size[[size_idx]][[rep_idx]]]
  }

  list(all_scores     = scores_df,
       top_panels     = top_panels,
       optimal_panels = optimal_panels)
}

# This script only DEFINES `build_optimal_panels()`.
# To format the input data and execute the full PSS workflow, run `run_PSS.R`.
