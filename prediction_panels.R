# =============================================================================
# PLASMA SPECIFICITY SCANNING (PSS) — STEP 4: PREDICT PLASMA SPECIFICITY
# =============================================================================
# This script implements step 4 of the PSS method. Using the optimal virus
# sub-panels selected in steps 1-3 (create_virus_panels.R), it predicts the
# neutralization specificity of each plasma/serum donor by the Spearman
# fingerprinting method and combines the predictions across all sub-panels.
#
# Two complementary outputs are produced:
#   (a) Individual bnAb-based predictions: the distribution of plasma-vs-bnAb
#       correlations across sub-panels, highlighting bnAbs with mean Spearman
#       correlation > 0.40 (see `cor_distributions`).
#   (b) bnAb class-based predictions: the number/percentage of sub-panel
#       predictions per epitope class for each donor (see `epitope_counts`).
#
# INPUT DATA FORMAT
#   `donor_data` (serum/plasma neutralization), one row per donor:
#       - column 1: unique donor ID (e.g. a combination of donor ID and date)
#       - remaining columns: 1/NT50 measurements against the virus panel
#         IMPORTANT: invert the NT50 measurements (use 1/NT50).
#       - virus column names MUST match those used in `mab_data`.
#   `mab_data` (reference bnAb neutralization): same file used in steps 1-3.
#
# Runtime: ~1-2 min.
#
# This file only DEFINES the function. Use `run_PSS.R` to format the data and
# run the full workflow.
# =============================================================================

# ---- Required packages (loaded by run_PSS.R) --------------------------------
# dplyr, tidyr, parallel
library(dplyr)
library(tidyr)
library(parallel)


# =============================================================================
# MAIN FUNCTION
# =============================================================================
# Arguments:
#   donor_data      : serum/plasma neutralization data (see format above)
#   optimal_panels  : list of optimal virus sub-panels from create_virus_panels.R
#   mab_data        : reference bnAb neutralization data
#   n_cores         : number of cores for the parallel correlation step
#   cor_min         : minimum mean correlation to call a confident prediction
#   cor_min2        : tolerance subtracted from cor_min when discarding weak hits
prediction_panels <- function(donor_data, optimal_panels, mab_data,
                               n_cores = 1L, cor_min = 0.4, cor_min2 = 0) {

  # The first column of donor_data is the donor identifier
  colnames(donor_data)[1] <- "uniqueID"

  # Lookup linking each donor ID to its row index (useful for downstream plots)
  id_index <- as.data.frame(cbind(donor_data$uniqueID, seq_len(nrow(donor_data))))
  colnames(id_index) <- c("uniqueID", "index")

  n_panels <- length(optimal_panels)

  # ---- Spearman correlation: each donor vs each reference bnAb --------------
  #      computed separately on every optimal virus sub-panel, in parallel.

  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(
    cl,
    c("mab_data", "donor_data", "optimal_panels"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, library(dplyr))

  cor_list <- parallel::parLapply(
    cl,
    seq_len(n_panels),
    function(k) {
      panel_viruses <- optimal_panels[[k]]
      cor(
        mab_data   %>% dplyr::select(all_of(panel_viruses)) %>% t(),
        donor_data %>% dplyr::select(all_of(panel_viruses)) %>% t(),
        method = "spearman",
        use    = "pairwise.complete.obs"
      )
    }
  )
  parallel::stopCluster(cl)

  # Stack the per-panel correlation matrices into a 3-D array:
  # dimensions are bnAb x donor x sub-panel. Missing correlations become 0.
  cor_array <- cor_list %>%
    unlist() %>%
    array(dim = c(dim(cor_list[[1]]), n_panels))
  cor_array[is.na(cor_array)] <- 0

  # ---- For each (donor, sub-panel): rank reference bnAbs by correlation -----
  best_mab_idx        <- apply(cor_array, c(2, 3), which.max)
  second_best_mab_idx <- apply(cor_array, c(2, 3), function(x) order(x, decreasing = TRUE)[2])
  best_cor            <- apply(cor_array, c(2, 3), max, na.rm = TRUE)

  # Discard the best match when it is too weakly correlated
  best_mab_idx[best_cor < cor_min - cor_min2] <- NA

  # Translate bnAb indices into their epitope classes
  best_epitope        <- matrix(mab_data$Epitope[best_mab_idx], ncol = ncol(best_mab_idx))
  second_best_epitope <- array(mab_data$Epitope[second_best_mab_idx], dim = dim(second_best_mab_idx))

  # If the top match is weakly correlated AND disagrees with the runner-up,
  # treat the prediction as undefined (NA).
  best_epitope[best_epitope != second_best_epitope & best_cor < cor_min] <- NA

  # ---- (b) bnAb class-based prediction: epitope counts per donor ------------
  epitope_classes <- sort(unique(mab_data$Epitope))

  # For each donor, count how many sub-panels predicted each epitope class
  epitope_counts <- apply(best_epitope, 1, function(row) {
    observed <- table(row)
    counts   <- setNames(rep(0, length(epitope_classes)), epitope_classes)
    counts[names(observed)] <- as.integer(observed)
    counts
  }) %>%
    t() %>%
    as_tibble() %>%
    mutate(uniqueID = donor_data$uniqueID)

  # Sub-panels with no confident prediction count toward an "unknown" class
  epitope_counts$unknown <-
    n_panels - apply(epitope_counts[, seq_along(epitope_classes)], 1, sum)

  # ---- Average correlation per donor across all sub-panels ------------------
  mean_cor <- apply(cor_array, c(2, 1), mean) %>% as.data.frame()
  names(mean_cor) <- mab_data$mAb
  mean_cor$ID <- donor_data$uniqueID

  # ---- (a) Individual bnAb-based prediction: distribution per donor ---------
  cor_distributions <- vector("list", nrow(donor_data))

  for (i in seq_len(nrow(donor_data))) {
    # Correlations for donor i: sub-panels (rows) x bnAbs (columns)
    panel_cor <- cor_array[, i, ] %>% t() %>% as.data.frame()
    names(panel_cor) <- mab_data$mAb
    panel_cor$Panel  <- seq_len(n_panels)
    n_col <- ncol(panel_cor)

    # Long format: one row per (sub-panel, bnAb) correlation
    cor_long <- panel_cor %>%
      pivot_longer(colnames(panel_cor)[-n_col], names_to = "mAb", values_to = "COR")
    cor_long <- merge(cor_long,
                      mab_data %>% dplyr::select(mAb, Epitope),
                      by = "mAb")

    # Flag bnAbs whose mean correlation with this donor exceeds 0.4
    strong_mabs <- cor_long %>%
      group_by(mAb) %>%
      summarise(meanCOR = mean(COR)) %>%
      filter(meanCOR > 0.4)

    cor_long$Epitopecolor <- cor_long$Epitope
    cor_long$Epitopecolor[!cor_long$mAb %in% strong_mabs$mAb] <- "under threshold"
    cor_distributions[[i]] <- cor_long
  }

  # ---- Return ---------------------------------------------------------------
  list(epitope_counts    = epitope_counts,    # (b) counts per epitope class
       cor_array         = cor_array,         # full bnAb x donor x panel array
       mean_cor          = mean_cor,          # mean correlation per donor/bnAb
       cor_distributions = cor_distributions, # (a) per-donor distributions
       id_index          = id_index)          # donor ID <-> row index lookup
}

# This script only DEFINES `prediction_panels()`.
# To format the input data and execute the full PSS workflow, run `run_PSS.R`.
