# =============================================================================
# PLASMA SPECIFICITY SCANNING (PSS) — MAIN DRIVER SCRIPT
# =============================================================================
# This is the single entry point for the PSS analysis. It:
#   1. Reads the paper's supplementary tables.
#   2. Formats them into the two data frames used by the method:
#        - mab_data   (reference bnAb panel)  <- Table S3_01
#        - donor_data (plasma/serum data)      <- Table S1_02
#      and saves them as data_mAb_panel.xlsx and data_donor.xlsx.
#   3. Sources the two function scripts and runs the full workflow:
#        - build_optimal_panels()  (steps 1-3, create_virus_panels.R)
#        - prediction_panels()     (step 4,    prediction_panels.R)
#   4. Writes all outputs, automatically date-stamped with the run date.
#
# HOW TO RUN
#   - Put this script, create_virus_panels.R, prediction_panels.R and the two
#     supplementary Excel files in the same folder.
#   - Set that folder as the working directory (see below), then run:
#         source("run_PSS.R")
# =============================================================================

# ---- Packages ---------------------------------------------------------------
library(readxl)     # read the supplementary .xlsx tables
library(openxlsx)   # write .xlsx outputs
library(dplyr)
library(tidyr)
library(parallel)

# ---- Working directory ------------------------------------------------------
# Run from the folder containing all scripts and data files.
# In RStudio you can set it automatically with:
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# =============================================================================
# CONFIGURATION  (edit here if needed)
# =============================================================================
# Input supplementary tables (place them in the working directory)
file_S3 <- "../Table S3.xlsx"   # contains sheet "Table S3_01" (bnAb panel)
file_S1 <- "../Table S1.xlsx"            # contains sheet "Table S1_02" (donor data)

# Where to write outputs (created if it does not exist)
out_dir <- "../results_test"

# Number of CPU cores for the parallel steps
n_cores <- max(1L, parallel::detectCores() - 1L)

# Panel-building settings (paper: sizes 15-35, 3,000 panels/size, keep top 100)
min_panel_size    <- 15
max_panel_size    <- 35
n_panels_per_size <- 100
n_top_panels      <- 100

# Execution date, used to stamp every output file name
run_date <- format(Sys.Date(), "%Y%m%d")

# Helper to build a date-stamped output path, e.g. out_path("prediction", ".xlsx")
out_path <- function(name, ext) file.path(out_dir, paste0(name, "_", run_date, ext))

dir.create(out_dir, showWarnings = FALSE)

# =============================================================================
# STEP 0 — FORMAT THE SUPPLEMENTARY TABLES
# =============================================================================
# Both tables carry three header rows (title / measurement / subtype) above the
# real column names, so we skip the first 3 rows when reading.

# ---- Reference bnAb panel: Table S3_01 --------------------------------------
# Layout after skipping 3 rows:
#   col 1  = bnAb              -> renamed "mAb"
#   col 2  = Epitope
#   cols 3-5 = reference set / PubMed annotations (dropped)
#   cols 6-46 = 41 reference viruses (IC50, ug/mL; lower = more potent)
# Censored values ">25" mean no neutralization at the highest tested
# concentration and are set to 25.
mab_raw <- read_xlsx(file_S3, sheet = "Table S3_01", skip = 3, .name_repair = "minimal")

mab_data <- mab_raw[, c(1, 2, 6:46)]
names(mab_data)[1:2] <- c("mAb", "Epitope")
names(mab_data)      <- trimws(names(mab_data))
mab_data <- mab_data[!is.na(mab_data$mAb), ]                 # drop trailing empty rows

clean_ic50 <- function(x) suppressWarnings(as.numeric(sub("^>", "", as.character(x))))
virus_cols <- names(mab_data)[-(1:2)]
mab_data[virus_cols] <- lapply(mab_data[virus_cols], clean_ic50)

# ---- Donor (plasma/serum) data: Table S1_02 ---------------------------------
# Layout after skipping 3 rows:
#   col 1   = donor ID         -> renamed "uniqueID"
#   cols 2-6 = breadth / potency / inducer annotations (dropped)
#   cols 7-47 = 41 viruses (NT50; higher = more potent)
# The bottom rows of the sheet are footnotes (their first cell starts with a
# digit) and are removed. NT50 values are inverted to 1/NT50 so that, like the
# bnAb IC50 data, *lower = more potent*. Censoring:
#   "<40"      -> 40        (no neutralization up to the 1/40 dilution)
#   ">1000000" -> 1000000
#   "NA"/blank -> NA
donor_raw <- read_xlsx(file_S1, sheet = "Table S1_02", skip = 3, .name_repair = "minimal")

donor_data <- donor_raw[, c(1, 7:47)]
names(donor_data)[1] <- "uniqueID"
names(donor_data)    <- trimws(names(donor_data))
donor_data <- donor_data[!is.na(donor_data$uniqueID) &
                           !grepl("^[0-9]", donor_data$uniqueID), ]   # drop footnotes

clean_nt50 <- function(x) {
  x <- gsub("^[<>]", "", as.character(x))   # strip "<" / ">" censoring marks
  x[x %in% c("", "NA")] <- NA
  1 / suppressWarnings(as.numeric(x))       # invert NT50
}
donor_virus_cols <- names(donor_data)[-1]
donor_data[donor_virus_cols] <- lapply(donor_data[donor_virus_cols], clean_nt50)

# ---- Sanity checks ----------------------------------------------------------
# In the corrected tables the 41 viruses share identical names AND order between
# the bnAb panel and the donor data; we check this holds.
stopifnot(ncol(mab_data) == 43)                          # mAb + Epitope + 41 viruses
stopifnot(identical(virus_cols, setdiff(names(donor_data), "uniqueID")))
message(sprintf("Formatted %d bnAbs and %d donors over %d viruses.",
                nrow(mab_data), nrow(donor_data), length(virus_cols)))

# ---- Save the formatted inputs (for transparency / reuse) -------------------
write.xlsx(mab_data,   file.path(out_dir, "data_mAb_panel.xlsx"))
write.xlsx(donor_data, file.path(out_dir, "data_donor.xlsx"))

# =============================================================================
# STEPS 1-3 — BUILD AND SELECT THE OPTIMAL VIRUS SUB-PANELS
# =============================================================================
source("create_virus_panels.R")

panels_out <- build_optimal_panels(
  mab_data,
  min_panel_size    = min_panel_size,
  max_panel_size    = max_panel_size,
  n_panels_per_size = n_panels_per_size,
  n_cores           = n_cores,
  n_top_panels      = n_top_panels
)

all_scores     <- panels_out$all_scores
top_panels     <- panels_out$top_panels
optimal_panels <- panels_out$optimal_panels

write.table(all_scores, out_path("allpanels", ".txt"), row.names = FALSE)
write.table(top_panels, out_path("panelopt",  ".txt"), row.names = FALSE)
save(optimal_panels, file = out_path("panelopt", ".RData"))

# =============================================================================
# STEP 4 — PREDICT PLASMA SPECIFICITY
# =============================================================================
source("prediction_panels.R")

prediction_out <- prediction_panels(
  donor_data,
  optimal_panels,
  mab_data,
  n_cores = n_cores
)

epitope_counts    <- prediction_out$epitope_counts
cor_array         <- prediction_out$cor_array
mean_cor          <- prediction_out$mean_cor
cor_distributions <- prediction_out$cor_distributions
id_index          <- prediction_out$id_index

write.xlsx(epitope_counts, out_path("prediction", ".xlsx"))
write.xlsx(mean_cor,       out_path("averagecor", ".xlsx"))
write.xlsx(id_index,       out_path("ID_index",   ".xlsx"))
save(cor_array,         file = out_path("cormatrix", ".RData"))
save(cor_distributions, file = out_path("cordata",   ".RData"))

message("PSS workflow complete. Outputs written to '", out_dir, "/'.")
