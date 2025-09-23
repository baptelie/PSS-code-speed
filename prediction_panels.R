# ### CODE FOR CREATING THE VIRUS PANELS PREDICTION
# ### INPUT data should have unique ID
# ########
# 
# ### PACKAGES
# library(readxl)
# library(dplyr)
# library(magrittr)
# library(tidyr)
# library(openxlsx)
# 
# ### DIRECTORY
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 
# 
# # LOAD THE SERUM NEUTRALIZATION DATA
# # NEEDS TO BE IN FOLLOWING FORMAT FOR RUNNING THE create_panels FUNCTION:
# # FIRST COLUMN: uniqueID = donor ID (can be a combination of ID and DATE)
# # REST OF THE COLUMNS = 1/NT50 measurements on virus panel (ENSURE THAT YOU INVERT THE MEASUREMENTS)
# # ENSURE THAT YOU REMOVE CNE40_BC
# # VIRUS NAMES SHOULD BE THE SAME IN data.donor AS IN data.mAb.panel
# 
# data.donor <- read_xlsx("/Volumes/Research/Group-AT/Projects/1 BOOK Swiss4.5K_XBnAb/XbnAb neut database July 2025/1_Master file 304 XbnAb users NT and prediction status Multi paper submitted July 9 2025.xlsx",
#                         sheet = "XbNAb NT50 working sheet",
#                         skip = 3
# ) %>% 
#   filter(!if_all(8:48, is.na)) %>% 
#   mutate(across(8:48, ~ {
#   pmin(as.double(.), 1e6)                    # cap values at 1e6
# })) %>% 
#   select(1,8:48) %>% 
#   rename(BG505_W6M_C2_T332N = BG505_W6M_ENV_A5_T332N)
# 
# 
# 
# # LOAD OUTCOMES FROM THE create_virus_panels FUNCTION
# panel_opt = readRDS("../results/allpanels_basic_220925.RDS")
# paneltest <- read.table("../results/panelopt_basic_220925.txt")
# 
# # LOAD mAb panel data (should be the same as input for the create_virus_panels FUNCTION)
# data.mAb.panel_full <- read_xlsx("/Volumes/Research/Group-AT/Projects/1 BOOK Swiss4.5K_XBnAb/XbnAb neut database July 2025/5_Working file SHCS bnAbs with final neut data sets_July 2025 onwards.xlsx",
#                                  skip = 7
# ) %>%
#   select(-`bnAb...14`) %>%
#   rename(mAb = `bnAb...1`) %>%
#   mutate(across(14:54, ~ as.double(gsub(">25", "50", .))))
# 
# data.mAb.panel <- data.mAb.panel_full %>%
#   filter(`Reference bnAb           Set Basic (n=43)` == "x") %>%
#   select(
#     -`Epitope Region`, -`Reference/PubMed ID`, -Comments, -starts_with("Reference"),
#     -starts_with("Set Multi")
#   )
# 
# data.mAb.panel.nona <- data.mAb.panel %>% mutate(across(3:ncol(.), ~ replace(., is.na(.), 50)))


# FUNCTION FOR RUNNING THE VIRUS PANELS METHOD
# ARGUMENTS OF THE FUNCTION:
# neut data from donors
# outcomes from create_virus_panels function
# neut data from mAb

prediction_panels <- function(data.donor, panel_opt, data.mAb.panel, mcores = 1L) {
  colnames(data.donor)[1] <- "uniqueID"
  
  # ID MATCHING FOR PLOT
  id.match <- as.data.frame(cbind(data.donor$`uniqueID`, seq(1, nrow(data.donor))))
  colnames(id.match) <- c("uniqueID", "index")
  
  # PREDICTION AND CORRELATION MATRIX
  Npan <- length(panel_opt)
  
  # COMPUTE CORRELATION FOR EACH PANEL
  cor_matrix_all1 = parallel::mclapply(1:length(panel_opt),
                     function(k){
                       # name.panel.temp <- panel_opt[[k]]
                       # panel.temp <- data.mAb.panel %>% dplyr::select(mAb, all_of(name.panel.temp))
                       # data.patient.temp <- data.patient %>% dplyr::select(`uniqueID`, all_of(name.panel.temp))
                       # cor.matrix <- as.data.frame(matrix(NA, nrow = nrow(data.patient.temp), ncol = nrow(panel.temp) + 1))
                       # colnames(cor.matrix) <- c("uniqueID", panel.temp$mAb)
                       # cor.matrix$uniqueID <- data.patient$uniqueID
                       cor(data.mAb.panel %>% dplyr::select(all_of(panel_opt[[k]])) %>% t,
                           data.donor %>% dplyr::select(all_of(panel_opt[[k]])) %>% t,
                           method = "spearman",
                           use = "pairwise.complete.obs")
                     }, mc.cores = mcores)
  array_correlations1 = cor_matrix_all1 %>% unlist %>% array(dim = c(dim(cor_matrix_all1[[1]]), Npan))
  array_correlations1[is.na(array_correlations1)] = 0

  which_max_spearman = apply(array_correlations1, c(2,3), which.max)
  max_spearman = apply(array_correlations1, c(2,3), max, na.rm=T)
  which_max_spearman[max_spearman < .4] = NA
  max_spearman_epitope = data.mAb.panel$Epitope[which_max_spearman] %>% matrix(ncol = ncol(which_max_spearman))
  # Unique characters across the whole matrix
  chars <- sort(unique(data.mAb.panel$Epitope))
  
  # Count per row and build tibble
  max_spearman_epitope_counts <- apply(max_spearman_epitope, 1, function(row) {
    tab <- table(row)
    counts <- setNames(rep(0, length(chars)), chars) # start with zeros
    counts[names(tab)] <- as.integer(tab)            # fill observed counts
    counts
  }) %>%
    t() %>%
    as_tibble() %>%
    mutate(uniqueID = data.donor$uniqueID)
  
  max_spearman_epitope_counts$unknown = Npan - apply(max_spearman_epitope_counts[,1:length(chars)],1, sum)
  
  # max_spearman_epitope_counts %<>% 
  #   mutate(
  #     CD4bs_pct = CD4bs / Npan,
  #     IFFP_pct = `Interface/Fusion Peptide` / Npan,
  #     MPER_pct = MPER / Npan,
  #     V1V2_pct = `V2-Apex` / Npan,
  #     V3Gly_pct = `V3-Glycan` / Npan,
  #     V3GlyT2_pct = ifelse(any(chars == "Inter-V3"), `Inter-V3` / Npan, 0),
  #     SF_pct = `Silent Face` / Npan,
  #     nopred_pct = unknown / Npan
  #   )
  
  # CORRELATION AVERAGE VALUES
  data.average.cor1 = apply(array_correlations1, c(2,1), mean)
  data.average.cor1 %<>% as.data.frame()
  names(data.average.cor1) = c(data.mAb.panel$mAb)
  data.average.cor1$ID = data.donor$uniqueID
  
  # CORRELATION DISTRIBUTION PER donor
  data_long <- data.average.cor1 %>%
    pivot_longer(1:(ncol(.)-1)) %>%
    as.data.frame()
  
  cor.data <- vector(mode = "list", length = nrow(data.donor))
  
  for (i in 1:nrow(data.donor)) {
    data.temp = array_correlations1[,i,] %>% t %>% as.data.frame()
    data.temp$Panel <- seq(1:Npan)
    NCOL <- ncol(data.temp)
    
    ## COR DISTRIBUTION
    data.temp.temp <- data.temp %>%
      pivot_longer(colnames(data.temp)[-NCOL], names_to = "mAb", values_to = "COR")
    data.temp.temp <- merge(data.temp.temp,
                            data.mAb.panel %>% select(mAb, Epitope),
                            by = "mAb"
    )
    test.temp <- data.temp.temp %>%
      group_by(mAb) %>%
      summarise(meanCOR = mean(COR)) %>%
      filter(meanCOR > 0.4)
    data.temp.temp$Epitopecolor <- data.temp.temp$Epitope
    data.temp.temp$Epitopecolor[!data.temp.temp$mAb %in% test.temp$mAb] <- "under threshold"
    cor.data[[i]] <- data.temp.temp
  }
  
  ### OUTCOME
  list(max_spearman_epitope_counts, array_correlations1, data.average.cor1, cor.data, id.match)
}

# RUN THE FUNCTION
# outcome <- prediction_panels(data.donor, panel_opt, data.mAb.panel)

