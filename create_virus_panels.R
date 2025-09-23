# ### CODE FOR CREATING THE VIRUS PANELS
# ### INDICATE UPDATES HERE:
# 
# ### PACKAGES
# library(readxl)
# library(dplyr)
# library(tidyr)
# library(clevr)
# library(foreach)
# library(doParallel)
# library(doSNOW)
# 
# ### DIRECTORY
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 
# # LOAD THE NEUTRALIZATION DATA FOR KNOWN mAb
# # NEEDS TO BE IN FOLLOWING FORMAT FOR RUNNING THE create_panels FUNCTION:
# # FIRST COLUMN: "mAb"
# # SECOND COLUMN: "Epitope"
# # REST OF THE COLUMN = VIRUS PANEL MEASUREMENTS (ONE COLUMN = ONE VIRUS)
# # ENSURE THAT YOU REMOVE CNE40_BC
# 
# ## mAb
# data.mAb.panel_full <- read_xlsx("/Volumes/Research/Group-AT/Projects/1 BOOK Swiss4.5K_XBnAb/XbnAb neut database July 2025/5_Working file SHCS bnAbs with final neut data sets_July 2025 onwards.xlsx",
#   skip = 7
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


create_panels <- function(data.mAb.panel, min.panel, max.panel, n.panel, n.cores, Npanel = 40, use.na.cor = "pairwise.complete.obs") {
  
  # browser()
  # Prepare data
  data.forcor <- t(as.matrix(data.mAb.panel[, -c(1, 2)]))
  colnames(data.forcor) <- data.mAb.panel$mAb
  data.panel <- data.forcor

  # Map epitopes to numeric categories once
  mAb_existingcl <- data.mAb.panel %>%
    dplyr::select(mAb, Epitope)
  epitope_map <- setNames(mAb_existingcl$Epitope, mAb_existingcl$mAb)

  # Generate candidate panels
  num.panel <- max.panel - min.panel + 1
  seq.length <- ncol(data.mAb.panel) - 2
  x <- seq_len(seq.length)

  set.seed(98)
  list_panel <- vector("list", num.panel)
  for (j in seq_len(num.panel)) {
    N <- min.panel - 1 + j
    # sample panels more efficiently
    out <- list()
    while (length(out) < n.panel) {
      out <- unique(c(
        out,
        replicate(n.panel - length(out), sort(sample.int(seq.length, N)), simplify = FALSE)
      ))
    }
    list_panel[[j]] <- out
    
    # list_panel[[j]] <- unique(replicate(n.panel, sort(sample.int(seq.length, N)), simplify = FALSE))
    # while (length(list_panel[[j]]) < n.panel) {
    #   list_panel[[j]] <- unique(c(
    #     list_panel[[j]],
    #     unique(sort(sample.int(seq.length, N)))
    #   ))
    # }
  }
  total_length <- num.panel * n.panel
  NCLUSTER <- length(unique(mAb_existingcl$Epitope))

  # Start parallel backend
  cl <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  # browser()

  homogen_vector <- foreach(
    j = seq_len(total_length),
    .combine = "rbind",
    .packages = c("clevr")
  ) %dopar% {
    nlist.temp <- (j - 1) %/% n.panel + 1
    klist.temp <- j - ((nlist.temp - 1) * n.panel)
    combvirus.temp <- list_panel[[nlist.temp]][[klist.temp]]

    data.red <- data.panel[combvirus.temp, , drop = FALSE]
    cor.matrix <- suppressWarnings(cor(data.red, method = "spearman", use = use.na.cor))

    if (!anyNA(cor.matrix)) {
      hier.cl <- hclust(as.dist(1 - cor.matrix))
      cl.attr <- cutree(hier.cl, NCLUSTER)

      # Homogeneity metrics
      epitope_vec <- epitope_map[names(cl.attr)]
      c(
        j,
        homogeneity(epitope_vec, cl.attr),
        completeness(epitope_vec, cl.attr),
        v_measure(epitope_vec, cl.attr),
        rand_index(epitope_vec, cl.attr)
      )
    }
  }

  parallel::stopCluster(cl)

  # Process results
  homogen_data <- as.data.frame(homogen_vector)
  homogen_order <- homogen_data %>%
    arrange(desc(V2), desc(V3), desc(V4))

  panel.tokeep <- head(homogen_order, Npanel)
  Npanel <- nrow(panel.tokeep)

  # Recover selected panels
  panel_opt <- vector("list", Npanel)
  for (i in seq_len(Npanel)) {
    j <- panel.tokeep$V1[i]
    nlist.temp <- (j - 1) %/% n.panel + 1
    klist.temp <- j - ((nlist.temp - 1) * n.panel)
    panel_opt[[i]] <- rownames(data.panel)[list_panel[[nlist.temp]][[klist.temp]]]
  }

  list(homogen_data, panel.tokeep, panel_opt)
}

# library(bench)
# 
# bench::mark(
#   # original = create_panels2(data.mAb.panel,min.panel=10,max.panel=35,n.panel=500),
#   optimized = create_panels(data.mAb.panel.nona, min.panel = 10, max.panel = 35, n.panel = 100, n.cores = 1L),
#   optimized_nona = create_panels(data.mAb.panel.nona, min.panel = 10, max.panel = 35, n.panel = 100, n.cores = 1L, use.na.cor = "everything"),
#   iterations = 1
# )
# 
# library(profvis)
# 
# profvis({
#   outcome <- create_panels(data.mAb.panel, min.panel = 10, max.panel = 35, n.panel = 500, n.cores = 1L)
# })

# # RUN THE FUNCTION
# outcome <- create_panels(data.mAb.panel, min.panel = 10, max.panel = 35, n.panel = 2500, n.cores = 4L, use.na.cor = "everything")
# 
# # INDICATE HERE WHERE YOU WANT TO SAVE THE OUTCOME OF THE FUNCTION
# write.table(outcome[[1]], "../results/allpanels_basic_220925.txt", row.names = FALSE)
# 
# write.table(outcome[[2]], "../results/panelopt_basic_220925.txt", row.names = FALSE)
# 
# panel_opt = outcome[[3]]
# saveRDS(panel_opt, file = "../results/allpanels_basic_220925.RDS")
