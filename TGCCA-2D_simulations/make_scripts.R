rm(list = ls())

### Get args from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Number of subjects per fold/Total number of subjects/Root path/Experiment")
} else {
  N <- as.numeric(args[1])
  N_total <- as.numeric(args[2])
  root_path <- args[3]
  experiment <- args[4]
}


SNR_levels <- c(0.1, 0.3, 0.5, 1)
if (experiment == "2B") {
  shapes <- c("gas", "cross_little")
  block_names <- c("Gas", "Cross (small)")
  models <- c(
    "twodcca", "tcca", "tgcca", "tcca", "sptgcca",
    "rtwodcca", "tcca", "tgcca", "tcca", "sptgcca",
    "rgcca", "svd"
  )
  names <- c(
    "2DCCA1", "TCCA1", "TGCCA1", "spTCCA1", "MGCCA",
    "2DCCA3", "TCCA3", "TGCCA3", "spTCCA3", "spTGCCA3",
    "RGCCA", "SVD"
  )
  params <- vector("list", length = length(models))
  params[[1]] <- list("1", "0.001", "T", "F") # 2DCCA1
  params[[2]] <- list("1,1", "0.001", "T", "F", "none", "10", "F") # TCCA1
  params[[3]] <- list("1", "0.001", "T", "horst", "F", "10", "10") # TGCCA1
  params[[4]] <- list("1,1", "0.001", "T", "F", "sep", "10", "F") # spTCCA1
  params[[5]] <- list("1", "0.001", "T", "horst", "F", "10", "10") # MGCCA
  params[[6]] <- list("3", "0.001", "T", "F") # 2DCCA3
  params[[7]] <- list("3,3", "0.001", "T", "F", "none", "10", "F") # TCCA3
  params[[8]] <- list("3", "0.001", "T", "horst", "F", "10", "10") # TGCCA3
  params[[9]] <- list("3,3", "0.001", "T", "F", "sep", "10", "F") # spTCCA3
  params[[10]] <- list("3", "0.001", "T", "horst", "F", "10", "10") # spTGCCA3
  params[[11]] <- list("0.001", "T", "horst", "F", "10", "10") # RGCCA
  params[[12]] <- list() # SVD
} else {
  shapes <- c("square", "gas", "cross", "cross_little", "vector")
  block_names <- c("Square", "Gas", "Cross", "Cross (small)", "Vector")
  models <- c(
    "tgcca", "sptgcca", "tgcca", "sptgcca", "rgcca", "svd"
  )
  names <- c(
    "TGCCA1", "MGCCA", "TGCCA3", "spTGCCA3", "RGCCA", "SVD"
  )
  params <- vector("list", length = length(models))
  params[[1]] <- list("1", "0.001", "T", "horst", "F", "10", "10") # TGCCA1
  params[[2]] <- list("1", "0.001", "T", "horst", "F", "10", "10") # MGCCA
  params[[3]] <- list("3", "0.001", "T", "horst", "F", "10", "10") # TGCCA3
  params[[4]] <- list("3", "0.001", "T", "horst", "F", "10", "10") # spTGCCA3
  params[[5]] <- list("0.001", "T", "horst", "F", "10", "10") # RGCCA
  params[[6]] <- list() # SVD
}

ref_shapes <- c("square", "gas", "cross", "cross_little", "vector")
blocks <- match(shapes, ref_shapes)

fileConn <- file(paste0(root_path, "/cmd/", experiment, "_", N, ".txt"), "w")
full_names <- list()

for (i in 1:length(models)) {
  full_names[[i]] <- paste0(names[i], "_B", paste(blocks, collapse = ""), "_", N)
  for (SNR in SNR_levels) {
    writeLines(paste0(
      "Rscript ",
      root_path,
      "/run.R ",
      N,
      " ",
      N_total,
      " ",
      SNR,
      " ",
      paste(blocks, collapse = ","),
      " ",
      models[i],
      " ",
      root_path,
      " ",
      full_names[[i]],
      " ",
      paste(params[[i]], collapse = " ")
    ), fileConn)
  }

  writeLines(paste0(
    "Rscript ",
    root_path,
    "/analyse.R ",
    paste(blocks, collapse = ","),
    " ",
    root_path,
    " ",
    full_names[[i]],
    "\n"
  ), fileConn)
}

writeLines(paste0(
  "Rscript ",
  root_path,
  "/compare.R ",
  paste(shapes, collapse = ","),
  " ",
  paste(SNR_levels, collapse = ","),
  " ",
  root_path,
  " ",
  paste(full_names, collapse = ","),
  " ",
  paste(names, collapse = ","),
  " ",
  "\"", paste(block_names, collapse = ","), "\""
), fileConn)

close(fileConn)
