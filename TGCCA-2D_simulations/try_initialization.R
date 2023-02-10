rm(list = ls())

library(RGCCA)
library(abind)

### Get args from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Number of subjects per fold/Number of repetitions/SNR/Blocks/method/Root path/Name of the experiment")
} else {
  N_fold <- as.numeric(args[1])
  N_rep <- as.numeric(args[2])
  SNR <- as.numeric(args[3])
  blocks <- as.numeric(strsplit(args[4], ",")[[1]])
  method <- args[5]
  root_path <- args[6]
  name <- args[7]
  if (length(args) > 7) {
    # Get model params
    params <- args[8:length(args)]
  } else {
    params <- NULL
  }
}

file.sources <- list.files(path = paste0(root_path, "/utils/."), pattern = "*.R", full.names = T)
sapply(file.sources, source, .GlobalEnv)

### Create output_dir
output_dir <- paste0(root_path, "/initialization/", name, "/SNR_", SNR)
dir.create(file.path(output_dir), showWarnings = TRUE, recursive = TRUE)

### Save log file
fileConn <- file(paste0(output_dir, "/log.txt"))
writeLines(args, fileConn)
close(fileConn)

### Run experiment on each fold
set.seed(0)
next_fold <- 1

### Load data
N_loaded <- 0
complete_dataset <- list(X = NULL, Xu = NULL, a = NULL, V = NULL)
while (N_loaded < N_fold) {
  data_path <- paste0(root_path, "/data/SNR_", SNR, "/fold_", next_fold, ".Rdata")
  load(data_path)
  complete_dataset$X <- lapply(1:length(blocks), function(j) abind(complete_dataset$X[[j]], dataset$X[[blocks[j]]], along = 1))
  complete_dataset$Xu <- lapply(1:length(blocks), function(j) abind(complete_dataset$Xu[[j]], dataset$Xu[[blocks[j]]], along = 1))
  complete_dataset$a <- dataset$a[blocks]
  complete_dataset$V <- dataset$V[blocks]
  next_fold <- next_fold + 1
  N_loaded <- N_loaded + nrow(dataset$X[[1]])
}

for (i in 1:N_rep) {
  ### Run model
  start <- Sys.time()
  func <- get(paste0("run_", tolower(method)))
  a <- func(complete_dataset, params, root_path, seed = i, random_init = TRUE)
  end <- Sys.time()
  duration <- end - start

  ### Save results
  save(a, duration, file = paste0(output_dir, "/fold_", i, ".Rdata"))
}
