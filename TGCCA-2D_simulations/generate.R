rm(list = ls())

library(RGCCA)

### Get args from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Number of folds/Number of subjects per fold/SNR/Correlation factors/Shapes/Root path")
} else {
  N <- as.numeric(args[1])
  n <- as.numeric(args[2])
  SNR <- as.numeric(args[3])
  correlation_factors <- as.numeric(strsplit(args[4], ",")[[1]])
  shapes <- strsplit(args[5], ",")[[1]]
  root_path <- args[6]
  if (length(args) > 6) {
    noise_shapes <- strsplit(args[7], ",")[[1]]
  } else {
    noise_shapes <- NULL
  }
}

file.sources <- list.files(path = paste0(root_path, "/utils/."), pattern = "*.R", full.names = T)
sapply(file.sources, source, .GlobalEnv)

### Create output_dir
output_dir <- paste0(root_path, "/data/SNR_", SNR)
dir.create(file.path(output_dir), showWarnings = TRUE, recursive = TRUE)

### Save log file
fileConn <- file(paste0(output_dir, "/log.txt"))
writeLines(args, fileConn)
close(fileConn)

### Create N folds of data
set.seed(0)
for (i in 1:N) {
  if (i %% 10 == 1) print(paste0("Creating fold ", i))
  dataset <- generate_data(shapes, noise_shapes, shape_getter, correlation_factors, n, SNR, root_path)

  save(dataset, file = paste0(output_dir, "/fold_", i, ".Rdata"))
}
