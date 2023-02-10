rm(list = ls())

### Get args from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Blocks/Root path/Name of the experiment have to be specified")
} else {
  blocks <- as.numeric(strsplit(args[1], ",")[[1]])
  root_path <- args[2]
  name <- args[3]
  if (length(args) > 3) {
    folder <- args[4]
  } else {
    folder <- "results"
  }
}

file.sources <- list.files(path = paste0(root_path, "/utils/."), pattern = "*.R", full.names = T)
sapply(file.sources, source, .GlobalEnv)

### Set output_dir
output_dir <- paste0(root_path, "/", folder, "/", name)

### Load ground truth
shapes <- c(
  "2D_factors/square", "3D_factors/gas", "2D_factors/cross",
  "3D_factors/cross_little", "1D_factors/vector"
)
name_shapes <- c("square", "gas_3D", "cross", "cross_little_3D", "vector")
shapes <- shapes[blocks]
name_shapes <- name_shapes[blocks]
a_GT <- lapply(shapes, function(x) shape_getter(root_path, x))
V_GT <- lapply(a_GT, function(x) matrix(data.matrix(x), ncol = 1))
names(a_GT) <- names(V_GT) <- name_shapes

### Get SNR directories
directories <- list.dirs(paste0(root_path, "/", folder, "/", name), recursive = F, full.names = F)
SNR_dirs <- directories[grepl("SNR_", directories)]

SNR_levels <- sapply(SNR_dirs, function(x) as.numeric(gsub("SNR_", "", x)))
SNR_levels <- unname(SNR_levels)

results <- list()
times <- list()

### Compute angles and collect times for all folds
for (SNR in SNR_levels) {
  SNR_path <- paste0(root_path, "/", folder, "/", name, "/SNR_", SNR)
  SNR_results <- list()
  SNR_times <- c()

  files <- list.files(path = SNR_path, pattern = "fold_", full.names = T)
  n_folds <- length(files)

  for (i in 1:n_folds) {
    load(files[i])
    SNR_results[[i]] <- angles(V_GT, a)$angles
    SNR_times <- c(SNR_times, as.numeric(duration, units = "secs"))
  }
  results[[paste0("SNR_", SNR)]] <- SNR_results
  times[[paste0("SNR_", SNR)]] <- SNR_times
}

### Save results
save(results, times, file = paste0(output_dir, "/analysis.Rdata"))
