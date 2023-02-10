rm(list = ls())

### Get args from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Number of subjects/SNR/Correlation factors")
} else {
  N <- as.numeric(args[1])
  SNR <- as.numeric(args[2])
  correlation_factors <- as.numeric(strsplit(args[3], ",")[[1]])
  root_path <- args[4]
}

shapes <- c(
  "2D_factors/square", "3D_factors/gas", "2D_factors/cross",
  "3D_factors/cross_little", "1D_factors/vector"
)
noise_shapes <- c(
  "2D_factors/information", "3D_factors/cup",
  "2D_factors/parking", "2D_factors/restaurant", "null"
)

source(paste0(root_path, "/utils/generate_data.R"), .GlobalEnv)

### Create output_dir
output_dir <- paste0(root_path, "/data/SNR_", SNR)
dir.create(file.path(output_dir), showWarnings = TRUE, recursive = TRUE)

### Save log file
fileConn <- file(paste0(output_dir, "/log.txt"))
writeLines(args, fileConn)
close(fileConn)

### Create 100 folds of data
set.seed(0)
for (i in 1:100) {
  if (i %% 10 == 1) print(paste0("Creating fold ", i))
  dataset <- generate_data(
    shapes, noise_shapes, shape_getter,
    correlation_factors, N, SNR, root_path
  )

  saveRDS(dataset, file = paste0(output_dir, "/fold_", i))
}
