### Get working directory from command line
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

dataset <- readRDS("./data/data.rds")
# 2 blocks:
#   dataset$Y: (28 x 13324 x 8) NMR block   (mixtures x chemical shift x gradient levels)
#   dataset$Z: (28 x 168)       LC-MS block (mixtures x features)

X <- dataset$Y$data
Z <- dataset$Z$data

# Load concentrations and remove mixture 27
concentrations <- readRDS("./data/concentrations.rds")
concentrations <- concentrations[-27, ]
normalized_concentrations <- apply(concentrations, 2, function(x) x / norm(x, type = "2"))

### Preprocess according to the paper procedure and center
X <- X / sqrt(sum(X^2))
X <- array(scale(matrix(as.vector(X), nrow = nrow(X)), center = T, scale = F), dim = dim(X))

Z <- scale(Z, center = F, scale = apply(Z, 2, sd, na.rm = TRUE))
Z <- Z / sqrt(sum(Z^2))
Z <- scale(Z, center = T, scale = F)

### Apply TGCCA
library(RGCCA)
n_repeat <- 100

### Run 100 times TGCCA rank 2 on the first component and rank 1 after
ranks <- matrix(c(2, rep(1, 7)), 4, 2)
for (i in 1:n_repeat) {
  start_time <- Sys.time()
  fit_mgcca <- rgcca(
    blocks = list(EEM = X, LCMS = Z), type = "mgcca",
    tau = 1, ranks = ranks, verbose = F, prescaling = T, ncomp = 4,
    init = "random"
  )
  end_time <- Sys.time()
  duration <- end_time - start_time
  saveRDS(fit_mgcca, file = paste0("./models/TGCCA/", i, ".rds"))
  saveRDS(duration, file = paste0("./models/TGCCA/times/", i, ".rds"))
}

### Run 100 times TGCCA rank 1
for (i in 1:n_repeat) {
  start_time <- Sys.time()
  fit_mgcca <- rgcca(
    blocks = list(EEM = X, LCMS = Z), type = "mgcca",
    tau = 1, ranks = 1, verbose = F, prescaling = T, ncomp = 5,
    init = "random"
  )
  end_time <- Sys.time()
  duration <- end_time - start_time
  saveRDS(fit_mgcca, file = paste0("./models/TGCCA1/", i, ".rds"))
  saveRDS(duration, file = paste0("./models/TGCCA1/times/", i, ".rds"))
}
