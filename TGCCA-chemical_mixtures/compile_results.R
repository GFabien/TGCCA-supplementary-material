library(RGCCA)
library(ggplot2)
library(viridis)
library(R.matlab)

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

N <- 28

### Preprocess according to the paper procedure
X <- X / sqrt(sum(X^2))
X_centered <- array(scale(matrix(as.vector(X), nrow = nrow(X)), center = T, scale = F), dim = dim(X))

Z <- scale(Z, center = F, scale = apply(Z, 2, sd, na.rm = TRUE))
Z <- Z / sqrt(sum(Z^2))
Z_centered <- scale(Z, center = T, scale = F)

### Utility functions
khatri_rao <- function(x, y) {
  if ((length(dim(x)) != 2) || (length(dim(y)) != 2)) {
    stop_rgcca("x and y must be matrices")
  }
  ncol <- dim(x)[2]
  if (ncol != dim(y)[2]) {
    stop_rgcca("x and y must have the same number of columns")
  }
  res <- sapply(1:ncol, function(i) x[, i] %x% y[, i])
  return(res)
}

list_khatri_rao <- function(factors) {
  Reduce("khatri_rao", rev(factors))
}

kron_sum <- function(factors) {
  apply(list_khatri_rao(factors), 1, sum)
}

weighted_kron_sum <- function(factors, weights) {
  list_khatri_rao(factors) %*% weights
}

### Align factors with ground truth
align <- function(Y, Y_true) {
  idx <- c()
  Y_final <- c()
  for (i in 1:ncol(Y_true)) {
    idx <- which.max(abs(cor(Y_true[, i], Y)))
    if (cor(Y_true[, i], Y[, idx]) < 0) Y[, idx] <- -Y[, idx]
    Y_final <- cbind(Y_final, Y[, idx])
    Y <- Y[, -idx, drop = F]
  }
  return(Y_final)
}

### Load data and compute factors
dir_path <- "./models/"
Y_true <- apply(normalized_concentrations, 2, function(x) x / norm(x, type = "2"))

load_Y <- function(i, dir_path, name = "TGCCA") {
  if (grepl("TGCCA", name, fixed = T)) {
    model <- readRDS(paste0(dir_path, "/", name, "/", i, ".rds"))
    tmp <- list_khatri_rao(lapply(model$factors[[1]], function(x) x))
    Y <- matrix(as.vector(X), nrow = nrow(X)) %*% tmp
  } else {
    model <- readMat(paste0(dir_path, "/", name, "/", i, ".mat"))
    Y <- model$Zhat[[1]][[1]]$u[[1]][[1]]
  }
  Y <- align(Y, Y_true)
  return(Y)
}

load_time <- function(i, dir_path, name = "TGCCA") {
  if (grepl("TGCCA", name, fixed = T)) {
    time <- as.numeric(readRDS(paste0(dir_path, "/", name, "/times/", i, ".rds")))
  } else {
    time <- drop(readMat(paste0(dir_path, "/", name, "/times/", i, ".mat"))[[1]])
  }
  return(time)
}

compute_model_crit <- function(i, dir_path, name = "TGCCA") {
  if (grepl("TGCCA", name, fixed = T)) {
    model <- readRDS(paste0(dir_path, "/", name, "/", i, ".rds"))
    crit <- -Reduce("+", lapply(model$crit, function(x) x[length(x)]))
  } else {
    model <- readMat(paste0(dir_path, "/", name, "/", i, ".mat"))
    A <- model$Zhat[[1]][[1]]$u[[1]][[1]]
    B <- model$Zhat[[1]][[1]]$u[[2]][[1]]
    C <- model$Zhat[[1]][[1]]$u[[3]][[1]]
    V <- model$Zhat[[2]][[1]]$u[[2]][[1]]
    w1 <- drop(model$Zhat[[1]][[1]]$lambda)
    w2 <- drop(model$Zhat[[2]][[1]]$lambda)
    crit <- sum((as.vector(X) - weighted_kron_sum(list(A, B, C), w1))^2) + sum((as.vector(Z) - weighted_kron_sum(list(A, V), w2))^2)
  }
  return(crit)
}

get_best_model <- function(dir_path, name = "TGCCA", n_max = 100) {
  best_crit <- Inf
  best_i <- 0
  for (i in 1:n_max) {
    crit <- compute_model_crit(i, dir_path, name)
    if (crit < best_crit) {
      best_crit <- crit
      best_i <- i
    }
  }
  return(best_i)
}

performances <- function(dir_path, name = "TGCCA", n_max = 100) {
  cosines <- c()
  times <- c()
  for (i in 1:n_max) {
    Y <- load_Y(i, dir_path, name)
    time <- load_time(i, dir_path, name)
    cosines <- rbind(cosines, abs(diag(crossprod(apply(Y, 2, function(x) x / norm(x, type = "2")), Y_true))))
    times <- c(times, time)
  }
  colnames(cosines) <- colnames(Y_true)
  return(list(
    cosines    = cosines,
    times      = times,
    avg_cosine = apply(cosines, 2, mean),
    sd_cosine  = apply(cosines, 2, sd),
    avg_time   = mean(times),
    sd_time    = sd(times)
  ))
}

n_max <- 100

res_TGCCA <- performances(dir_path, "TGCCA", n_max)
res_TGCCA_1 <- performances(dir_path, "TGCCA1", n_max)
res_CMTF <- performances(dir_path, "CMTF", n_max)
res_ACMTF <- performances(dir_path, "ACMTF", n_max)

saveRDS(res_TGCCA, paste0(dir_path, "/TGCCA/results.rds"))
saveRDS(res_TGCCA_1, paste0(dir_path, "/TGCCA1/results.rds"))
saveRDS(res_CMTF, paste0(dir_path, "/CMTF/results.rds"))
saveRDS(res_ACMTF, paste0(dir_path, "/ACMTF/results.rds"))

# res_TGCCA = readRDS(paste0(dir_path, "/TGCCA/results.rds"))
# res_TGCCA_1 = readRDS(paste0(dir_path, "/TGCCA1/results.rds"))
# res_CMTF = readRDS(paste0(dir_path, "/CMTF/results.rds"))
# res_ACMTF = readRDS(paste0(dir_path, "/ACMTF/results.rds"))

names <- c("TGCCA", "CMTF", "ACMTF")

### Create latex code to make tables from results
write_latex_table <- function() {
  cols <- paste(rep("l", ncol(Y_true)), collapse = "")
  name_cols <- paste(colnames(Y_true), collapse = " & ")
  header <- paste0(
    "\\begin{tabular}{l", cols, "l} \n",
    "\\toprule \n",
    "Model & ", name_cols, " & Computation time \\\\ \n",
    "\\midrule \n"
  )
  bottom <- "\\bottomrule \n\\end{tabular}"
  lines <- ""
  nb_lines <- length(names)
  for (i in 1:nb_lines) {
    lines <- paste0(lines, names[i], " & ")
    res <- get(paste0("res_", names[i]))
    for (j in 1:ncol(Y_true)) {
      lines <- paste0(
        lines,
        format(round(res$avg_cosine[j], 3), nsmall = 3), " $\\pm$ ",
        format(round(res$sd_cosine[j], 3), nsmall = 3), " & "
      )
    }
    lines <- paste0(
      lines,
      format(round(res$avg_time, 3), nsmall = 3), " $\\pm$ ",
      format(round(res$sd_time, 3), nsmall = 3),
      " \\\\ \n"
    )
  }
  return(paste0(header, lines, bottom))
}

fileConn <- file(paste0(dir_path, "/table.txt"))
writeLines(write_latex_table(), fileConn)
close(fileConn)


### Boxplot
source(file = "./shift_legend.R")

names <- c("TGCCA", "CMTF", "ACMTF")
df <- data.frame(cosine = c(
  as.vector(res_TGCCA$cosines),
  as.vector(res_CMTF$cosines),
  as.vector(res_ACMTF$cosines)
))
df$mol <- factor(rep(colnames(Y_true), length(names), each = n_max), levels = colnames(Y_true))
df$names <- rep(names, each = n_max * ncol(Y_true))

p <- ggplot(df, aes(
  x = mol, y = cosine,
  group = interaction(names, factor(mol)), fill = names
)) +
  geom_boxplot() +
  scale_fill_manual(values = cividis(length(names))) +
  xlab("Chemicals") +
  ylab("Cosine between estimated and real concentrations") +
  theme_light() +
  labs(fill = "Model") +
  theme(text = element_text(size = 20))
ggsave(filename = paste0(dir_path, "/boxplot.pdf"), plot = p)

### Represent the reconstructed chemicals for the best models
i_TGCCA <- get_best_model(dir_path, "TGCCA", n_max)
Y <- align(load_Y(i_TGCCA, dir_path, "TGCCA"), Y_true)
Y_TGCCA <- data.frame(apply(Y, 2, function(x) x / norm(x, type = "2")))

i_TGCCA_1 <- get_best_model(dir_path, "TGCCA1", n_max)
Y <- align(load_Y(i_TGCCA_1, dir_path, "TGCCA1"), Y_true)
Y_TGCCA_1 <- data.frame(apply(Y, 2, function(x) x / norm(x, type = "2")))

i_CMTF <- get_best_model(dir_path, "CMTF", n_max)
Y <- align(load_Y(i_CMTF, dir_path, "CMTF"), Y_true)
Y_CMTF <- data.frame(apply(Y, 2, function(x) x / norm(x, type = "2")))

i_ACMTF <- get_best_model(dir_path, "ACMTF", n_max)
Y <- align(load_Y(i_ACMTF, dir_path, "ACMTF"), Y_true)
Y_ACMTF <- data.frame(apply(Y, 2, function(x) x / norm(x, type = "2")))

GT <- data.frame(Y_true)

## Y_TGCCA vs Y_TGCCA_1
names <- c("GT", "TGCCA2", "TGCCA1")
df <- data.frame(Y = c(
  as.vector(Y_true),
  as.vector(data.matrix(Y_TGCCA)),
  as.vector(data.matrix(Y_TGCCA_1))
))
df$names <- factor(rep(names, each = N * ncol(Y_true)), levels = names)
df$mol <- factor(rep(colnames(Y_true), each = N, length(names)), levels = colnames(Y_true))
df$i <- rep(1:N, length(names) * ncol(Y_true))
p <- ggplot(data = df, aes(x = i, y = Y, color = names, group = names, linetype = names)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = cividis(length(names))) +
  xlab("Mixtures") +
  ylab("Normalized concentration") +
  theme_light() +
  facet_wrap(~mol, ncol = 2, scales = "free_y", strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(colour = "black"), text = element_text(size = 20)) +
  labs(color = "Model", linetype = "Model") +
  scale_x_continuous(breaks = seq(0, 30, 5))
p <- shift_legend(p)
ggsave(filename = paste0(dir_path, "/TGCCA2_vs_TGCCA1.pdf"), plot = p)

## Y_TGCCA vs CMTF vs ACMTF
names <- c("GT", "TGCCA", "CMTF", "ACMTF")
df <- data.frame(Y = c(
  as.vector(Y_true),
  as.vector(data.matrix(Y_TGCCA)),
  as.vector(data.matrix(Y_CMTF)),
  as.vector(data.matrix(Y_ACMTF))
))
df$names <- factor(rep(names, each = N * ncol(Y_true)), levels = names)
df$mol <- factor(rep(colnames(Y_true), each = N, length(names)), levels = colnames(Y_true))
df$i <- rep(1:N, length(names) * ncol(Y_true))
p <- ggplot(data = df, aes(x = i, y = Y, color = names, group = names, linetype = names)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = cividis(length(names))) +
  xlab("Mixtures") +
  ylab("Normalized concentration") +
  theme_light() +
  facet_wrap(~mol, ncol = 2, scales = "free_y", strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(colour = "black"), text = element_text(size = 20)) +
  labs(color = "Model", linetype = "Model") +
  scale_x_continuous(breaks = seq(0, 30, 5))
p <- shift_legend(p)
ggsave(filename = paste0(dir_path, "/TGCCA_vs_CMTF_vs_ACMTF.pdf"), plot = p)
