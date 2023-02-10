rm(list = ls())

library(ggplot2)
library(viridis)

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
output_dir <- paste0(root_path, "/comparisons/", gsub(pattern = "  | |:", replacement = "_", x = Sys.time()))
dir.create(file.path(output_dir), showWarnings = TRUE, recursive = TRUE)

shapes <- c("Square", "Gas", "Cross", "Cross (small)", "Vector")
shapes <- shapes[blocks]

### Get SNR directories
directories <- list.dirs(paste0(root_path, "/", folder, "/", name), recursive = F, full.names = F)
SNR_dirs <- directories[grepl("SNR_", directories)]

SNR_levels <- sapply(SNR_dirs, function(x) as.numeric(gsub("SNR_", "", x)))
SNR_levels <- unname(SNR_levels)

SNR_mapping <- c(0.1, 0.3, 0.5, 1)
SNR_mapping_db <- c(-20, -10, -6, 0)
SNR_db_levels <- vapply(SNR_levels, function(x) SNR_mapping_db[match(x, SNR_mapping)], FUN.VALUE = 1.)

weights <- list()
complete_df <- c()

### Collect scores for all folds and models
for (SNR in SNR_levels) {
  SNR_path <- paste0(root_path, "/", folder, "/", name, "/SNR_", SNR)
  SNR_weights <- list()

  files <- list.files(path = SNR_path, pattern = "fold_", full.names = T)
  n_folds <- length(files)

  for (i in 1:n_folds) {
    load(files[i])
    SNR_weights[[i]] <- res$weights
  }
  weights[[paste0("SNR_", SNR)]] <- SNR_weights

  n_SNR <- length(SNR_levels)
  n_shapes <- length(shapes)
  n_rank <- 3
}


df <- data.frame(
  rank = rep(1:3, n_folds * n_shapes * n_SNR),
  SNR = rep(SNR_db_levels, each = n_folds * n_shapes * n_rank),
  shape = factor(rep(shapes, each = n_rank, n_folds * n_SNR), levels = shapes),
  y = unlist(weights)^2
)

df <- aggregate(. ~ rank + shape + SNR, df, mean)

width <- 60
height <- 42

fig_final <- ggplot(
  data = df, aes(x = factor(SNR), y = y, fill = factor(rank))
) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("SNR (in dB)") +
  ylab("Weights") +
  guides(fill = guide_legend(title = "Factor")) +
  scale_fill_manual(values = cividis(3), guide = FALSE) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.text.x = element_text(size = 40, angle = 0),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 50),
    axis.title.y = element_text(size = 50),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 50),
    legend.position = "right",
    strip.text.x = element_text(size = 40),
    strip.background = element_blank(),
    legend.background = element_rect(colour = NA, fill = NA)
  ) +
  facet_wrap(~shape, strip.position = "top")
fig_final
ggsave(
  filename = paste0(output_dir, "/comparison.pdf"),
  plot = fig_final, width = width, height = height, units = "cm"
)
