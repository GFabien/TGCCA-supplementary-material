rm(list = ls())

library(ggplot2)

### Get args from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Blocks/SNR levels/Root path/Names of the experiments/Short names/Block names have to be specified")
} else {
  blocks <- strsplit(args[1], ",")[[1]]
  SNR_levels <- as.numeric(strsplit(args[2], ",")[[1]])
  root_path <- args[3]
  names <- strsplit(args[4], ",")[[1]]
  short_names <- strsplit(args[5], ",")[[1]]
  block_names <- strsplit(args[6], ",")[[1]]
  if (length(args) > 6) {
    output_dir <- paste0(args[7])
  } else {
    output_dir <- paste0(root_path, "/comparisons/", gsub(pattern = "  | |:", replacement = "_", x = Sys.time()))
  }
  if (length(args) > 7) {
    folder <- args[8]
  } else {
    folder <- "results"
  }
}

file.sources <- list.files(path = paste0(root_path, "/utils/."), pattern = "*.R", full.names = T)
sapply(file.sources, source, .GlobalEnv)

### Create output_dir
dir.create(file.path(output_dir), showWarnings = TRUE, recursive = TRUE)

### Save log file
fileConn <- file(paste0(output_dir, "/log.txt"))
writeLines(args, fileConn)
close(fileConn)

### Gather all analysis in a single data.frame
complete_df <- c()

for (i in 1:length(names)) {
  name <- names[i]
  short_name <- short_names[i]

  load(paste0(root_path, "/", folder, "/", name, "/analysis.Rdata"))
  results <- results[names(results) %in% paste0("SNR_", SNR_levels)]
  times <- times[names(times) %in% paste0("SNR_", SNR_levels)]

  results <- lapply(results, function(x) lapply(x, function(y) y[blocks]))
  # shapes   = names(results[[1]][[1]])
  shapes <- block_names

  n_SNR <- length(SNR_levels)
  n_simu <- length(results[[1]])
  n_shapes <- length(shapes)

  df <- data.frame(
    name = rep(short_name, n_simu * n_SNR * n_shapes),
    SNR = rep(SNR_levels, each = n_simu * n_shapes),
    shape = rep(shapes, n_simu * n_SNR),
    y = unlist(results),
    t = rep(unlist(times), each = n_shapes)
  )
  complete_df <- rbind(complete_df, df)
}

complete_df$shape <- factor(complete_df$shape, levels = block_names)
complete_df$name <- factor(complete_df$name, levels = short_names)

### Display and save results
width <- 60
height <- 42

fig_final <- ggplot(
  complete_df, aes(
    x = 20 * log10(SNR), y = data.matrix(y),
    group = interaction(factor(SNR), name), fill = name
  )
) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_grey(start = 0.5, end = 1) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("") +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.text.x = element_text(size = 40, angle = 0),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 60),
    axis.title.y = element_text(size = 60),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 0),
    legend.position = "right",
    strip.text.x = element_text(size = 40),
    legend.background = element_rect(colour = NA, fill = NA)
  ) +
  facet_wrap(~shape)
ggsave(
  filename = paste0(output_dir, "/comparison.pdf"),
  plot = fig_final, width = width, height = height, units = "cm"
)

### Create latex code to make tables from results
write_latex_table <- function(df, blocks, SNR) {
  median_info <- aggregate(. ~ name + shape + SNR, df[df$SNR == SNR, ], median)
  q1_info <- aggregate(. ~ name + shape + SNR, df[df$SNR == SNR, ], quantile, 0.025)
  q3_info <- aggregate(. ~ name + shape + SNR, df[df$SNR == SNR, ], quantile, 0.975)
  cols <- paste(rep("l", length(blocks)), collapse = "")
  names <- paste(blocks, collapse = " & ")
  header <- paste0(
    "\\begin{tabular}{l", cols, "l} \n",
    "\\toprule \n",
    "Model & ", names, " & Computation time \\\\ \n",
    "\\midrule \n"
  )
  bottom <- "\\bottomrule \n\\end{tabular}"
  lines <- ""
  nb_lines <- nrow(median_info)
  for (i in 1:(nb_lines / length(blocks))) {
    lines <- paste0(lines, median_info[["name"]][i], " & ")
    for (j in 0:(length(blocks) - 1)) {
      lines <- paste0(
        lines,
        format(round(median_info[["y"]][i + j * nb_lines / length(blocks)], 2), nsmall = 2), " (",
        format(round(q1_info[["y"]][i + j * nb_lines / length(blocks)], 2), nsmall = 2), ", ",
        format(round(q3_info[["y"]][i + j * nb_lines / length(blocks)], 2), nsmall = 2), ") & "
      )
    }
    lines <- paste0(
      lines,
      format(round(median_info[["t"]][i], 2), nsmall = 2), " (",
      format(round(q1_info[["t"]][i], 2), nsmall = 2), ", ",
      format(round(q3_info[["t"]][i], 2), nsmall = 2),
      ") \\\\ \n"
    )
  }
  return(paste0(header, lines, bottom))
}

for (SNR in SNR_levels) {
  fileConn <- file(paste0(output_dir, "/table_SNR_", SNR, ".txt"))
  writeLines(write_latex_table(complete_df, block_names, SNR), fileConn)
  close(fileConn)
}
