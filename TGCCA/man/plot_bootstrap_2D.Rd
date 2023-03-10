% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bootstrap_2D.R
\name{plot_bootstrap_2D}
\alias{plot_bootstrap_2D}
\title{Plot a bootstrap in 2D}
\usage{
plot_bootstrap_2D(
  b = NULL,
  df_b = NULL,
  x = "bootstrap_ratio",
  y = "occurrences",
  title = paste("Variable selection \\nby", attributes(b)$n_boot, "bootstraps"),
  colors = NULL,
  cex = 1,
  cex_main = 14 * cex,
  cex_sub = 12 * cex,
  cex_point = 3 * cex,
  cex_lab = 10 * cex,
  comp = 1,
  i_block = NULL,
  collapse = FALSE,
  n_cores = parallel::detectCores() - 1
)
}
\arguments{
\item{b}{A boostrap object \code{\link[RGCCA]{bootstrap}}}

\item{df_b}{A get_bootstrap object \code{\link[RGCCA]{get_bootstrap}}}

\item{x}{A character for the index to plot in x-axis (see details).}

\item{y}{A character for the index to plot in y-axis (see details).}

\item{title}{A character giving the title of the plot}

\item{colors}{Either a vector of integers (each integer corresponding to a color) or of characters corresponding to names of colors (as "blue",see colors()) or RGB code ("#FFFFFF").}

\item{cex}{An integer for the size of the plot parameters}

\item{cex_main}{An integer for the size of the title}

\item{cex_sub}{An integer for the size of the subtitle}

\item{cex_point}{An integer for the size of the points or the text in the plot}

\item{cex_lab}{An integer for the size of the axis titles}

\item{comp}{An integer giving the index of the analysis components}

\item{i_block}{An integer giving the index of a list of blocks}

\item{collapse}{A boolean to combine the variables of each block as result}

\item{n_cores}{Number of cores for parallelization}
}
\description{
Graph of the best variables from a bootstrap with, in x-axis, the number of
non-zero occurences (SGCCA) or the significant 95% bootstrap
intervals (RGCCA; '*' or 'ns'; see 'p.vals' in details for
\code{\link[RGCCA]{get_bootstrap}}). In in y-axis are the bootstrap-ratios (mean/sd) .
Negative weights are colored in red and the positive ones are in green.
}
\details{
\itemize{
\item 'estimate' for RGCCA weights
\item 'bootstrap_ratio' for the mean of the bootstrap weights / their standard error
\item 'sign' for significant 95% bootstrap interval
\item 'occurrences' for non-zero occurences
\item 'mean' for the mean of the bootstrap weights
}
}
\examples{
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
rgcca_out <- rgcca(blocks, sparsity = 0.75, type = "sgcca")
boot <- bootstrap(rgcca_out, 2, n_cores = 1)
plot_bootstrap_2D(boot, n_cores = 1)
rgcca_out <- rgcca(blocks)
boot <- bootstrap(rgcca_out, 2, n_cores = 1)
selected.var <- get_bootstrap(boot, n_cores = 1, display_order = TRUE)
plot_bootstrap_2D(boot, n_cores = 1)
plot_bootstrap_2D(df_b = selected.var, n_cores = 1)
}
\seealso{
\code{\link[RGCCA]{bootstrap}}, \code{\link[RGCCA]{get_bootstrap}}
}
