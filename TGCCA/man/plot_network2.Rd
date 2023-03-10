% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_network2.R
\name{plot_network2}
\alias{plot_network2}
\title{Plot the connection between blocks (dynamic plot)}
\usage{
plot_network2(
  rgcca_res,
  title = paste0("Common rows between blocks : ", NROW(rgcca_res$blocks[[1]])),
  cex_main = 20,
  colors = c("#eee685", "gray")
)
}
\arguments{
\item{rgcca_res}{A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})}

\item{title}{A character giving the title of the plot}

\item{cex_main}{An integer for the size of the title}

\item{colors}{Either a vector of integers (each integer corresponding to a color) or of characters corresponding to names of colors (as "blue",see colors()) or RGB code ("#FFFFFF").}
}
\value{
A dataframe with tuples of connected blocks
}
\description{
Plot the connection between blocks (dynamic plot)
}
\examples{
library(visNetwork)
data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
rgcca_out <- rgcca(blocks)
plot_network2(rgcca_out)
}
