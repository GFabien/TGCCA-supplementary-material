#' Plot the connection between blocks
#'
#' @inheritParams plot_ind
#' @inheritParams plot2D
#' @return A dataframe with tuples of connected blocks
#' @examples
#' library(igraph)
#' data("Russett")
#' blocks <- list(
#'   agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#' rgcca_out <- rgcca(blocks)
#' plot_network(rgcca_out)
#' @export
plot_network <- function(rgcca_res,
                         title = paste0(
                           "Common rows between blocks : ",
                           NROW(rgcca_res$blocks[[1]])
                         ),
                         cex_main = 2,
                         cex_point = 23,
                         colors = c("#eee685", "gray60")) {
  stopifnot(is(rgcca_res, "rgcca"))
  title <- paste0(title, collapse = " ")
  check_colors(colors)

  load_libraries("igraph")

  # Avoid random
  set.seed(1)
  `V<-` <- igraph::`V<-`
  `E<-` <- igraph::`E<-`
  V <- E <- NULL

  nodes <- get_nodes(rgcca_res)
  edges <- get_edges(rgcca_res)

  par <- ifelse("sparsity" %in% names(nodes), "sparsity", "tau")

  net <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )

  V(net)$color <- as.vector(colors[1])
  V(net)$label <- paste(
    nodes$id,
    "\nP =",
    nodes$P,
    paste0("\n", par, " ="),
    nodes[, par],
    "\nN =",
    nodes$nrow,
    sep = " "
  )
  V(net)$shape <- "square"
  E(net)$width <- E(net)$weight * 2
  plot(
    net,
    edge.color = colors[2],
    edge.lty = 2,
    vertex.frame.color = colors[2],
    vertex.label.color = "black",
    vertex.label.dist = 6,
    vertex.label.degree = 1.5,
    vertex.size = cex_point
  )
  title(title, cex.main = cex_main)
}
