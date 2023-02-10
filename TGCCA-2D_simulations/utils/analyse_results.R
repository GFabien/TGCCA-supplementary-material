angle <- function(v, w) {
  abs(drop(crossprod(v, w)) / (norm(as.vector(v), type = "2") * norm(as.vector(w), type = "2")))
}

angles <- function(true, pred) {
  J <- length(true)
  res <- sapply(1:J, function(j) angle(true[[j]], pred[[j]]))
  names(res) <- names(true)
  return(list(angles = res, mean = mean(res), sd = sd(res)))
}

plot_result <- function(true, pred, plot_true = F) {
  J <- length(true)
  for (j in 1:J) {
    if (plot_true) {
      pheatmap::pheatmap(true[[j]],
        cluster_rows = F, cluster_cols = F,
        main = paste0(names(true)[[j]], " real canonical factor")
      )
    }
    pheatmap::pheatmap(matrix(pred[[j]], nrow = nrow(true[[j]])),
      cluster_rows = F, cluster_cols = F,
      main = paste0(names(true)[[j]], " predicted canonical factor")
    )
  }
}
