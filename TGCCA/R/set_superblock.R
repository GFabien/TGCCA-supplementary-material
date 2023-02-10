# #' @export
set_superblock <- function(blocks, superblock = FALSE, type = "rgcca", verbose = TRUE) {
  if (superblock | tolower(type) == "pca") {
    blocks[["superblock"]] <- Reduce(cbind, blocks)
  }
  return(blocks)
}
