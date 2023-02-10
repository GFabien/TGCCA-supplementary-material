# rand_mat <- function(x) matrix(runif(9), 3, 3)
# A = lapply(1:3, rand_mat)
# check_blocks(A)
# names(A) <- LETTERS[1:3]
# check_blocks(A[1])
# check_blocks(A)
# row.names(A[[1]]) <- letters[1:3]
# check_blocks(A)
# for(i in 1:3)
#   row.names(A[[i]]) <- letters[(0+i*3):(2+i*3)]
# check_blocks(A)
# for(i in 1:3)
#   row.names(A[[i]]) <- letters[1:3]
# A[[1]][2, 3] <- NA
# for(i in 1:3)
#   colnames(A[[i]]) <- letters[(0+i*3):(2+i*3)]
# check_blocks(A,add_NAlines=TRUE)
# A[[1]][2, 3] <- "character"
# check_blocks(A)
# A[[1]][2, 3] <- runif(1)
# init : boolean (FALSE by default) for the first block checking

check_blocks <- function(blocks, init = FALSE, n = 2,
                         add_NAlines = FALSE, allow_unnames = TRUE,
                         quiet = FALSE, no_character = FALSE) {
  msg <- ""
  if (is.matrix(blocks)) blocks <- list(blocks)
  if (!is.list(blocks)) stop_rgcca(paste(msg, "is not a list."))
  if (!init && length(blocks) < n) {
    stop_rgcca(paste(msg, "should at least have two elements."))
  }

  # Completing block names
  if (is.null(names(blocks))) {
    names(blocks) <- paste0("block", 1:length(blocks))
    message("Blocks are unnamed and automatically labeled as block1, ..., blockJ")
  }
  # Gestion of the case of one variable only
  vector_idx <- sapply(1:length(blocks), function(x) length(dim(blocks[[x]])) < 2)
  blocks[vector_idx] <- lapply(blocks[vector_idx], as.matrix)
  nameBlocks <- names(blocks)

  # Dealing with rownames (if they are all missing)
  if (all(sapply(blocks, function(x) is.null(row.names(x))))) {
    if (sd(sapply(blocks, function(x) NROW(x))) == 0 && allow_unnames) {
      blocks <- lapply(
        blocks,
        function(x) {
          rownames(x) <- paste0("S", 1:(NROW(x)))
          return(x)
        }
      )
      message("Rownames are missing and automatically labeled as S1, ..., Sn \n")
    } else {
      stop_rgcca(paste(msg, "Blocks should have rownames.\n "))
    }
  }


  # Dealing with colnames for univariate block
  blocks1 <- lapply(
    1:length(blocks),
    function(i) {
      if (NCOL(blocks[[i]]) == 1 & is.null(colnames(blocks[[i]]))) {
        colnames(blocks[[i]]) <- nameBlocks[i]
        return(blocks[[i]])
      } else {
        return(blocks[[i]])
      }
    }
  )

  blocks <- blocks1
  names(blocks) <- nameBlocks

  # Dealing with other dimnames
  if (any(unlist(sapply(blocks, function(x) sapply(dimnames(x), is.null))))) {
    if (allow_unnames) {
      message("Some dimnames are missing and automatically labeled as
              block1_1_1, ..., block1_2_Im1, ..., blockJ_mJ_ImJ \n")
      for (i in 1:length(blocks)) {
        for (mode in seq(1, length(dim(blocks[[i]])))[-1]) {
          if (is.null(dimnames(blocks[[i]])[[mode]])) {
            dimnames(blocks[[i]])[[mode]] <- paste(
              names(blocks)[i], mode - 1, 1:(dim(blocks[[i]])[mode]),
              sep = "_"
            )
          }
        }
      }
    } else {
      stop_rgcca(paste(msg, "Blocks should have dimnames.\n "))
    }
  }

  # Dealing with duplicated dimnames (except rownames)
  if (any(duplicated(unlist(sapply(blocks, function(x) dimnames(x)[-1]))))) {
    if (!quiet) {
      message("Dimnames are duplicated and modified to avoid confusion \n")
    }

    for (i in 1:length(blocks)) {
      for (mode in seq(1, length(dim(blocks[[i]])))[-1]) {
        if (length(dim(blocks[[i]])) > 2) {
          dimnames(blocks[[i]])[[mode]] <- paste(
            names(blocks)[i], mode - 1,
            dimnames(blocks[[i]])[[mode]],
            sep = "_"
          )
        } else {
          dimnames(blocks[[i]])[[mode]] <- paste(
            names(blocks)[i], dimnames(blocks[[i]])[[mode]],
            sep = "_"
          )
        }
      }
    }
  }

  # If one rownames is missing but the size of blocks is correct
  if (any(sapply(blocks, function(x) is.null(row.names(x))))) {
    matrixOfRownames <- Reduce(cbind, lapply(blocks, row.names))
    if (sum(!apply(
      matrixOfRownames, 2,
      function(x) x == matrixOfRownames[, 1]
    )) == 0) {
      blocks <- lapply(
        blocks,
        function(x) {
          row.names(x) <- matrixOfRownames[, 1]
          return(x)
        }
      )
    }
  }

  lapply(
    blocks,
    function(x) {
      resdup <- duplicated(rownames(x))
      if (sum(resdup) != 0) {
        if (!quiet) {
          warning(paste0(
            "Duplicated rownames were removed: ",
            rownames(x)[resdup], "\n"
          ))
        }
      }
    }
  )
  inters_rows <- Reduce(intersect, lapply(blocks, row.names))

  if (length(inters_rows) == 0) {
    stop_rgcca(paste(msg, "elements of the list should have at least
                         one common rowname.\n "))
  }

  equal_rows <- Reduce(identical, lapply(blocks, row.names))

  # If add_NAlines=FALSE, taking the intersection_list
  if (!add_NAlines) {
    if (length(blocks) > 1 && !equal_rows) blocks <- common_rows(blocks)
  }

  if (init) {
    blocks <- remove_null_sd(blocks)
    for (i in seq(length(blocks))) {
      attributes(blocks[[i]])$nrow <- nrow(blocks[[i]])
    }
  }

  if (no_character) {
    if (any(sapply(blocks, is.character2))) {
      stop(paste(msg, "Blocks contain non-numeric values."))
    }

    for (i in seq(length(blocks))) {
      if (is.character(blocks[[i]])) {
        blocks[[i]] <- to_numeric(blocks[[i]])
      }
    }
  }

  # Add lines if subjects are missing
  if (add_NAlines) {
    union_rows <- Reduce(union, lapply(blocks, row.names))
    blocks2 <- lapply(nameBlocks, function(name) {
      # if some subjects are missing (in the rownames)
      # TODO: Adapt this to handle missing subjects when there are tensor blocks
      if (sum(!union_rows %in% rownames(blocks[[name]])) != 0) {
        message("Some subjects are not present in all blocks. NA rows were added to get blocks with appropriate dimensions")
        y <- matrix(NA,
          length(union_rows),
          ncol = ifelse(is.null(dim(blocks[[name]])),
            1,
            NCOL(blocks[[name]])
          )
        )
        if (is.null(dim(blocks[[name]]))) {
          colnames(y) <- name
        } else {
          colnames(y) <- colnames(blocks[[name]])
        }
        rownames(y) <- union_rows
        y[rownames(blocks[[name]]), ] <- unlist(blocks[[name]])
        return(y)
      } else {
        if (NCOL(blocks[[name]]) == 1) {
          y <- matrix(blocks[[name]][union_rows, ], ncol = 1)
          rownames(y) <- union_rows
          colnames(y) <- name
        } else {
          # Select the lines of the arrays with 2 or more dimensions
          y <- apply(blocks[[name]], -1, "[", union_rows)
        }
        return(y)
      }
    })
    names(blocks2) <- nameBlocks
    blocks <- blocks2
  }
  invisible(blocks)
}
