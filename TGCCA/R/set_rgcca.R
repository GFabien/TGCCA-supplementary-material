# inds : individuals removed from the blocks (for crossvalidation)
# affects same parameters as rgcca_res (or other ones if specified) on a subset of individuals determined by inds (individuals to remove)
set_rgcca <- function(rgcca_res,
                      blocks = NULL,
                      connection = NULL,
                      tau = 1,
                      sparsity = 1,
                      ncomp = NULL,
                      scheme = NULL,
                      init = NULL,
                      bias = TRUE,
                      tol = 1e-03,
                      type = NULL,
                      scale = NULL,
                      scale_block = NULL,
                      superblock = NULL,
                      response = NULL,
                      method = NULL,
                      boot = FALSE,
                      inds = NULL) {
  if (is.null(connection)) {
    connection <- rgcca_res$call$connection
  }
  if (is.null(scale)) {
    scale <- rgcca_res$call$scale
  }
  if (is.null(scale_block)) {
    scale_block <- rgcca_res$call$scale_block
  }
  if (is.null(superblock)) {
    superblock <- rgcca_res$call$superblock
  }
  if (is.null(method)) {
    method <- rgcca_res$call$method
  }
  if (is.null(scheme)) {
    scheme <- rgcca_res$call$scheme
  }
  if (is.null(bias)) {
    bias <- rgcca_res$call$bias
  }
  if (is.null(type)) {
    type <- rgcca_res$call$type
  }
  if (is.null(init)) {
    init <- rgcca_res$call$init
  }
  if (is.null(ncomp)) {
    ncomp <- rgcca_res$call$ncomp
  }
  if (is.null(blocks)) {
    #  blocks <- rgcca_res$blocks
    blocks <- rgcca_res$call$raw
    blocks <- descale(blocks)
    if (superblock) {
      for (i in c("tau", "sparsity", "ncomp")) {
        if (class(rgcca_res$call[[i]]) %in% c("matrix", "data.frame")) {
          rgcca_res$call[[i]] <- rgcca_res$call[[i]]
        } else {
          rgcca_res$call[[i]] <- rgcca_res$call[[i]]
        }
      }

      connection <- NULL
    }

    sparsity <- rgcca_res$call$sparsity
    tau <- rgcca_res$call$tau


    if (!is.null(rgcca_res$call$response)) {
      response <- rgcca_res$call$response
    }
  } else
  #        blocks <- scaling(blocks, scale, scale_block = scale_block)

  if (!boot) {
    blocks <- intersection_list(blocks)
  }

  if (tolower(type) %in% c("sgcca", "spca", "spls")) {
    if (!is.null(blocks) && !missing(tau) && missing(sparsity)) {
      stop_rgcca(paste0("sparsity parameter required for ", tolower(type), "instead of tau."))
    }

    par <- "sparsity"
    penalty <- sparsity
  } else {
    if (!is.null(blocks) && !missing(sparsity) && missing(tau)) {
      stop_rgcca(paste0("tau parameter required for ", tolower(type), "instead of sparsity."))
    }

    par <- "tau"
    penalty <- tau
  }

  if (boot) {
    boot_blocks <- list(NULL)
    while (any(sapply(boot_blocks, function(x) length(x)) == 0)) {
      id_boot <- sample(NROW(blocks[[1]]), replace = TRUE)
      boot_blocks <- lapply(
        blocks,
        function(x) {
          y <- x[id_boot, , drop = FALSE]
          rownames(y) <- paste("S", 1:length(id_boot))
          return(y)
        }
      )
      # TODO : to be replaced by something else
      boot_blocks <- remove_null_sd(boot_blocks)
    }
  } else {
    if (length(inds) == 0) {
      boot_blocks <- blocks
    } else {
      boot_blocks <- lapply(blocks, function(x) x[-inds, , drop = FALSE])
      if ("character" %in% class(boot_blocks[[response]])) {
        if (length(unique(boot_blocks[[response]])) == 1) {
          warning("One sample has no variablity. Resulted rgcca can not be run")
          return(NULL)
        }
      }
    }
  }


  func <- quote(
    rgcca(
      boot_blocks,
      connection = connection,
      superblock = superblock,
      response = response,
      ncomp = ncomp,
      scheme = scheme,
      scale = scale,
      scale_block = scale_block,
      type = type,
      verbose = FALSE,
      init = init,
      bias = bias,
      method = method,
      tol = tol
    )
  )

  func[[par]] <- penalty

  res <- eval(as.call(func))
  #  attributes(res)$bigA_scaled <- blocks
  return(res)
}
