parse_param <- function(param, type = NULL) {
  if (type == "bool") {
    return(param == "T" || param == "TRUE")
  }
  if (type == "num_vec") {
    return(as.numeric(strsplit(param, ",")[[1]]))
  }
  return(param)
}

run_tgcca <- function(dataset, raw_params = NULL, root_path = NULL, seed = 0, random_init = FALSE) {
  params <- list(ranks = 1, tau = 0.001, scale_block = T, scheme = "horst", scale = F, n_run = 10, n_cores = 1)
  mapping <- c("ranks", "tau", "scale_block", "scheme", "scale", "n_run", "n_cores")
  param_type <- c("num_vec", "num_vec", "bool", "string", "bool", "num_vec", "num_vec")
  if (!is.null(raw_params)) {
    for (i in 1:length(raw_params)) {
      params[[mapping[i]]] <- parse_param(raw_params[i], param_type[i])
    }
  }
  params[["blocks"]] <- dataset$X
  params[["type"]] <- "ns_mgcca"

  fit_mgcca <- do.call(rgcca, params)
  return(list(a = fit_mgcca$astar))
}

run_sptgcca <- function(dataset, raw_params = NULL, root_path = NULL, seed = 0, random_init = FALSE) {
  params <- list(ranks = 1, tau = 0.001, scale_block = T, scheme = "horst", scale = F, n_run = 10, n_cores = 1)
  mapping <- c("ranks", "tau", "scale_block", "scheme", "scale", "n_run", "n_cores")
  param_type <- c("num_vec", "num_vec", "bool", "string", "bool", "num_vec", "num_vec")
  if (!is.null(raw_params)) {
    for (i in 1:length(raw_params)) {
      params[[mapping[i]]] <- parse_param(raw_params[i], param_type[i])
    }
  }
  params[["blocks"]] <- dataset$X
  params[["type"]] <- "mgcca"

  fit_mgcca <- do.call(rgcca, params)
  return(list(a = fit_mgcca$astar, weights = fit_mgcca$weights))
}

run_rgcca <- function(dataset, raw_params = NULL, root_path = NULL, seed = 0) {
  params <- list(tau = 0.001, scale_block = T, scheme = "horst", scale = F, n_run = 10, n_cores = 1)
  mapping <- c("tau", "scale_block", "scheme", "scale", "n_run", "n_cores")
  param_type <- c("num_vec", "bool", "string", "bool", "num_vec", "num_vec")
  if (!is.null(raw_params)) {
    for (i in 1:length(raw_params)) {
      params[[mapping[i]]] <- parse_param(raw_params[i], param_type[i])
    }
  }
  params[["blocks"]] <- lapply(dataset$X, function(x) {
    matrix(as.vector(x), nrow = nrow(x))
  })
  params[["type"]] <- "rgcca"
  fit_mgcca <- do.call(rgcca, params)
  return(list(a = fit_mgcca$astar))
}

run_svd <- function(dataset, raw_params = NULL, root_path = NULL, seed = 0) {
  a <- lapply(lapply(dataset$X, function(x) {
    matrix(as.vector(x), nrow = nrow(x))
  }), function(x) svd(x, nv = 1)$v)
  return(list(a = a))
}

run_tcca <- function(dataset, raw_params = NULL, root_path = NULL, seed = 0, random_init = FALSE) {
  require(R.matlab)

  params <- list(ranks = c(1, 1), tau = 0.001, scale_block = T, scale = F, covstr = "none", n_run = 10, compute_score = FALSE)
  mapping <- c("ranks", "tau", "scale_block", "scale", "covstr", "n_run", "compute_score")
  param_type <- c("num_vec", "num_vec", "bool", "bool", "string", "num_vec", "bool")
  if (!is.null(raw_params)) {
    for (i in 1:length(raw_params)) {
      params[[mapping[i]]] <- parse_param(raw_params[i], param_type[i])
    }
  }

  if (random_init) {
    params$n_run <- 1
  }

  blocks <- scaling(dataset$X, scale_block = params$scale_block, scale = params$scale)
  mx <- length(dim(blocks[[1]]))
  X <- aperm(blocks[[1]], c(2:mx, 1))
  my <- length(dim(blocks[[2]]))
  Y <- aperm(blocks[[2]], c(2:my, 1))

  output_dir <- paste0(
    "/tmp/R/MGCCA/",
    gsub(pattern = "  | |:", replacement = "_", x = Sys.time()),
    paste0(sample(LETTERS, 15, TRUE), collapse = "")
  )
  dir.create(file.path(output_dir), showWarnings = TRUE, recursive = TRUE)
  writeMat(paste0(output_dir, "/params.mat"),
    X = X,
    Y = Y,
    rx = params$ranks[1],
    ry = params$ranks[2],
    tau = params$tau,
    covstr = params$covstr,
    replicates = params$n_run,
    seed = seed,
    compute_score = params$compute_score
  )

  script_path <- paste0(root_path, "/TCCA-minchizhou/run_tcca.m")
  cmd <- paste0("matlab  -nodisplay -nosplash -nodesktop -r \"data_path='", output_dir, "';run('", script_path, "')\"")
  print(cmd)
  system(paste0("matlab  -nodisplay -nosplash -nodesktop -r \"data_path='", output_dir, "';run('", script_path, "')\""), wait = TRUE)

  a <- list()
  res <- readMat(paste0(output_dir, "/results.mat"))
  a[[1]] <- res$betax
  a[[2]] <- res$betay

  ### Remove tmp files
  unlink(output_dir, recursive = T)

  return(list(a = a, scorex = res$scorex, scorey = res$scorey))
}
