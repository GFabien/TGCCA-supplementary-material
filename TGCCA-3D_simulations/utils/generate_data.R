library(abind)

shape_getter <- function(root_path, name) {
  if (name == "null") {
    return(matrix(0, 1, 1))
  }
  file_path <- paste0(root_path, name)
  if (file.exists(paste0(file_path, ".csv"))) {
    x <- read.csv(file = paste0(file_path, ".csv"), row.names = NULL, header = FALSE)
    x <- data.matrix(x)
    dimnames(x) <- NULL
  } else {
    x <- readRDS(paste0(file_path, ".rds"))
  }
  return(x)
}

trim_pad_rep <- function(a, b) {
  dim_a <- dim(a)
  dim_b <- dim(b)
  # If object to resize has more dimensions, unfold it
  if (length(dim_b) > length(dim_a)) {
    b <- array(b, dim = c(dim_b[1:(length(dim_a) - 1)], prod(dim_b[length(dim_a):length(dim_b)])))
    dim_b <- dim(b)
  }
  # If object to resize has less dimensions, repeat it
  if (length(dim_b) < length(dim_a)) {
    b <- array(b, dim = c(dim_b, dim_a[(length(dim_b) + 1):length(dim_a)]))
    dim_b <- dim(b)
  }
  # For all dimensions, pad or trim according to size of a
  for (d in 1:length(dim_a)) {
    if (dim_a[[d]] < dim_b[[d]]) {
      q <- (dim_b[[d]] - dim_a[[d]]) %/% 2
      b <- apply(b, -d, "[", c(q:(q + dim_a[[d]] - 1)))
      if (d > 1) b <- aperm(b, c(2:d, 1, (d + 1):length(dim_b))[1:length(dim_b)])
    }
    if (dim_a[[d]] > dim_b[[d]]) {
      q <- (dim_a[[d]] - dim_b[[d]]) %/% 2
      r <- (dim_a[[d]] - dim_b[[d]]) %% 2
      dim1 <- dim_b
      dim1[[d]] <- q
      dim2 <- dim_b
      dim2[[d]] <- q + r
      b <- abind(array(0, dim = dim1), b, array(0, dim = dim2), along = d)
    }
    dim_b <- dim(b)
  }
  return(b)
}

generate_block <- function(z, cor_factor, shared, unshared, SNR, unshared_shape) {
  N <- length(z)
  p <- NROW(shared)
  shared <- shared / norm(shared, type = "2")
  unshared <- unshared / norm(unshared, type = "2")

  ### Construct deterministic part
  signal <- shared
  x <- sqrt(cor_factor) * signal %*% t(z)

  ### Construct part related to variance
  # Construct part related to the signal (shared component)
  norm_signal <- norm(signal, type = "2")

  # Construct part related to the structured noise (unshared component)
  if (unshared_shape == "null") {
    str_noise <- rep(0, p)
  } else {
    str_noise <- unshared - shared %*% crossprod(shared, unshared)
    str_noise <- str_noise / norm(str_noise, type = "2")
    str_noise <- str_noise %*% t(rnorm(N))
  }

  # Construct part related to the unstructured noise, divide by the frobenius
  # norm of the projector (i.e. sqrt(p - 1))
  unstr_noise <- matrix(rnorm(p * N), p) / sqrt(p - 1)
  unstr_noise <- unstr_noise - shared %*% crossprod(shared, unstr_noise)

  # Add parts together
  x <- x + signal %*% t(rnorm(N)) + norm_signal / SNR * (str_noise + unstr_noise)

  return(t(x))
}

# Assumes there is only one canonical component per block
generate_data <- function(shapes, noise_shapes, shape_getter,
                          correlation_factors, N, SNR,
                          root_path, dims = NULL) {
  # Get canonical components
  J <- length(shapes)
  a <- lapply(shapes, function(x) shape_getter(root_path, x))
  if (!is.null(dims)) {
    a <- lapply(1:length(a), function(j) trim_pad_rep(array(0, dim = dims[[j]]), a[[j]]))
  }
  a_noise <- lapply(noise_shapes, function(x) shape_getter(root_path, x))
  a_noise <- lapply(1:J, function(j) trim_pad_rep(a[[j]], a_noise[[j]]))
  V <- lapply(a, function(x) matrix(data.matrix(x), ncol = 1))
  V_noise <- lapply(a_noise, function(x) matrix(data.matrix(x), ncol = 1))

  names(a) <- names(V) <- shapes
  names(a_noise) <- names(V_noise) <- noise_shapes

  # Generate data
  z <- rnorm(N)
  X <- lapply(1:J, function(j) {
    drop(array(
      generate_block(
        z, correlation_factors[j], V[[j]], V_noise[[j]],
        SNR, noise_shapes[j]
      ),
      dim = c(N, dim(a[[j]]))
    ))
  })
  return(list(X = X, a = a))
}
