var_matrix <- function(q, r, SNR, n, noise_shape) {
  signal <- q %*% (1 / r) %*% (1 / r) %*% t(q)
  p <- diag(nrow(q)) - q %*% t(q)
  tx <- matrix(rnorm(nrow(q)^2), nrow(q))
  tx[upper.tri(tx)] <- 0
  noise1 <- p %*% tcrossprod(tx) %*% p
  noise1 <- noise1 / norm(noise1, type = "F")
  if (noise_shape == "null") {
    noise2 <- 0
  } else {
    noise2 <- p %*% tcrossprod(n) %*% p
    noise2 <- noise2 / norm(noise2, type = "F")
  }
  res <- signal + (noise1 + noise2) * norm(signal, type = "F") / (norm(noise1 + noise2, type = "F") * SNR)
  return(list(Sigma = res, S = signal))
}

shape_getter <- function(root_path, name) {
  if (name == "null") {
    return(matrix(0, 1, 1))
  }
  read.csv(file = paste0(root_path, "/shapes/", name, ".csv"), row.names = NULL, header = FALSE)
}

trim_or_pad <- function(a, b) {
  if (nrow(a) < nrow(b)) {
    b <- b[1:nrow(a), ]
  } else {
    b <- rbind(b, matrix(0, nrow = nrow(a) - nrow(b), ncol = ncol(b)))
  }
  if (ncol(a) < ncol(b)) {
    b <- b[, 1:ncol(a)]
  } else {
    b <- cbind(b, matrix(0, nrow = nrow(b), ncol = ncol(a) - ncol(b)))
  }
  return(b)
}

# Assumes there is only one canonical component per block
generate_data <- function(shapes, noise_shapes, shape_getter,
                          correlation_factors, N, SNR, root_path) {
  # Get canonical components
  J <- length(shapes)
  a <- lapply(shapes, function(x) shape_getter(root_path, x))
  a_noise <- lapply(noise_shapes, function(x) shape_getter(root_path, x))
  a_noise <- lapply(1:J, function(j) trim_or_pad(a[[j]], a_noise[[j]]))
  V <- lapply(a, function(x) matrix(data.matrix(x), ncol = 1))
  V_noise <- lapply(a_noise, function(x) matrix(data.matrix(x), ncol = 1))

  names(a) <- names(V) <- shapes
  names(a_noise) <- names(V_noise) <- noise_shapes

  # Generate covariance matrices
  QR_decompositions <- lapply(V, qr)
  Q <- lapply(QR_decompositions, qr.Q)
  R <- lapply(QR_decompositions, qr.R)
  M <- lapply(correlation_factors, matrix)
  pjs <- sapply(Q, nrow)
  mu <- matrix(0, sum(pjs), 1)
  S <- lapply(1:J, function(j) var_matrix(Q[[j]], R[[j]], SNR, V_noise[[j]], noise_shapes[j]))
  Sigma <- matrix(0, nrow = sum(pjs), ncol = sum(pjs))
  for (i in 0:(J - 1)) {
    for (j in i:(J - 1)) {
      ridx <- seq(1 + sum(pjs[0:i]), sum(pjs[0:(i + 1)]))
      cidx <- seq(1 + sum(pjs[0:j]), sum(pjs[0:(j + 1)]))
      if (i == j) {
        Sigma[ridx, cidx] <- S[[i + 1]]$Sigma
      } else {
        Sigma[ridx, cidx] <- S[[i + 1]]$S %*% V[[i + 1]] %*% M[[i + 1]] %*%
          t(M[[j + 1]]) %*% t(V[[j + 1]]) %*% S[[j + 1]]$S
        Sigma[cidx, ridx] <- t(Sigma[ridx, cidx])
      }
    }
  }

  # Generate data
  x <- mvrnorm(n = N, mu = mu, Sigma = Sigma)
  Xu <- lapply(0:(J - 1), function(j) x[, seq(1 + sum(pjs[0:j]), sum(pjs[0:(j + 1)]))])
  X <- lapply(1:J, function(j) drop(array(data = Xu[[j]], dim = c(N, dim(a[[j]])))))
  return(list(X = X, Xu = Xu, a = a, V = V))
}
