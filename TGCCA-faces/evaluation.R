### Goal of the analysis:
# From a set of images of individuals under different illumination
# conditions, captured in two different poses, learn a latent representation.
# This representation is then used to match subjects under the two different
# poses under randomly chosen (and potentially different) illumination
# conditions.

### Get working directory from command line
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

library(dplyr)
library(tidyr)
library(RGCCA)
library(abind)
library(imager)
library(ggplot2)
library(viridis)
library(lpSolve)
source("./data_loader.R")

poses <- c(51, 50)
illuminations <- c(seq(2, 6), seq(10, 19))

### Define utility functions
cosine <- function(x, y) {
  -drop(crossprod(x, y)) / (norm(x, type = "2") * norm(y, type = "2"))
}

compute_matching_accuracy <- function(X_test) {
  # Compute distances between projected points in the two modalities of the test set
  dist_matrix <- sapply(seq(nrow(X_test[[1]])), function(i) {
    sapply(seq(nrow(X_test[[2]])), function(j) {
      cosine(X_test[[1]][i, ], X_test[[2]][j, ])
    })
  })

  # Assign pairs that maximize the sum of cosines between pairs
  n <- nrow(dist_matrix)
  constraints <- rbind(
    t(rep(1, n) %x% diag(n)), t(diag(n) %x% rep(1, n))
  )
  directions <- rep("=", nrow(constraints))
  right_hand_side <- rep(1, nrow(constraints))
  matches <- lp("min", as.vector(dist_matrix), constraints, directions,
    right_hand_side,
    int.vec = seq(n^2), all.bin = TRUE
  )$solution
  matches <- matrix(matches, n, n)
  matches <- apply(matches, 2, which.max)

  # Compute accuracy metric
  accuracy <- mean(seq_along(matches) == matches)

  return(accuracy)
}

evaluate_TGCCA <- function(X_test, fit, n_var = 40, test_ill = list(
                             seq(length(illuminations)), seq(length(illuminations))
                           ), projector = "astar") {
  # Apply preprocessing to test blocks
  X_test <- lapply(seq_along(X_test), function(j) {
    aperm(abind(lapply(seq(nrow(X_test[[j]])), function(i) {
      # Create a tensor with illumination slice from the data and others at zero
      y <- array(0, dim = c(dim(X_test[[j]])[seq(2, 3)], length(illuminations)))
      ill <- as.matrix(test_ill[[j]])[i, ]
      x <- array(X_test[[j]][i, , , ], dim = dim(X_test[[j]])[-1])
      DIM <- c(dim(X_test[[j]])[seq(2, 3)], length(illuminations))
      x <- x - array(fit$center[[j]], dim = DIM)[, , ill, drop = FALSE]
      x <- x / array(fit$scaling[[j]], dim = DIM)[, , ill, drop = FALSE]
      y[, , ill] <- x
      return(y)
    }), along = 4), c(4, 1, 2, 3))
  })

  # Project test data
  X_test <- lapply(seq_along(X_test), function(j) {
    x <- matrix(as.vector(X_test[[j]]), nrow = nrow(X_test[[j]]))
    if (projector == "astar") {
      x %*% fit$astar[[j]][, seq(n_var)]
    } else {
      x %*% fit$a[[j]][, seq(n_var)]
    }
  })

  compute_matching_accuracy(X_test)
}

##### Training on the whole set of training images
individuals <- seq(1, 100)
X_train <- lapply(poses, function(pose) {
  aperm(drop(abind(
    lapply(individuals, function(ind) {
      abind(lapply(
        illuminations, function(illumination) {
          data_loader(ind, pose, illumination)
        }
      ), along = 3)
    }),
    along = 4
  )), c(4, 1, 2, 3))
})

# Train TGCCA1
file_name <- "./TGCCA1"
if (!file.exists(file_name)) {
  fit <- rgcca(
    blocks = X_train, type = "mgcca", ranks = 1,
    tau = 1, scale = FALSE, scale_block = TRUE,
    verbose = TRUE, scheme = "horst", ncomp = 25,
    init = "random", tol = 1e-4
  )
  saveRDS(list(
    a = fit$a,
    astar = fit$astar,
    weights = fit$weights,
    factors = fit$factors,
    center = lapply(fit$blocks, attr, "scaled:center"),
    scaling = lapply(fit$blocks, attr, "scaled:scale")
  ), file = file_name)
  rm(fit)
}

# Train TGCCA3
file_name <- "./TGCCA3"
if (!file.exists(file_name)) {
  fit <- rgcca(
    blocks = X_train, type = "mgcca", ranks = 3,
    tau = 1, scale = FALSE, scale_block = TRUE,
    verbose = TRUE, scheme = "horst", ncomp = 25,
    init = "random", tol = 1e-4
  )
  saveRDS(list(
    a = fit$a,
    astar = fit$astar,
    weights = fit$weights,
    factors = fit$factors,
    center = lapply(fit$blocks, attr, "scaled:center"),
    scaling = lapply(fit$blocks, attr, "scaled:scale")
  ), file = file_name)
  rm(fit)
}

# Train RGCCA
file_name <- "./RGCCA"
if (!file.exists(file_name)) {
  fit <- rgcca(
    blocks = lapply(X_train, function(x) matrix(as.vector(x), nrow = nrow(x))),
    type = "rgcca", tau = 1, scale = FALSE, scale_block = TRUE,
    verbose = TRUE, scheme = "horst", ncomp = 25,
    init = "random", tol = 1e-4
  )
  saveRDS(list(
    a = fit$a,
    astar = fit$astar,
    center = lapply(fit$blocks, attr, "scaled:center"),
    scaling = lapply(fit$blocks, attr, "scaled:scale")
  ), file = file_name)
  rm(fit)
}

# Train LDA model
X_lda_train <- data.frame(do.call(rbind, lapply(X_train, function(x) {
  x <- aperm(x, c(4, 1, 2, 3))
  x <- matrix(as.vector(x), nrow = dim(x)[1] * dim(x)[2])
  x <- t(apply(x, 1, function(y) {
    as.vector(resize(as.cimg(matrix(y, 64)), 16, 16, interpolation_type = 3))
  }))
})))
y_lda_train <- as.factor(rep(
  rep(seq_along(illuminations), length(individuals)), 2
))
fit_lda <- lda(y_lda_train ~ ., data = X_lda_train)

##### Evaluation
# Sample illuminations
set.seed(0)
sample_illuminations <- function(n_ill = 1, seed = 0) {
  set.seed(seed)
  # For each view, sample n_ill illumination condtions per subject
  lapply(seq(2), function(j) {
    x <- vapply(seq(nrow(X_train[[j]])), function(i) {
      sample(length(illuminations), size = n_ill, replace = FALSE)
    }, FUN.VALUE = integer(n_ill))
    if (n_ill == 1) {
      x <- as.matrix(x)
    } else {
      x <- t(x)
    }
    return(x)
  })
}

evaluate_models <- function(X, n_ill, n_var = 10, seed = 0) {
  # Sample illumination conditions
  test_illuminations <- sample_illuminations(n_ill, seed)
  # Remove images that are not sampled
  if (n_ill == 1) {
    X_test <- lapply(seq_along(X), function(j) {
      x <- aperm(abind(lapply(seq(nrow(X[[j]])), function(k) {
        X[[j]][k, , , test_illuminations[[j]][k, ]]
      }), along = 3), c(3, 1, 2))
      x <- array(x, dim = c(dim(x), 1))
    })
  } else {
    X_test <- lapply(seq_along(X), function(j) {
      abind(lapply(seq(nrow(X[[j]])), function(k) {
        X[[j]][k, , , test_illuminations[[j]][k, ], drop = FALSE]
      }), along = 1)
    })
  }

  # Predict illuminations using the LDA classifier
  X_test_lda <- lapply(X_test, function(x) {
    do.call(rbind, lapply(seq(nrow(x)), function(k) {
      t(vapply(seq(n_ill), function(ill) {
        imager::resize(as.cimg(x[k, , , ill]), 16, 16, interpolation_type = 3)
      }, FUN.VALUE = double(256L)))
    }))
  })

  predictions <- lapply(X_test_lda, function(x) {
    y <- predict(fit_lda, data.frame(x))$class
    matrix(as.integer(y), ncol = n_ill, byrow = TRUE)
  })

  # Load fit and evaluate accuracy
  accuracy <- lapply(models, function(model) {
    fit <- readRDS(paste0("./", model))
    evaluate_TGCCA(X_test, fit, n_var = n_var, test_ill = predictions)
  })
  names(accuracy) <- models
  return(accuracy)
}

### Repeat experiments several times to draw confidence intervals
# Load test data
test_individuals <- seq(101, 200)
X_test <- lapply(poses, function(pose) {
  aperm(drop(abind(
    lapply(test_individuals, function(ind) {
      abind(lapply(
        illuminations, function(illumination) {
          data_loader(ind, pose, illumination)
        }
      ), along = 3)
    }),
    along = 4
  )), c(4, 1, 2, 3))
})

models <- c("TGCCA1", "TGCCA3", "RGCCA")

# Run 100 times the experiment
file_name <- "./results.rds"
if (file.exists(file_name)) {
  complete_df <- readRDS(file_name)
} else {
  n_rep <- 5
  n_vars <- c(10, 15, 20, 25)
  complete_df <- c()
  for (n in seq(n_rep)) {
    print(paste0("Repetition ", n))
    for (n_var in n_vars) {
      results_test <- lapply(models, function(x) NULL)
      names(results_test) <- models
      for (n_ill in seq(length(illuminations))) {
        accuracy <- evaluate_models(X_test, n_ill, n_var = n_var, seed = n)

        for (m in models) {
          results_test[[m]] <- c(results_test[[m]], accuracy[[m]])
        }
      }

      df <- data.frame(do.call(cbind, results_test))
      df$n_var <- paste(n_var, "components")
      df$illuminations <- seq_along(illuminations)
      df$fold <- n
      complete_df <- rbind(complete_df, df)
    }
  }

  saveRDS(complete_df, file_name)
}

complete_df <- complete_df %>%
  pivot_longer(cols = c("RGCCA", "TGCCA1", "TGCCA3"))

summary_df <- complete_df %>%
  group_by(n_var, illuminations, name) %>%
  summarise(across(
    .cols = c("value"),
    .fns = list(
      median = median,
      q025 = function(x) quantile(x, 0.025),
      q975 = function(x) quantile(x, 0.975)
    ),
    .names = "{fn}"
  ))

p <- ggplot(data = summary_df, aes(
  x = illuminations, y = median, color = name, linetype = name
)) +
  geom_line(size = 1.5) +
  geom_ribbon(
    data = summary_df, aes(ymin = q025, ymax = q975, fill = name),
    alpha = 0.2
  ) +
  scale_color_manual(values = cividis(length(models))) +
  scale_fill_manual(values = cividis(length(models))) +
  facet_wrap(~n_var, strip.position = "top", ncol = 4) +
  theme_light() +
  labs(color = "Model", linetype = "Model", fill = "Model") +
  xlab("Number of available illumination conditions") +
  ylab("Matching accuracy") +
  theme(
    strip.background = element_blank(), strip.placement = "outside",
    strip.text = element_text(colour = "black"), text = element_text(size = 12)
  )
ggsave(filename = "./accuracy_vs_illumination.pdf", plot = p)
