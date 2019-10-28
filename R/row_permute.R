#' Estimate the rank of a graph adjacency matrix with row permutation
#'
#' @param X A matrix. Sparse matrices are automatically
#'   coerces to dense matrices, so be careful.
#'
#' @return The estimated rank of the matrix, as an integer.
#' @export
#'
#' @examples
#'
#' set.seed(27)
#'
#' n <- 1000
#' k <- 5
#'
#' B <- matrix(runif(k * k), nrow = k, ncol = k) # mixing probabilities
#'
#' theta <- round(rlnorm(n, 2)) # degree parameter
#' pi <- c(1, 2, 4, 1, 1) / 5 # community membership
#'
#' A <- dcsbm(theta, pi, B, avg_deg = 50)
#' kaiser(A)
#'
#' X <- iris[, 1:4]
#' kaiser(X)
#'
row_permute <- function(X, k_max = min(ncol(X), 1000), ...) {

  # traditionally this is supposed to work on correlation matrices
  as.integer(sum(svd(X)$d) > 1)
}

#' Sampling distribution of row permutation test lambda
#'
#' @param X
#' @param k
#' @param B
#'
#' @return
#' @export
#'
#' @examples
#'
#' library(Matrix)
#' library(tidyverse)
#'
#' A <- rsparsematrix(1000, 1000, 0.1)
#'
#' # symmetrize
#' A <- sign(t(A) + A)
#'
#' sampling_dist <- rp_distribution(A, 10, B = 200)
#'
#' sampling_dist %>%
#'   as_tibble() %>%
#'   mutate(B = row_number()) %>%
#'   gather(lambda, value, -B) %>%
#'   ggplot(aes(value)) +
#'   geom_histogram() +
#'   facet_wrap(~lambda) +
#'   theme_minimal()
#'
rp_distribution <- function(X, k, B = 1000) {

  stopifnot(isSymmetric(X))

  splits <- split_edges(X, symmetric = TRUE)
  eig_train <- RSpectra::eigs(splits$train, k = k)

  lambda <- matrix(0, nrow = B, ncol = k)
  colnames(lambda) <- paste0("lambda", 1:k)

  for (i in 1:B) {

    # irritating transpose in apply, but avoid correcting with
    # t() because that is expensive
    U <- apply(eig_train$vectors, 1, sample)

    for (j in 1:k) {
      quad_form <- t(U[j, ]) %*% splits$test %*% U[j, ]
      lambda[i, j] <- quad_form[1, 1]
    }

  }

  lambda
}

