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
#' sampling_dist <- sf_distribution(A, 10, B = 200)
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
sf_distribution <- function(X, k, B = 1000) {

  stopifnot(isSymmetric(X))

  splits <- split_edges(X, symmetric = TRUE)
  eig_train <- RSpectra::eigs(splits$train, k = k)

  lambda <- matrix(0, nrow = B, ncol = k)
  colnames(lambda) <- paste0("lambda", 1:k)

  for (i in 1:B) {

    flips <- matrix(extraDistr::rsign(nrow(X) * k), nrow = nrow(X), ncol = k)

    # independently flip signs of eigenvector elements
    U <- eig_train$vectors * flips

    for (j in 1:k) {
      quad_form <- t(U[, j]) %*% splits$test %*% U[, j]
      lambda[i, j] <- quad_form[1, 1]
    }

  }

  lambda
}

