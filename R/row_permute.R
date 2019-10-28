#' Estimate the rank of a graph adjacency matrix with row permutation
#'
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

# undirected graphs

alpha = 0.05

summary(A)

nnzero(A)
A[1:5]
