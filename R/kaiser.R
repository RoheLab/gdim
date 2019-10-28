#' Estimate the rank of a matrix with the Kaiser-Guttman Rule
#'
#' That is, set `k` to the number of components of the
#' correlation matrix with eigenvalue greater than one.
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
kaiser <- function(X, ...) {
  # traditionally this is supposed to work on correlation matrices
  as.integer(sum(svd(X)$d) > 1)
}
