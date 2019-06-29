#' Calculate the degree normalization of an adjacency matrix
#'
#' @param A An adjacency matrix. Either a [matrix()] or a [Matrix::Matrix()].
#' @param tau_row Out degree regularization term. Should be positive.
#'   Defaults to `NULL`, in which case we use the mean out degree.
#' @param tau_col In degree regularization term. Should be positive.
#'   Defaults to `NULL`, in which case we use the mean in degree.
#'
#' @return
#' @export
#'
#' @examples
normalize <- function(A, tau_row = NULL, tau_col = NULL) {

  n <- nrow(A)
  d <- ncol(A)

  default_row <- is.null(tau_row)
  default_col <- is.null(tau_col)

  rsA <- Matrix::rowSums(A)  # out-degree
  csA <- Matrix::colSums(A)  # in-degree

  tau_r <- if (default_row) mean(rsA) else tau_row
  tau_c <- if (default_col) mean(csA) else tau_col

  D_row <- Matrix::Diagonal(n = n, x = 1 / sqrt(rsA + tau_r))
  D_col <- Matrix::Diagonal(n = d, x = 1 / sqrt(csA + tau_c))

  # note: no identity matrix in the graph Laplacian here
  D_row %*% A %*% D_col
}
