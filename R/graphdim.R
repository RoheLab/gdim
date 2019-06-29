#' Estimate the dimension of a graph
#'
#' @param x A graph.
#' @param ... Unused.
#'
#' @return
#' @export
#'
#' @examples
graphdim <- function(x, ...) {
  UseMethod("graphdim")
}


graphdim.default <- function(x, ...) {

  # step 1: build A_train and A_test

  # step 2: optionally normalize both A_train and A_test

  # step 3: compute the eigenvectors of A_train

  # TODO: what to do for directed graphs? SVD instead?
  # https://see.stanford.edu/materials/lsoeldsee263/15-symm.pdf

  # step 4: compute the bootstrapped eigenvalues


}

graphdim.igrph <- function(x, ...) {
  x <- igraph::get.adjacency(x, sparse = TRUE)
  NextMethod()
}

# add horn's method for comparison. paran::paran() does PCA
# but always on a correlation matrix, which is not what you want
# in the case of adjacency matrices IIRC
