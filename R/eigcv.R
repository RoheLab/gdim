
# assumes elements of A are non-negative integers (and Poisson)
split_graph <- function(A, test_portion = 0.1) {

  A <- methods::as(A, "dMatrix")
  A <- methods::as(A, "TsparseMatrix")

  stopifnot(isTRUE(all.equal(A@x, as.integer(A@x))))
  test_edges <- stats::rbinom(length(A@x), size = A@x, prob = test_portion)

  train <- A
  test <- A

  # use as.numeric() to match `dMatrix` data type
  train@x <- as.numeric(train@x - test_edges)
  test@x <- as.numeric(test_edges)

  train <- methods::as(drop0(train), "dgCMatrix")
  test <- methods::as(drop0(test), "dgCMatrix")

  list(train = train, test = test)
}

glaplacian <- function(A, regularize = TRUE) {
  deg_row <- Matrix::rowSums(A)
  deg_col <- Matrix::colSums(A)
  if (regularize) {
    tau_row <- mean(deg_row)
    tau_col <- mean(deg_col)
  } else {
    if (any(c(deg_row, deg_row) == 0)) {
      stop(
        "Cannot use Laplacian because some nodes are isolated. ",
        "Set either `regularize=TRUE` or `laplacian=FALSE` option."
      )
    }
    tau_row <- tau_col <- 0
  }
  D_row <- Diagonal(nrow(A), 1 / sqrt(deg_row + tau_row))
  D_col <- Diagonal(ncol(A), 1 / sqrt(deg_col + tau_col))
  L <- D_row %*% A %*% D_col
  L
}

# given the training singular vectors, compute the test statistic for
# graph dimension.
gdstat <- function(full, test, u, v, test_portion) {

  # standard error calculation is different for directed and undirected graphs
  if (isSymmetric(full)) {
    se <- sqrt(2 * test_portion * as.numeric(t(u^2) %*% full %*% v^2) -
      test_portion * sum(diag(full) * u^2 * v^2))
  } ## standard error
  if (!isSymmetric(full)) {
    se <- sqrt(test_portion * as.numeric(t(u^2) %*% full %*% v^2))
  }
  lamL <- as.numeric(t(u) %*% glaplacian(test / test_portion) %*% v)
  lamA <- as.numeric(t(u) %*% test %*% v) / test_portion
  z <- as.numeric(t(u) %*% test %*% v) / se ## test stat
  c(cv_lambda_A = lamA, cv_lambda_L = lamL, z = z)
}



#' Compute cross-validate eigenvalues
#'
#' Estimate graph dimension via eigenvalue cross-validation (EigCV).
#' A graph has dimension `k` if the first `k` eigenvectors of its adjacency
#' matrix are correlated with its population eigenspace, and the others are not.
#' Edge bootstrapping sub-samples the edges of the graph (without replacement).
#' Edge splitting separates the edges into a training part and a testing part.
#'
#' @param A The adjacency matrix of graph. Must be non-negative and
#'  integer valued.
#' @param k_max The maximum dimension of the graph to consider. This many
#'   eigenvectors are computed. Should be a non-negative integer smallish
#'   relative the dimensions of `A`.
#' @param ... Ignored.
#' @param num_bootstraps The number of times to bootstrap the graph. Since
#'   cross-validated eigenvalues are based on a random graph split, they
#'   are themselves random. By repeatedly computing cross-validated eigenvalues
#'   for different sample splits, the idea is to smooth away some of the
#'   randomness due to the graph splits. A small number of bootstraps
#'   (3 to 10) usually suffices. Defaults to `10`. Test statistics (i.e.
#'   z-scores for cv eigenvalues) are averaged across bootstraps
#'   and the p-values will be calculated based on the averaged statistics.
#' @param test_portion The portion of the graph to put into the test graph,
#'   as opposed to the training graph. Defaults to `0.1`. Must be strictly
#'   between zero and one.
#' @param alpha Significance level for hypothesis tests. Each dimension
#'   `1, ..., k_max` is tested when estimating graph dimension, and the
#'   overall graph dimension is taken to be the smallest number of dimensions
#'   such that all the tests reject.
#' @param method Method to adjust p-values for multiple testing. Must be
#'   one of `"none"`, `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`,
#'   `"BH"`, `"BY"`, or `"fdr"`. Passed to [stats::p.adjust()]. Defaults to
#'   `"none"`.
#' @param laplacian Logical value indicating where to compute cross-validated
#'   eigenvalues for the degree-normalize graph Laplacian rather than the
#'   graph adjacency matrix. Experimental and should be used with caution.
#'   Defaults to `FALSE`.
#' @param regularize Only applicable when `laplacian == TRUE`, in which case
#'   this parameter controls whether or not the degree-normalized graph
#'   Laplacian is regularized. Defaults to `TRUE`.
#'
#' @return A `eigcv` object, which is a list with the following named
#'   elements.
#'
#'   - `estimated_dimension`: inferred graph dimension.
#'   - `summary`: summary table of the tests.
#'   - `num_bootstraps`: number of bootstraps performed.
#'   - `test_portion`: graph splitting probability used.
#'   - `alpha`: significance level of each test.
#'
#' @export
#'
#' @examples
#'
#' library(fastRG)
#'
#' set.seed(27)
#'
#' B <- matrix(0.1, 5, 5)
#' diag(B) <- 0.3
#'
#' model <- sbm(
#'   n = 1000,
#'   k = 5,
#'   B = B,
#'   expected_degree = 40,
#'   poisson_edges = FALSE,
#'   allow_self_loops = FALSE
#' )
#'
#' A <- sample_sparse(model)
#'
#' eigs<- eigcv(A, k_max = 10)
#' eigs
#'
#' plot(eigs, type = "z-score")    # default
#' plot(eigs, type = "adjacency")
#' plot(eigs, type = "laplacian")
#'
#'
eigcv <- function(A, k_max,
                  ...,
                  num_bootstraps = 10, test_portion = 0.1,
                  alpha = 0.05,
                  method = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                  laplacian = FALSE,
                  regularize = TRUE) {
  n <- min(dim(A))
  stopifnot("`k_max` is too large." = k_max <= n)
  stopifnot("`test_portion` must range between 0 and 1." = (test_portion > 0 && test_portion < 1))
  stopifnot("`num_bootstraps` must be a positive integer." = num_bootstraps >= 1)

  method <- rlang::arg_match(method)

  pb <- progress::progress_bar$new(
    total = num_bootstraps,
    format = "bootstrapping [:bar] :current/:total (:percent) eta: :eta"
  )
  pb$tick(0)

  ## full graph
  full <- A <- A * 1
  if (laplacian) {
    full <- glaplacian(A, regularize = regularize)
  }

  cv_stats <- tibble::tibble(
    .rows = k_max * num_bootstraps,
    boot = 0,
    k = 0,
    cv_lambda_A = 0,
    cv_lambda_L = 0,
    z = 0
  )

  tick <- 0
  for (boot in 1:num_bootstraps) {

    # must split on adjacency matrix A, even if computing cv-eigs for
    # graph Laplacian

    es <- split_graph(A, test_portion)
    train <- es$train
    test <- es$test

    if (laplacian) {
      train <- glaplacian(es$train, regularize = regularize)
    }

    # train should be dgcMatrix class for max computational efficiency
    s_train <- irlba::irlba(train, k_max)

    for (k in 1:k_max) {
      tick <- tick + 1
      gds <- gdstat(full = A, test = test, u = s_train$u[, k], v = s_train$v[, k], test_portion = test_portion)
      cv_stats[tick, 1:5] <- matrix(c(boot, k, gds), nrow = 1)
    }

    pb$tick()
  }

  cv_means <- cv_stats %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(
      cv_lambda_A = mean(cv_lambda_A),
      cv_lambda_L = mean(cv_lambda_L),
      z = mean(z)
    ) %>%
    dplyr::ungroup()


  cv_means <- cv_means %>%
    dplyr::mutate(
      pvals = stats::pnorm(z, lower.tail = FALSE)
    )

  cv_means$padj <- stats::p.adjust(cv_means$pvals, method = method)

  criteria <- cv_means$padj
  k_stop <- which(criteria > alpha)
  k_infer <- ifelse(length(k_stop), min(k_stop) - 1, k_max)
  res <- list(
    estimated_dimension = k_infer,
    summary = cv_means,
    num_bootstraps = num_bootstraps,
    test_portion = test_portion,
    alpha = alpha
  )
  res$stats <- cv_stats
  class(res) <- "eigcv"
  res
}


#' Print cross-validated eigenvalues
#'
#' @param x An `eigcv` object created by a call to [eigcv()].
#' @param ... Ignored.
#' @export
#'
#' @inherit eigcv examples
#' @return `x`, but invisibly.
#'
#' @method print eigcv
print.eigcv <- function(x,  ...) {
  cat("Estimated graph dimension:\t", x$estimated_dimension, fill = TRUE)
  cat("\nNumber of bootstraps:\t\t", x$num_bootstraps, fill = TRUE)
  cat("Edge splitting probabaility:\t", x$test_portion, fill = TRUE)
  cat("Significance level:\t\t", x$alpha, fill = TRUE)
  cat("\n ------------ Summary of Tests ------------\n")
  print(data.frame(x$summary[, -c(2, 3)]), row.names = FALSE)
  cat(fill = TRUE)
  invisible(x)
}

#' Plot cross-validated eigenvalues
#'
#' @param x An `eigcv` object created by a call to [eigcv()].
#' @param type Specifies what to plot. Must be one of the following options:
#'
#'    - `"z-score"`, in which case the Z-statistic test scores are plotted
#'      for each value of `k` (i.e. dimension of the eigenspace).
#'    - `"adjacency"` in which case the cross-validated eigenvalues of the
#'      adjacency matrix are plotted for each value of `k`.
#'    - `"laplacian"` in which case the cross-validated eigenvalues of the
#'      graph Laplacian matrix are plotted for each value of `k`.
#'
#' @param threshold Only used when `type == "z-score"`. Adds a horizontal
#'   line at the value of `threshold`, which should be a numeric of length
#'   one. Defaults to `2`.
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @export
#' @method plot eigcv
#'
#' @inherit eigcv examples
plot.eigcv <- function(x, type = c("z-score", "adjacency", "laplacian"), threshold = 2, ...) {
  stopifnot("Threshold of statistics must be greater than 0." = threshold > 0)

  type <- rlang::arg_match(type)

  type <- type[1]
  if (type == "z-score") {
    dat <- dplyr::select(x$summary, k, val = z)
    ylab <- "z score"
  }
  if (type == "adjacency") {
    dat <- dplyr::select(x$summary, k, val = cv_lambda_A)
    ylab <- "cross validated x' A x"
  }
  if (type == "laplacian") {
    dat <- dplyr::select(x$summary, k, val = cv_lambda_L)
    ylab <- "cross validated x' L x"
  }

  g <- ggplot2::ggplot(ggplot2::aes(k, val), data = dat) +
    ggplot2::geom_point(alpha = .8) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      breaks = function(x) {
        unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
      }
    )

  if (type == "z-score") {
    g <- g +
      ggplot2::geom_hline(
        yintercept = threshold, alpha = .8,
        linetype = 2, color = "grey60", show.legend = TRUE
      )
  }

  g
}
