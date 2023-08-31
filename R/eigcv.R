#' Graph Splitting
#'
#' Randomly split edges of graph into two, proportional to `1-split` and `split`.
#'
#' @param A the adjacency matrix, will be coerced to a `Matrix` class.
#' @param split `numeric(1)`, the proportion of leave-out edges (into the test set).
#'   Must take value from (0,1).
#' @return A list of two `Matrix` elements:
#' \item{train}{the subgraph with `1-split` portion of edges.}
#' \item{test}{the subgraph with `split` portion of edges.}
# @importFrom  Matrix sparseMatrix summary drop0 isTriangular
#' @importFrom stats rbinom runif p.adjust pnorm
esplit <- function(A, split = 0.1) {
  el <- Matrix::summary(A) ## edge list
  if (is.null(el$x)) {
    el$x <- 1
  }
  isSymmetric <- Matrix::isSymmetric(A)
  if (isSymmetric) {
    ## simplify edge list; this implies A is a square matrix
    el <- el[el$i <= el$j, ]
  }
  # stopifnot("A contains negative/fractional" = all(el$x %% 1 == 0 & el$x >= 0))

  ## splitting
  # leave <- (runif(nrow(el)) < split) * 1
  stopifnot(all.equal(el$x, as.integer(el$x)))
  test_edges <- rbinom(nrow(el), size = el$x, prob = split)

  train <- sparseMatrix(
    i = el$i,
    j = el$j,
    x = el$x - test_edges,
    # x = el$x * (1 - leave),
    dims = dim(A),
    dimnames = dimnames(A),
    symmetric = isSymmetric,
    triangular = isTriangular(A)
  )

  test <- sparseMatrix(
    i = el$i,
    j = el$j,
    x = test_edges,
    # x = el$x * leave,
    dims = dim(A),
    dimnames = dimnames(A),
    symmetric = isSymmetric,
    triangular = isTriangular(A)
  )

  list(train = drop0(train), test = drop0(test))
}

#' Graph Laplacian
#'
#' Given an `m` by `n` graph adjacency matrix `A`, calculate the symmetric
#' graph Laplacian.
#'
#' @inheritParams esplit
#' @param regularize `logical(1)`, return a regularized symmetric Laplacian.
#'   The default is `TRUE`. This is ignored if `laplacian=FALSE`.
# @importFrom Matrix %*% Diagonal rowSums colSums
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
  return(L)
}

#' Graph Spectral Decomposition
#'
#' Given ab `m` by `n` graph matrix `A`, find the largest `k` singular values
#' and the corresponding singular vectors.
#'
#' @inheritParams esplit
#' @param k `integer(1)`, number of eigenvalues.
#' @param ... additonal parameters to pass into [RSpectra::eigs] (for symmetric
#'   `A`) or [RSpectra::svds] (for asymmetric `A`).
#' @importFrom RSpectra eigs svds
#' @return a list (see also [RSpectra::svds]):
#'   \item{u}{an `m` by `k` matrix whose columns contain the left singular vectors.}
#'   \item{v}{an `n` by `k` matrix whose columns contain the right singular vectors.}
gspectral <- function(A, k, ...) {
  stopifnot("k too large" = k <= min(dim(A)))
  is_sym <- Matrix::isSymmetric(A)
  if (is_sym) {
    ei <- RSpectra::eigs(A, k, which = "LR", ...)
    return(list(
      u = ei$vect,
      v = ei$vect,
      d = ei$val
    ))
  }
  if (!is_sym) {
    S <- RSpectra::svds(A, k, ...)
  }
  return(list(u = S$u, v = S$v, d = S$d))
}


#' Graph Dimension Statistic
#'
#' Given the trained left/right singular vectors, compute the test statistic for
#' graph dimension.
#'
#' @param full,test `matrix` or `Matrix`, the adjacency or Laplacian matrix of
#'   the full and test graphs.
#' @param u,v `numeric` vector, the trained left and right singular vectors.
#' @inheritParams esplit
#' @return `numeric(3)`, test statistics
gdstat <- function(full, test, u, v, split) {
  if (isSymmetric(full)) {
    se <- sqrt(2 * split * as.numeric(t(u^2) %*% full %*% v^2) -
      split * sum(diag(full) * u^2 * v^2))
  } ## standard error
  if (!isSymmetric(full)) {
    se <- sqrt(split * as.numeric(t(u^2) %*% full %*% v^2))
  }
  lamL <- as.numeric(t(u) %*% glaplacian(test / split) %*% v)
  lamA <- as.numeric(t(u) %*% test %*% v) / split
  z <- as.numeric(t(u) %*% test %*% v) / se ## test stat
  return(c(cv_lambda_A = lamA, cv_lambda_L = lamL, z = z))
}



#' Edge Bootstrapping and Splitting
#'
#' Estimate the graph dimension via eigenvalue cross-validation (EigCV).
#' A graph has dimension `k` if the first `k` eigenvectors of its adjacency
#' matrix are correlated with its population eigenspace, and the others are not.
#' Edge bootstrapping sub-samples the edges of the graph (without replacement).
#' Edge splitting separates the edges into a training part and a testing part.
#'
#' @inheritParams esplit
#' @inheritParams glaplacian
#' @inheritParams gspectral
#' @param k_max `integer(1)`, number of eigenvectors to compute.
#' @param bootstrap `integer(1)`, number of graph bootstraps, default to 10.
#'   Graph bootstrapping is to account for the randomness in graph splitting,
#'   rather than obtaining any statistic (as a traditional bootstrap does).
#'   Hence, a small number (e.g., 3~10) of bootstraps usually suffices.
#'   If `bootstrap>1`, the test statistics will be averaged across bootstraps
#'   and the p-values will be calculated based on the averaged statistics.
#' @param alpha `numeric(1)`, significance level of each test, default to 0.05.
#'   This is used to cut off the dimension estimation.
#' @param ptol `numeric(1)`, the tolerance of minimal p-value.
#' @param correct `character(1)`, correction method for multiple testing.
#'   This is recommended if `k_max` is large. See [p.adjust] for the acceptable
#'   values of `method`. If `correct` is not "none", the returned summary table
#'   will have an extra column called `padj` that repeats the p-value in the
#'   previous column.
#' @param laplacian `logical(1)`, use the normalized and regularized adjacency
#'   matrix (i.e. L)
#'   This option is experimental and should be used with caution.
#' @param trace `logical(1)`, for diagnostic use. If `TRUE`, return additionally
#'   `cv_stats` containing the test statistics of every bootstrap, and if
#'   `n=TRUE`, then the eigs of the full graph.
#' @return A `eigcv` object, which contains:
#'   \item{inference}{inferred graph dimension.}
#'   \item{summary}{summary table of the tests.}
#'   \item{bootstrap}{number of bootstraps performed.}
#'   \item{split}{graph splitting probability used.}
#'   \item{alpha}{significance level of each test.}
#'   \item{trace}{only if `trace=TRUE`, all bootstrapped statistics.}
#' @importFrom tibble tibble
#' @importFrom dplyr summarize group_by ungroup mutate summarise
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
eigcv <- function(A, k_max,
                  bootstrap = 10, split = 0.1,
                  alpha = 0.05,
                  ptol = .Machine$double.eps,
                  correct = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
                              "none"),
                  laplacian = TRUE,
                  regularize = TRUE,
                  trace = FALSE) {
  n <- min(dim(A))
  stopifnot("`k_max` is too large." = k_max <= n)
  stopifnot("`split` must range between 0 and 1." = (split > 0 && split < 1))
  stopifnot("`bootstrap` must be a positive integer." = bootstrap >= 1)

  ## full graph
  full <- A <- A * 1
  if (laplacian) {
    full <- glaplacian(A, regularize = regularize)
  }

  cv_stats <- tibble::tibble(
    .rows = k_max * bootstrap,
    boot = 0,
    k = 0,
    cv_lambda_A = 0,
    cv_lambda_L = 0,
    z = 0
  )
  tick <- 0

  for (boot in 1:bootstrap) {
    ## edge splitting
    es <- esplit(A, split)
    train <- es$train
    test <- es$test
    if (laplacian) {
      train <- glaplacian(es$train, regularize = regularize)
    }

    ## graph spectral
    gs <- gspectral(train, k_max)
    U <- gs$u
    V <- gs$v


    ## sequential test statistics
    for (k in 1:k_max) {
      tick <- tick + 1
      gds <- gdstat(full = A, test = test, u = U[, k], v = V[, k], split = split)
      cv_stats[tick, 1:5] <- matrix(c(boot, k, gds), nrow = 1)
    }
  }

  ## summarize across CV/bootstrap
  if (bootstrap > 1) {
    cv_means <- cv_stats %>%
      group_by(.data$k) %>%
      summarise(
        cv_lambda_A = mean(.data$cv_lambda_A),
        cv_lambda_L = mean(.data$cv_lambda_L),
        z = mean(.data$z)
      ) %>%
      ungroup()
  } else {
    cv_means <- cv_stats
  }
  cv_means <- cv_means %>%
    mutate(
      pvals = pnorm(.data$z, lower.tail = FALSE),
      pvals = pmax(.data$pvals, ptol)
    ) ## avoid exact 0

  ## correct for multiplicity
  if (is.null(correct)) {
    correct <- "none"
  }
  cv_means <- mutate(cv_means, padj = p.adjust(.data$pvals, method = correct))

  ## inference
  criteria <- cv_means$padj
  k_stop <- which(criteria > alpha)
  k_infer <- ifelse(length(k_stop), min(k_stop) - 1, k_max)
  res <- list(
    inference = k_infer,
    summary = cv_means,
    bootstrap = bootstrap,
    split = split,
    alpha = alpha
  )

  if (trace) {
    res$stats <- cv_stats
  }
  class(res) <- "eigcv"
  return(res)
}


#' Print `eigcv`
#'
#' @method print eigcv
#'
#' @param x an `eigcv` object.
#' @param verbose `logical(1)`, whether to print out the summary table.
#' @param ... additional input to generic [print].
#' @return Print an `eigcv` object interactively.
#' @export
print.eigcv <- function(x,  ...) {
  cat("Estimated graph dimension:\t", x$inference, fill = TRUE)
    cat("\nNumber of bootstraps:\t\t", x$bootstrap, fill = TRUE)
    cat("Edge splitting probabaility:\t", x$split, fill = TRUE)
    cat("Significance level:\t\t", x$alpha, fill = TRUE)
    cat("\n ------------ Summary of Tests ------------\n")
    print(data.frame(x$summary[, -c(2, 3)]), row.names = FALSE)
    cat(fill = TRUE)
}

#' Plot `eigcv`
#'
#' @method plot eigcv
#'
#' @param x an `eigcv` object.
#' @param type either "z", "A", or "L" to specify the y-axis of the plot.
#'   If "z", plot the test statistics (asymptotic z score) for each k.
#'   If "A", plot x'Ax for each eigenvector x.
#'   If "L", plot x'Lx for each eigenvector x.
#' @param threshold `numeric(1)`, cut-off of p-value (in log10), default to 2.
#' @param ... ignored.
#' @return Plot an `eigcv` object.
#' @importFrom ggplot2 ggplot aes labs theme_bw theme scale_color_manual
#' @importFrom ggplot2 geom_hline geom_point geom_line scale_x_continuous
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
plot.eigcv <- function(x, type = c("z", "A", "L"), threshold = 2, ...) {
  stopifnot("x must contains an object called stats." = !is.null(x$summary))
  stopifnot("Threshold of statistics must be greater than 0." = threshold > 0)

  type <- type[1]
  if (type == "z") {
    dat <- x$summary %>% select(.data$k, val = .data$z)
    ylab <- "z score"
  }
  if (type == "A") {
    dat <- x$summary %>% select(.data$k, val = .data$cv_lambda_A)
    ylab <- "cross validated x' A x"
  }
  if (type == "L") {
    dat <- x$summary %>% select(.data$k, val = .data$cv_lambda_L)
    ylab <- "cross validated x' L x"
  }

  g <- ggplot(aes(.data$k, .data$val), data = dat) +
    geom_point(alpha = .8) +
    geom_line(color = "blue") +
    theme_bw() +
    ggplot2::scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
    })

  if (type == "z") {
    g <- g +
      geom_hline(
        yintercept = threshold, alpha = .8,
        linetype = 2, color = "grey60", show.legend = TRUE
      )
  }
  return(g)
}
