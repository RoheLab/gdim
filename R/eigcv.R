split_graph <- function(A, test_portion = 0.1) {

  A <- methods::as(A, "dMatrix")
  A <- methods::as(A, "TsparseMatrix")

  # assumes elements of A are Poisson (i.e. non-negative integers)
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

# Given the trained left/right singular vectors, compute the test statistic for
# graph dimension.
gdstat <- function(full, test, u, v, test_portion) {
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



#' Edge Bootstrapping and Splitting
#'
#' Estimate the graph dimension via eigenvalue cross-validation (EigCV).
#' A graph has dimension `k` if the first `k` eigenvectors of its adjacency
#' matrix are correlated with its population eigenspace, and the others are not.
#' Edge bootstrapping sub-samples the edges of the graph (without replacement).
#' Edge splitting separates the edges into a training part and a testing part.
#'
#' @param A The adjacency matrix of graph. Must be non-negative integer valued.
#' @param k_max `integer(1)`, number of eigenvectors to compute.
#' @param num_bootstraps `integer(1)`, number of graph bootstraps, default to 10.
#'   Graph bootstrapping is to account for the randomness in graph splitting,
#'   rather than obtaining any statistic (as a traditional num_bootstraps does).
#'   Hence, a small number (e.g., 3~10) of bootstraps usually suffices.
#'   If `num_bootstraps>1`, the test statistics will be averaged across bootstraps
#'   and the p-values will be calculated based on the averaged statistics.
#' @param alpha Significance level of each test, defaults to `0.05`.
#'   This is used to cut off the dimension estimation.
#' @param ptol `numeric(1)`, the tolerance of minimal p-value.
#' @param regularize TODO
#' @param test_portion TODO
#' @inheritParams stats::p.adjust
#' @param laplacian `logical(1)`, use the normalized and regularized adjacency
#'   matrix (i.e. L)
#'   This option is experimental and should be used with caution.
#' @return A `eigcv` object, which contains:
#'   \item{estimated_dimension}{inferred graph dimension.}
#'   \item{summary}{summary table of the tests.}
#'   \item{num_bootstraps}{number of bootstraps performed.}
#'   \item{test_portion}{graph splitting probability used.}
#'   \item{alpha}{significance level of each test.}
#'
#' @importFrom dplyr summarize group_by ungroup mutate summarise
#' @importFrom magrittr %>%
#' @importFrom rlang .data
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
#' eigcv_result <- eigcv(A, k_max = 10)
#' eigcv_result
#'
eigcv <- function(A, k_max,
                  ...,
                  num_bootstraps = 10, test_portion = 0.1,
                  alpha = 0.05,
                  ptol = .Machine$double.eps,
                  method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
                             "none"),
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

  ## summarize across CV/num_bootstraps
  if (num_bootstraps > 1) {
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
      pvals = stats::pnorm(.data$z, lower.tail = FALSE),
      pvals = pmax(.data$pvals, ptol)
    ) ## avoid exact 0

  ## correct for multiplicity
  cv_means$padj <- stats::p.adjust(cv_means$pvals, method = method)

  ## inference
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


#' Print `eigcv`
#'
#' @method print eigcv
#'
#' @param x an `eigcv` object.
#' @param ... Ignored.
#' @export
print.eigcv <- function(x,  ...) {
  cat("Estimated graph dimension:\t", x$estimated_dimension, fill = TRUE)
    cat("\nNumber of bootstraps:\t\t", x$num_bootstraps, fill = TRUE)
    cat("Edge splitting probabaility:\t", x$test_portion, fill = TRUE)
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
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
    })

  if (type == "z") {
    g <- g +
      geom_hline(
        yintercept = threshold, alpha = .8,
        linetype = 2, color = "grey60", show.legend = TRUE
      )
  }

  g
}
