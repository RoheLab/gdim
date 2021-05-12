## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## utility for graph dimensionality estimation 
## 9/25/2020
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  # library(irlba) ## this loads `Matrix`
  library(RSpectra)
  # library(igraph)
  library(reshape2)
  library(metap)
})


## ---- graph dim test ---- 

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
    el$x = 1
  }
  isSymmetric <- Matrix::isSymmetric(A)
  if (isSymmetric) {
    ## simplify edge list; this implies A is a square matrix 
    el <- el[el$i <= el$j, ]
  }
  # stopifnot("A contains negative/fractional" = all(el$x %% 1 == 0 & el$x >= 0))
  
  ## splitting
  leave <- (runif(nrow(el)) < split) * 1
  
  train <- sparseMatrix(
    i = el$i,
    j = el$j,
    x = el$x * (1 - leave),
    dims = dim(A), 
    dimnames = dimnames(A), 
    symmetric = isSymmetric, 
    triangular = isTriangular(A)
  )
  
  test <- sparseMatrix(
    i = el$i,
    j = el$j,
    x = el$x * leave,
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
      stop("Cannot use Laplacian because some nodes are isolated. ", 
           "Set either \"regularize=TRUE\" or \"laplacian=FALSE\" option.")
    }
    tau_row <- tau_col <- 0
  }
  D_row = Diagonal(nrow(A), 1 / sqrt(deg_row + tau_row))
  D_col = Diagonal(ncol(A), 1 / sqrt(deg_col + tau_col))
  L = D_row %*% A %*% D_col
  return(L)
}

#' Graph Spectral Decomposition
#' 
#' Given ab `m` by `n` graph matrix `A`, find the largest `k` singular values 
#' and the corresponding singular vectors. 
#' 
#' @inheritParams esplit 
#' @param k `integer(1)`, number of eigenvalues.
#' @param ... additonal parameters to pass into [RSpectra::svds].
#' @importFrom RSpectra eigs svds
#' @return a list (see also [RSpectra::svds]): 
#'   \item{u}{an `m` by `k` matrix whose columns contain the left singular vectors.}
#'   \item{v}{an `n` by `k` matrix whose columns contain the right singular vectors.}
gspectral <- function(A, k, ...) {
  stopifnot("k too large" = k <= min(dim(A))) 
  S <- svds(A, k, ...)
  return(list(u = S$u, v = S$v))
}

#' Graph Dimensionality Statistic
#' 
#' Given the trained left/right singular vectors, compute the test statistic for
#' graph dimensionality. 
#' 
#' @param full,test `matrix` or `Matrix`, the adjacency or Laplacian matrix of 
#'   the full and test graphs. 
#' @param u,v `numeric` vector, the trained left and right singular vectors. 
#' @inheritParams esplit
#' @return `numeric(1)`, test statistic
gdstat <- function(full, test, u, v, split) {
  se <- sqrt(as.numeric(split * t(u^2) %*% full %*% v^2)) ## standard error
  # se <- sqrt(as.numeric(t(u^2) %*% test %*% v^2)) ## standard error
  as.numeric(t(u) %*% test %*% v) / se ## test stat
}

#' Edge Bootstrapping and Splitting
#' 
#' Estimate the graph dimensionality via eigenvalue cross-validation (EigCV). 
#' A graph has dimensionality `k` if the first `k` eigenvectors of its adjacency 
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
#' @param meta either "z" or "p", how to aggregate bootstrap results. 
#'   If "z", will take the average of test statistics. 
#'   If "p", will perform meta analysis using [metap::sump]. 
#' @param ptol `numeric(1)`, the tolerance of minimal p-value. 
#'   This is only useful when `meta=="p"` (e.g., "bonferroni", "fdr"). 
#' @param correct `character(1)`, correction method for multiple testing. 
#'   This is recommended if `k_max` is large. See [p.adjust] for the acceptable 
#'   values of `method`. If `correct` is not "none", the returned summary table
#'   will have an extra column called `padj`.
#' @param laplacian `logical(1)`, use the (symmetric) graph Laplacian. 
#' @param align `logical(1)`, align the trained singular vectors with those of
#'   the full graph. The default is `FALSE`.
#'   This option is experimental and does not guarantee any improvement. 
#' @param trace `logical(1)`, for diagnostic use. If `TRUE`, return additionally 
#'   `stats`, the test statistics of every bootstrap, and `pvals`, the eigen 
#'   vectors of the (last, if `bootstrap>1`) training graph.
#' @return A `eigcv` object, which contains:
#'   \item{inference}{inferred graph dimensionality.} 
#'   \item{summary}{summary table of the tests.}
#'   \item{bootstrap}{number of bootstraps performed.}
#'   \item{split}{graph splitting probability used.} 
#'   \item{alpha}{significance level of each test.}
#'   \item{trace}{only if `trace=TRUE`, the spectral decomposition of the traning graph.}
#' @importFrom dplyr bind_rows
#' @importFrom metap logitp
#' @export
eigcv <- function(A, k_max, 
                  bootstrap = 10, split = 0.1, 
                  alpha = 0.05, meta = "z", 
                  ptol = 1e-32, correct = "none",
                  laplacian = TRUE, 
                  regularize = TRUE, 
                  align = FALSE, 
                  trace = FALSE) {
  ## check
  # stopifnot("A must be square"= nrow(A) == ncol(A))
  # stopifnot("A should be symmetric" = Matrix::isSymmetric(A))
  n = min(dim(A))
  stopifnot("\"k_max\" is too large." = k_max <= n)
  stopifnot("\"split\" must range between 0 and 1." = (split > 0 && split < 1))
  stopifnot("\"bootstrap\" must be a positive integer." = bootstrap >= 1)
  stopifnot("Unknown \"meta\" method." = meta %in% c("z", "p"))
  
  ## full graph
  full <- A <- A * 1
  if (laplacian) {
    full <- glaplacian(A, regularize = regularize)
  } 
  if (align) {
    gs_full <- gspectral(full, k_max)
  }
  
  ## bootstrap 
  stats <- sapply(seq_len(bootstrap), function(trial) {
    ## edge splitting 
    es <- esplit(A, split)
    train <- es$train
    test <- es$test
    if (laplacian) {
      train <- glaplacian(es$train, regularize = regularize)
      test <- glaplacian(es$test, regularize = regularize)
    } 
    
    ## graph spectral  
    gs <- gspectral(train, k_max)
    if (align) {
      ## apply Procrustes rotation
      Su = svd(crossprod(gs$u, gs_full$u))
      U <- gs$u %*% tcrossprod(Su$u, Su$v)
      Sv = svd(crossprod(gs$v, gs_full$v))
      V <- gs$v %*% tcrossprod(Sv$u, Sv$v)
    } else {
      U <- gs$u
      V <- gs$v
    }
    
    ## sequential test statistics
    mapply(gdstat, data.frame(U), data.frame(V), USE.NAMES = FALSE,
           MoreArgs = list(full = full, test = test, split = split))
  })
  dimnames(stats) <- list(k = seq_len(k_max), bootstrap = seq_len(bootstrap))
  pvals <- pmax(pnorm(stats, lower.tail = FALSE), ptol) ## one-sided, avoid exact 0 
  
  ## meta bootstraps
  zbar = rowMeans(stats)
  zmin = apply(stats, 1, min)
  if (meta == "z") {
    pval <- pnorm(zbar, lower.tail = FALSE) 
  } else if (meta == "p") {
    pval <- apply(pvals, 1, function(p) metap::sump(p)$p) ## one-sided
  }
  tbl = data.frame(k = seq_len(k_max), 
                   zbar = zbar, 
                   zmin = zmin,
                   pval = pval)
  if (!is.null(correct)) 
    tbl$padj <- p.adjust(p = tbl$pval, method = correct)
  
  ## inference
  criteria <- tbl[,ncol(tbl)]
  k_stop <- which(criteria > alpha)
  k_infer <- ifelse(length(k_stop), min(k_stop) - 1, k_max)
  res <- list(inference = k_infer, 
              summary = tbl, 
              bootstrap = bootstrap,
              split = split, 
              alpha = alpha)
  
  if (trace) {
    res$stats <- stats
    res$pvals <- pvals
    if (align) {
      res$spectral <- gs_full
    }
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
print.eigcv <- function(x, verbose = TRUE, ...) {
  cat("Estimated graph dimensionality:\t", x$inference, fill = TRUE)
  if (verbose) {
    cat("\nNumber of bootstraps:\t\t", x$bootstrap, fill = TRUE)
    cat("Edge splitting probabaility:\t", x$split, fill = TRUE)
    cat("Significance level:\t\t", x$alpha, fill = TRUE)
    cat("\n ------------ Summary of Tests ------------\n")
    print(x$summary, row.names = FALSE)
    cat(fill = TRUE)
  }
}

#' Plot `eigcv`
#' 
#' @method plot eigcv
#' 
#' @param x an `eigcv` object.
#' @param type either "z" or "p" to specify the y-axis of the plot. 
#'   If "z", plot the test statistics (asymptotic z score) for each k. 
#'   If "p", plot the log10 transformed p-values.
#' @param threshold `numeric(1)`, cut-off of p-value (in log10), default to 2. 
#' @param verbose `logical(1)`, whether to print out the summary table. 
#' @param ... ignored.
#' @return Plot an `eigcv` object.
#' @importFrom ggplot2 ggplot aes labs theme_bw theme scale_color_manual
#' @importFrom ggplot2 geom_hline geom_smooth geom_point geom_line geom_pointrange
#' @importFrom reshape2 dcast melt
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate group_by summarise
#' @importFrom rlang .data
#' @export
plot.eigcv <- function(x, type = "z", threshold = 2, 
                       verbose = FALSE, ...) {
  if (type == "z") {
    stopifnot("x must contains an object called stats." = !is.null(x$stats))
    stopifnot("Threshold of statistics must be greater than 0." = threshold > 0)
    ylab <- "z score"
    yintercept = threshold
    tbl <- x$stats %>% 
      melt(value.name = "stat") %>%
      filter(.data$k > 1) %>% 
      group_by(.data$k) %>%
      summarise(y = mean(.data$stat), 
                ymin = min(.data$stat), 
                ymax = max(.data$stat)) %>%
      mutate(signif = ifelse(y < min(yintercept), "no", "yes")) 
    g <- ggplot(aes(.data$k, .data$y), data = tbl) + 
      geom_hline(yintercept = yintercept, alpha = .8, 
                 linetype = 2, color = "grey60", show.legend = TRUE) + 
      geom_pointrange(aes(ymin = .data$ymin, 
                          ymax = .data$ymax, 
                          color = .data$signif), 
                      alpha = .8) + 
      scale_color_manual(values = c("yes" = "grey30", 
                                    "no" = "#ef3b2c")) + 
      labs(x = "graph dimensionality", y = ylab, color = "signif.") + 
      geom_smooth() + 
      # scale_y_continuous(trans = "log2") + 
      theme_bw() 
  } 
  
  if (type == "p") {
    stopifnot("Threshold of p-value must be between 0 and 1 (preferably > 1e-20)." = 
                threshold > 1e-20 && threshold < 1)
    yintercept = -log10(threshold)
    ylab <- "-log10(p-value)"
    y <- pmin(-log10(x$summary[,ncol(x$summary)]), 20)
    tbl <- data.frame(k = x$summary$k, y = y, 
                      signif = ifelse(y < min(yintercept), "no", "yes")) 
    g <- tbl %>% 
      filter(.data$k > 1) %>% 
      ggplot(aes(.data$k, y)) + 
      geom_hline(yintercept = yintercept, alpha = .8, 
                 linetype = 2, color = "grey60", show.legend = TRUE) + 
      geom_line(alpha = .8) +
      geom_point(aes(color = .data$signif), alpha = .8) + 
      scale_color_manual(values = c("yes" = "grey30", "no" = "#ef3b2c")) + 
      labs(x = "graph dimensionality", y = ylab, color = "signif.") + 
      theme_bw() 
  }
  
  return(g)
}

