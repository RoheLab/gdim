#' From cluster label to cluster membership matrix
#'
#' @param z a `integer` vector, cluster label 1,2,...
#' @return a binary `sparseMatrix` of cluster membership.
#' @importFrom Matrix sparseMatrix
z2Z <- function(z) {
  f <- as.factor(z)
  sparseMatrix(
    i = seq_along(f), j = as.numeric(f), x = 1,
    dimnames = list(names(f), levels(f))
  )
}



#' Spectral Clustering On Ratios-of-Eigenvectors
#'
#' Perform Spectral Clustering On Ratios-of-Eigenvectors (SCORE).
#'
#' @inheritParams StGoF
#' @param K `integer(1)`, number of dimension.
#' @param truncate `numeric(1)`, truncate the eigen-ratio (in absolute value)
#' matrix, default to `log(nrow(A))`. Set this to `0` if don't want truncation.
#' @inheritParams stats::kmeans
#' @return a vector of `1:K`, the community labels of all nodes.
#' @importFrom RSpectra eigs svds
#' @importFrom stats kmeans
#' @references Jin, Jiashun. Fast community detection by SCORE. Ann. Statist. 43 (2015), no. 1, 57--89. doi:10.1214/14-AOS1265. \url{https://projecteuclid.org/euclid.aos/1416322036}
SCORE <- function(A, K,
                  truncate = log(nrow(A)),
                  iter.max = 1000,
                  nstart = 1000,
                  algorithm = "Hartigan-Wong",
                  trace = FALSE) {
  # Initial checks
  stopifnot("K should be larger than 1.\n" = K > 1)

  ## SCORE
  if (!Matrix::isSymmetric(A)) {
    warning(
      "The adjacency matrix is not symmetric. ",
      "Use the singular vector version of the SCORE algorithm."
    )
    s <- svds(A, K)
    vb <- s$v
  } else {
    s <- eigs(A, K)
    vb <- s$vectors
  }
  vb2 <- vb[, 2:K, drop = FALSE]
  rb <- vb2 / vb[, 1]
  rb[is.nan(rb)] <- 0
  if (truncate) {
    rb[rb > truncate] <- truncate
    rb[rb < -truncate] <- -truncate
  }
  label <- kmeans(rb, K,
    iter.max = iter.max, nstart = nstart,
    algorithm = algorithm, trace = trace
  )$cluster
}


#' Refitted Quadrilateral
#'
#' Calculate the Refitted Quadrilateral (RQ) statistics.
#'
#' @inheritParams StGoF
#' @param label a vector of `1:K`, the community labels of all nodes.
#' @return The RQ statistics.
# @importFrom Matrix Diagonal
#' @references Jin, Jiashun. Fast community detection by SCORE. Ann. Statist. 43 (2015), no. 1, 57--89. doi:10.1214/14-AOS1265. \url{https://projecteuclid.org/euclid.aos/1416322036}
RQ <- function(A, label) {
  stopifnot(
    "The dimensions of `A` and `label` disagree." =
      nrow(A) == length(label)
  )
  ## should also check label input here

  n <- nrow(A)
  Pi_hat <- z2Z(label)
  k <- ncol(Pi_hat)
  stopifnot("label must have at least 2 levels" = k > 1)

  APi <- A %*% Pi_hat
  PiAPi <- t(Pi_hat) %*% APi
  P_diag_inv <- Diagonal(k, 1 / diag(PiAPi))
  d <- rowSums(A)
  IkAIn <- as.numeric(t(Pi_hat) %*% d)
  P_hat_tilde <- Diagonal(k, 1 / IkAIn) %*% PiAPi %*% Diagonal(k, 1 / IkAIn)

  Pi_Ptilde_Pi <- Pi_hat %*% P_hat_tilde %*% t(Pi_hat)
  A2 <- Pi_Ptilde_Pi * d
  Omega_hat <- t(A2) * d
  A_tilde <- A - Omega_hat

  S <- A_tilde - Diagonal(n, diag(A_tilde))
  S2 <- S %*% S
  S3 <- S2 %*% S
  s <- diag(S)

  Q <- as.numeric(sum(sum(S * S3)) -
    4 * sum(s * diag(S3)) +
    8 * sum((s^2) * diag(S2)) -
    6 * sum(s^4) -
    2 * sum(diag(S2)^2) +
    2 * t(s) %*% S^2 %*% s +
    sum(S^4))

  P_diag <- diag(PiAPi)
  P_delta <- sqrt(P_diag) / IkAIn
  theta_hat <- rowSums(A) * rowSums(Pi_hat %*% Diagonal(k, P_delta))
  Theta_hat <- Diagonal(n, theta_hat)
  APi <- A %*% Pi_hat
  PiAPi <- t(Pi_hat) %*% APi
  P_diag_inv_sqrt <- Diagonal(k, 1 / sqrt(pmax(1, diag(PiAPi))))
  P_hat <- P_diag_inv_sqrt %*% PiAPi %*% P_diag_inv_sqrt
  diag(P_hat) <- diag(P_hat) + (diag(P_hat) == 0)

  G <- t(Pi_hat) %*% Theta_hat %*% Pi_hat / sum(diag(Theta_hat))
  H2 <- t(Pi_hat) %*% Theta_hat^2 %*% Pi_hat
  PH2P <- P_hat %*% H2 %*% P_hat
  g <- diag(G)
  gtilde <- g / (P_hat %*% g)

  ## mean
  # theta_norm_2 = sum(diag(Theta_hat).^2)^0.5
  bn <- as.numeric(2 * t(gtilde) %*% PH2P^2 %*% gtilde)

  ## variance
  B <- A %*% A
  C4 <- sum(diag(B %*% B)) - 2 * sum(B) + sum(A)

  ## statistics
  SgnQstat <- (Q - bn) / sqrt(8 * C4)
  return(as.numeric(SgnQstat))
}

#' Stepwise Goodness-of-Fit
#'
#' Perform the stepwise goodness-of-fit (StGoF) tests of graph dimension.
#' The tests are run from k=2, 3, ..., K_max separately.
#' Given a significance level (`alpha`), the infered graph dimension is
#' determined by the first test (i.e., smallest k) that is accepted.
#'
#' @param A a `matrix` or `Matrix`, the adjacency matrix of an _undirected_ network.
#' @param K_max `integer(1)`, number of communities.
#' @param alpha `numeric(1)`, significance level, default to 0.05.
#'   Set `alpha=1` if you want to compute all the statistics from 1 to K_max.
#' @return a `data.table` with three columns:
#' \item{k}{tested graph dimension.}
#' \item{rq.stat}{the refitted quadrilateral test statistics.}
#' \item{p.value}{test p-value (two-sided, no adjustment).}
#' @references Jin, Jiashun. Fast community detection by SCORE. Ann. Statist. 43 (2015), no. 1, 57--89. doi:10.1214/14-AOS1265. \url{https://projecteuclid.org/euclid.aos/1416322036}
#' @importFrom dplyr bind_rows
#' @export
StGoF <- function(A, K_max, alpha = 0.05) {
  stopifnot(K_max > 1)
  tbl <- lapply(2:K_max, function(k) {
    label <- SCORE(A, k)
    z <- RQ(A, label)
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    data.frame(k = k, rq.stat = z, p.value = p)
  })
  tbl <- bind_rows(tbl)

  k_stop <- which(tbl$p.value > alpha)
  k_infer <- ifelse(length(k_stop), min(k_stop) + 1, K_max) ## k starts from 2
  res <- list(
    inference = k_infer,
    summary = tbl,
    alpha = alpha
  )
  class(res) <- c("stgof", class(res))
  return(res)
}

#' Print `stgof`
#'
#' @method print stgof
#'
#' @param x an `stgof` object.
#' @param verbose `logical(1)`, whether to print out the summary table.
#' @param ... additional input to generic [print].
#' @return Print an `stgof` object interactively.
#' @export
print.stgof <- function(x, verbose = TRUE, ...) {
  cat("Estimated graph dimension:\t", x$inference, fill = TRUE)
  if (verbose) {
    cat("Significance level:\t\t", x$alpha, fill = TRUE)
    cat("\n ------------ Summary of Tests ------------\n")
    print(x$summary, row.names = FALSE)
    cat(fill = TRUE)
  }
}
