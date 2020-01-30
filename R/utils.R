split_edges <- function(A, alpha = 0.05, symmetric = FALSE, check_symmetric = TRUE) {

  s <- as.data.frame(Matrix::summary(A))

  if (!symmetric) {

    index <- sample.int(nrow(s), size = round(nrow(s) * alpha))

    train <- sparseMatrix(s[-index, ]$i, s[-index, ]$j, x = s[-index, ]$x)
    test <- sparseMatrix(s[index, ]$i, s[index, ]$j, x = s[index, ]$x)

  } else {

    if (check_symmetric) {
      stopifnot(isSymmetric(A))
    }

    # only consider the lower triangle
    s <- dplyr::filter(s, i <= j)

    # NOTE: only count symmetric edges *once*
    index <- sample.int(nrow(s), size = round(nrow(s) * alpha))

    train <- sparseMatrix(
      s[-index, ]$i, s[-index, ]$j, x = s[-index, ]$x, symmetric = TRUE
    )

    test <- sparseMatrix(
      s[index, ]$i, s[index, ]$j, x = s[index, ]$x, symmetric = TRUE
    )
  }

  list(train = train, test = test)
}

# TODO: allow symmetric edge splitting
sample_edges <- function(A, alpha = 0.05) {

  s <- as.data.frame(summary(A))

  num_train <- round(nrow(s) * (1 - alpha))
  num_test <- nrow(s) - num_train

  itrain <- sample.int(nrow(s), size = num_train)
  itest <- sample.int(nrow(s), size = num_test)

  list(
    train = sparseMatrix(s[itrain, ]$i, s[itrain, ]$j, x = s[itrain, ]$x),
    test = sparseMatrix(s[itest, ]$i, s[itest, ]$j, x = s[itest, ]$x)
  )
}

permute_rows <- function(U) {
  stopifnot(ncol(U) > 1)  # avoid calling sample on scalars
  t(apply(U, 1, sample))
}
