library(fastRG)

test_that("symmetric eigcv, calculations on A", {

  set.seed(27)

  B <- matrix(0.1, 5, 5)
  diag(B) <- 0.3

  model <- sbm(
    n = 1000,
    k = 5,
    B = B,
    expected_degree = 40,
    poisson_edges = FALSE,
    allow_self_loops = FALSE
  )

  A <- sample_sparse(model)

  eigcv_result <- eigcv(A, k_max = 10)

  # expect_silent(
  #   eigcv_result <- eigcv(A, k_max = 10)
  # )

  expect_equal(
    eigcv_result$estimated_dimension,
    5
  )
})

test_that("asymmetric eigcv, calculations on A", {

  set.seed(27)

  B <- matrix(0.1, nrow = 5, ncol = 8)
  diag(B) <- 0.9

  n <- 1000

  model <- directed_dcsbm(
    theta_in = rep(1, n),
    theta_out = rep(1, n),
    B = B,
    k_in = 5,
    k_out = 8,
    expected_density = 0.05
  )

  A <- sample_sparse(model)

  eigcv_result <- eigcv(A, k_max = 10)

  # expect_silent(
  #   eigcv_result <- eigcv(A, k_max = 10)
  # )

  expect_equal(
    eigcv_result$estimated_dimension,
    5
  )
})

library(fastRG)

test_that("symmetric eigcv, calculations on L", {

  set.seed(27)

  B <- matrix(0.1, 5, 5)
  diag(B) <- 0.3

  model <- sbm(
    n = 1000,
    k = 5,
    B = B,
    expected_degree = 40,
    poisson_edges = FALSE,
    allow_self_loops = FALSE
  )

  A <- sample_sparse(model)

  eigcv_result <- eigcv(A, k_max = 10, laplacian = TRUE)

  # expect_silent(
  #   eigcv_result <- eigcv(A, k_max = 10)
  # )

  expect_equal(
    eigcv_result$estimated_dimension,
    5
  )
})

test_that("asymmetric eigcv, calculations on A", {

  set.seed(27)

  B <- matrix(0.1, nrow = 5, ncol = 8)
  diag(B) <- 0.9

  n <- 1000

  model <- directed_dcsbm(
    theta_in = rep(1, n),
    theta_out = rep(1, n),
    B = B,
    k_in = 5,
    k_out = 8,
    expected_density = 0.05
  )

  A <- sample_sparse(model)

  eigcv_result <- eigcv(A, k_max = 10, laplacian = TRUE)

  # expect_silent(
  #   eigcv_result <- eigcv(A, k_max = 10)
  # )

  expect_equal(
    eigcv_result$estimated_dimension,
    5
  )
})

