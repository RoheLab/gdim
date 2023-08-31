test_that("multiplication works", {

  set.seed(17)

  M <- rsparsematrix(8, 12, nnz = 30, rand.x = NULL)

  graph_parts <- split_graph(M)

  expect_equal(
    graph_parts$train + graph_parts$test,
    methods::as(M, "dMatrix")
  )
})
