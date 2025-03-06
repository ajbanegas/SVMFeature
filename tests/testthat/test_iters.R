library(testthat)
library(SVMFeature)

test_that("SVMFeature initializes correctly", {
  #data <- data.frame(x = 1:4, y = 5:8)
  #grid <- GRID(data, inputs = c("x"), outputs = c("y"), d = 2)

  #expect_s3_class(grid, "GRID")  # Verifica que sea de clase GRID
  #expect_equal(grid$data, data)  # Verifica que los datos coincidan
  #expect_equal(grid$d, 2)        # Verifica que d se asignó correctamente
  #expect_null(grid$data_grid)    # Verifica que data_grid es NULL inicialmente
  #expect_null(grid$knot_list)    # Verifica que knot_list es NULL inicialmente

  datos <- data.frame(x1 = 1:4, x2 = 11:14, y = 5:8)

  inputs <- c('x1', 'x2')
  output <- 'y'
  costs <- rep(1, length(inputs))
  pop_size <- 50
  num_fea <- 1
  n_iter <- 10
  max_time <- 300
  mode <- "iters"

  svm <- SVMFeature(datos, inputs, output, costs, pop_size, num_fea, n_iter, max_time, mode)

  expect_s3_class(svm, "SVMFeature")
  expect_equal(svm$mode, "iters")
  expect_equal(svm$num_fea, num_fea)
})
