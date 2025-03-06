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

  datos <- data.frame(x1 = 1:4, x2 = 10:14, y = 5:8)

  # Definir variables iniciales
  inputs <- c('x1', 'x2')  # Columnas de entrada
  output <- 'y'  # Columna de salida
  costs <- rep(1, length(inputs))  # Asignar un coste de 1 a cada característica
  tam_pob <- 50
  num_fea <- 5
  n_iter <- 10  # Número de iteraciones
  tiempo_max <- 300  # Tiempo máximo en segundos
  modo <- "iters"  # Modo de ejecución: "iteraciones" o "tiempo"

  # Crear el objeto SVMFeature
  svm_feature <- SVMFeature$new(datos, inputs, output, costs, tam_pob, num_fea, n_iter, tiempo_max, modo)

  expect_s3_class(grid, "SVMFeature")  # Verifica que sea de clase SVMFeature

  # Ejecutar el algoritmo
  #svm_feature$run()
})
