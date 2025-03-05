source("R/SVMFeature.R")
archivo <- "data/prueba.txt"
inputs <- c('x1', 'x2', 'x3')
output <- 'y'
costes <- c(5, 27, 10)
tam_pob <- 10
num_fea <- 2
n_iter <- 10  # Número de iteraciones
tiempo_max <- 300  # Tiempo máximo en segundos
modo <- "tiempo"  # Modo de ejecución: "iteraciones" o "tiempo"
svm_feature <- SVMFeature$new(archivo, inputs, output, costes, tam_pob, num_fea, n_iter, tiempo_max, modo)
svm_feature$run()
