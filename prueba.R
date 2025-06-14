source("R/SVMFeature.R")

datos <- read.table("data/prueba.txt", sep = ";", header = TRUE)

inputs <- c('x1', 'x2', 'x3')
output <- 'y'
costs <- c(5, 27, 10)
tam_pob <- 10
num_fea <- 2
n_iter <- 10  # Número de iteraciones
tiempo_max <- 10  # Tiempo máximo en segundos (300)
modo <- "time"  # Modo de ejecución: "iters" o "time"
objective <- "confusion-matrix" # "distance-epsilon"

svm_feature <- SVMFeature(datos, inputs, output, costs, tam_pob, num_fea,
                          n_iter, tiempo_max, modo, objective)
svm_feature <- run.SVMFeature(svm_feature)
