source("R/SVMFeature.R")

datos <- read.table("data/prueba2.txt", sep = "\t", header = TRUE)
print(datos)

# Definir variables iniciales
inputs <- paste0('x', 0:277)  # Columnas de entrada
output <- 'y1'  # Columna de salida
costs <- rep(1, length(inputs))  # Asignar un coste de 1 a cada característica
tam_pob <- 50
num_fea <- 5
n_iter <- 10  # Número de iteraciones
tiempo_max <- 300  # Tiempo máximo en segundos
modo <- "iters"  # Modo de ejecución: "iteraciones" o "tiempo"

# Crear el objeto SVMFeature
svm_feature <- SVMFeature$new(datos, inputs, output, costs, tam_pob, num_fea, n_iter, tiempo_max, modo)

# Ejecutar el algoritmo
svm_feature$run()

