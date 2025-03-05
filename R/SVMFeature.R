source("R/solucion.R")
source("R/poblacion.R")

library(dplyr)
library(magrittr)
library(ggplot2)

#' SVMFeature Class
#'
#' @description
#' This class implements a feature selection algorithm using Support Vector Machines (SVM).
#' It initializes with user-provided data, runs the optimization process, and prints the final population.
#'
#' @param archivo Character. The file path of the input data.
#' @param inputs Character vector. The names of the input features.
#' @param output Character. The name of the output variable.
#' @param costes Numeric vector. The costs associated with the features.
#' @param tam_pob Numeric. The population size.
#' @param num_fea Numeric. The number of features.
#' @param n_iter Numeric. The number of iterations.
#' @return An instance of the SVMFeature class.
#' @export
#'
#' @examples
#' archivo <- "data/prueba.txt"
#' inputs <- c('x1', 'x2', 'x3')
#' output <- 'y'
#' costes <- c(5, 27, 10)
#' tam_pob <- 10
#' num_fea <- 2
#  n_iter <- 10  # Número de iteraciones
#  tiempo_max <- 300  # Tiempo máximo en segundos
#  modo <- "tiempo"  # Modo de ejecución: "iteraciones" o "tiempo"
#' svm_feature <- SVMFeature$new(archivo, inputs, output, costes, tam_pob, num_fea, n_iter, tiempo_max, modo)
#' svm_feature$run()
SVMFeature <- setRefClass(
  "SVMFeature",
  fields = list(
    archivo = "character",
    inputs = "character",
    output = "character",
    costes = "numeric",
    tam_pob = "numeric",
    num_fea = "numeric",
    n_iter = "numeric",
    tiempo_max = "numeric",
    modo = "character",
    poblacion = "ANY",
    mejor_poblacion = "ANY"
  ),

  methods = list(
    #' Initialize the SVMFeature object
    #'
    #' @param archivo Character. The file path of the input data.
    #' @param inputs Character vector. The names of the input features.
    #' @param output Character. The name of the output variable.
    #' @param costes Numeric vector. The costs associated with the features.
    #' @param tam_pob Numeric. The population size.
    #' @param num_fea Numeric. The number of features.
    #' @param n_iter Numeric. The number of iterations.
    #' @return An initialized SVMFeature object.
    initialize = function(archivo, inputs, output, costes, tam_pob, num_fea, n_iter = 10, tiempo_max = 300, modo = "iteraciones") {
      .self$archivo <- archivo
      .self$inputs <- inputs
      .self$output <- output
      .self$costes <- costes
      .self$tam_pob <- tam_pob
      .self$num_fea <- num_fea
      .self$n_iter <- n_iter
      .self$tiempo_max <- tiempo_max
      .self$modo <- modo

      datos <- read.table(.self$archivo, sep = "\t", header = TRUE)

      # Imputar valores faltantes con la media de la columna
      for (col in names(datos)) {
        if (any(is.na(datos[[col]]))) {
          datos[[col]][is.na(datos[[col]])] <- mean(datos[[col]], na.rm = TRUE)
        }
      }

      # Normalización de los datos
      scaler <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (rng[1] == rng[2]) {
          return(rep(0.5, length(x)))  # Si todos los valores son iguales, asignar 0.5 a todos
        } else {
          return((x - rng[1]) / (rng[2] - rng[1]))
        }
      }

      datos_norm <- as.data.frame(lapply(datos[, inputs], scaler))
      datos_norm[output] <- datos[[output]]

      # Depuramos los datos cargados
      cat("Datos cargados y normalizados:\n")
      print(head(datos_norm))

      .self$poblacion <- Poblacion(data = datos_norm, costes = .self$costes, tam_pob = .self$tam_pob, inputs = .self$inputs, output = .self$output, num_features = .self$num_fea)
      .self$mejor_poblacion <- NULL
    },

    #' Run the SVMFeature optimization process
    #'
    #' @return The final population after the optimization process.
    run = function() {
      .self$poblacion <- generar_poblacion_inicial(.self$poblacion)
      .self$poblacion <- fnds(.self$poblacion)
      .self$actualizar_df_soluciones()

      .self$mejor_poblacion <- .self$poblacion

      tiempo_inicial <- Sys.time()

      ejecutar_iteracion <- function() {
        .self$poblacion <- nueva_poblacion(.self$poblacion)
        .self$poblacion <- fnds(.self$poblacion)

        poblacion_reducida <- reducir_poblacion(.self$poblacion)
        poblacion_reducida <- fnds(poblacion_reducida)

        .self$poblacion$lista_soluciones <- poblacion_reducida$lista_soluciones
        .self$actualizar_df_soluciones()
        rm(poblacion_reducida)

        # Comparar y mantener la mejor población basada en el número de soluciones en el frente 1
        soluciones_actuales <- nrow(filter(.self$poblacion$df_soluciones, FRONT == 1))
        soluciones_mejor <- nrow(filter(.self$mejor_poblacion$df_soluciones, FRONT == 1))

        if (soluciones_actuales > soluciones_mejor) {
          .self$mejor_poblacion <- .self$poblacion
        }
      }

      if (.self$modo == "iteraciones") {
        iter <- 0
        while (iter < .self$n_iter) {
          cat(sprintf("ITERACION %d: %f\n", iter, as.numeric(difftime(Sys.time(), tiempo_inicial, units = "secs"))))
          ejecutar_iteracion()
          iter <- iter + 1
        }
      } else if (.self$modo == "tiempo") {
        while (as.numeric(difftime(Sys.time(), tiempo_inicial, units = "secs")) < .self$tiempo_max) {
          cat(sprintf("TIEMPO TRANSCURRIDO: %f\n", as.numeric(difftime(Sys.time(), tiempo_inicial, units = "secs"))))
          ejecutar_iteracion()
        }
      } else {
        stop("Modo de ejecución no válido. Use 'iteraciones' o 'tiempo'.")
      }

      imprimir_poblacion(.self$mejor_poblacion)

      if (is.null(.self$mejor_poblacion$df_soluciones)) {
        stop("Error: df_soluciones es NULL.")
      }

      poblacion_front1 <- filter(.self$mejor_poblacion$df_soluciones, FRONT == 1)

      if (nrow(poblacion_front1) == 0) {
        cat("No hay soluciones en el frente 1 para graficar.\n")
      } else {
        poblacion_front1$DIST <- as.numeric(poblacion_front1$DIST)
        poblacion_front1$EPS <- as.numeric(poblacion_front1$EPS)

        plot <- ggplot(poblacion_front1, aes(x = DIST, y = EPS)) +
          geom_point(color = 'blue') +
          labs(x = 'Distancia', y = 'Epsilon', title = 'Gráfico de Soluciones en el Frente 1') +
          theme_minimal() +
          theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray"))

        print(plot)
      }
    },

    actualizar_df_soluciones = function() {
      df_soluciones <- do.call(rbind, lapply(.self$poblacion$lista_soluciones, function(solucion) {
        solucion_dict <- to_dict(solucion)
        solucion_dict <- lapply(solucion_dict, function(item) if (is.list(item)) unlist(item) else item)
        solucion_df <- as.data.frame(t(solucion_dict), stringsAsFactors = FALSE)
        return(solucion_df)
      }))
      df_soluciones$FRONT <- as.numeric(as.character(df_soluciones$FRONT))
      df_soluciones <- df_soluciones[order(df_soluciones$FRONT), ]
      .self$poblacion$df_soluciones <- df_soluciones

      # Imprimir número de soluciones por frente
      fronts_count <- df_soluciones %>% group_by(FRONT) %>% summarize(count = n())
      cat("Número de soluciones por frente:\n")
      print(fronts_count)
    }
  )
)
