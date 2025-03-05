source("R/Solution.R")
source("R/Population.R")

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
#  modo <- "time"  # Modo de ejecución: "iters" o "time"
#' svm_feature <- SVMFeature$new(archivo, inputs, output, costes, tam_pob, num_fea, n_iter, tiempo_max, modo)
#' svm_feature$run()
SVMFeature <- setRefClass(
  "SVMFeature",
  fields = list(
    datos = "ANY",
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
    #' @param datos DataFrame. The input data.
    #' @param inputs Character vector. The names of the input features.
    #' @param output Character. The name of the output variable.
    #' @param costes Numeric vector. The costs associated with the features.
    #' @param tam_pob Numeric. The population size.
    #' @param num_fea Numeric. The number of features.
    #' @param n_iter Numeric. The number of iterations.
    #' @return An initialized SVMFeature object.
    initialize = function(datos, inputs, output, costes, tam_pob, num_fea, n_iter = 10, tiempo_max = 300, modo = "iters") {
      .self$datos <- datos
      .self$inputs <- inputs
      .self$output <- output
      .self$costes <- costes
      .self$tam_pob <- tam_pob
      .self$num_fea <- num_fea
      .self$n_iter <- n_iter
      .self$tiempo_max <- tiempo_max
      .self$modo <- modo

      # Normalización de los datos
      scaler <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (rng[1] == rng[2]) {
          return(rep(0.5, length(x)))  # Si todos los valores son iguales, asignar 0.5 a todos
        } else {
          return((x - rng[1]) / (rng[2] - rng[1]))
        }
      }

      datos_norm <- as.data.frame(lapply(.self$datos[, inputs], scaler))
      datos_norm[output] <- .self$datos[[output]]

      # Depuramos los datos cargados
      cat("Normalized data:\n")
      print(head(datos_norm))

      .self$poblacion <- Population(data = datos_norm, costes = .self$costes, tam_pob = .self$tam_pob, inputs = .self$inputs, output = .self$output, num_features = .self$num_fea)
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

      if (.self$modo == "iters") {
        iter <- 0
        while (iter < .self$n_iter) {
          cat(sprintf("ITERATION %d: %f\n", iter, as.numeric(difftime(Sys.time(), tiempo_inicial, units = "secs"))))
          ejecutar_iteracion()
          iter <- iter + 1
        }
      } else if (.self$modo == "time") {
        while (as.numeric(difftime(Sys.time(), tiempo_inicial, units = "secs")) < .self$tiempo_max) {
          cat(sprintf("TIME: %f\n", as.numeric(difftime(Sys.time(), tiempo_inicial, units = "secs"))))
          ejecutar_iteracion()
        }
      } else {
        stop("Invalid execution mode. Usage: 'iters' or 'time'.")
      }

      imprimir_poblacion(.self$mejor_poblacion)

      if (is.null(.self$mejor_poblacion$df_soluciones)) {
        stop("Error: df_soluciones is NULL.")
      }

      poblacion_front1 <- filter(.self$mejor_poblacion$df_soluciones, FRONT == 1)

      if (nrow(poblacion_front1) == 0) {
        cat("There is no solution in front 1 to plot.\n")
      } else {
        poblacion_front1$DIST <- as.numeric(poblacion_front1$DIST)
        poblacion_front1$EPS <- as.numeric(poblacion_front1$EPS)

        plot <- ggplot(poblacion_front1, aes(x = DIST, y = EPS)) +
          geom_point(color = 'blue') +
          labs(x = 'Distance', y = 'Epsilon', title = 'Solutions in Front 1') +
          theme_minimal() +
          theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray"))

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
      cat("Number of solutions by front:\n")
      print(fronts_count)
    }
  )
)
