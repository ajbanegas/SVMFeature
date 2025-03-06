#source("R/Solution.R")
#source("R/Population.R")

#usethis::use_package("dplyr")
#usethis::use_package("ggplot2")
#usethis::use_package("magrittr")
library("dplyr")
library("ggplot2")
library("magrittr")


#' SVMFeature Class
#'
#' @description
#' This class implements a feature selection algorithm using Support Vector Machines (SVM).
#' It initializes with user-provided data, runs the optimization process, and prints the final population.
#'
#' @param data DataFrame. The input data.
#' @param inputs Character vector. The names of the input features.
#' @param output Character. The name of the output variable.
#' @param costs Numeric vector. The costs associated with the features.
#' @param pop_size Numeric. The population size.
#' @param num_fea Numeric. The number of features.
#' @param n_iter Numeric. The number of iterations.
#' @return An instance of the SVMFeature class.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize n
#'
#' @examples
#' data <- read.table("data/prueba.txt", sep = ";", header = TRUE)
#' inputs <- c('x1', 'x2', 'x3')
#' output <- 'y'
#' costs <- c(5, 27, 10)
#' pop_size <- 10
#' num_fea <- 2
#  n_iter <- 10  # Número de iteraciones
#  max_time <- 300  # Tiempo máximo en segundos
#  mode <- "time"  # Modo de ejecución: "iters" o "time"
#' svm_feature <- SVMFeature$new(archivo, inputs, output, costs, pop_size, num_fea, n_iter, max_time, mode)
#' svm_feature$run()
SVMFeature <- setRefClass(
  "SVMFeature",
  fields = list(
    data = "ANY",
    inputs = "character",
    output = "character",
    costs = "numeric",
    pop_size = "numeric",
    num_fea = "numeric",
    n_iter = "numeric",
    max_time = "numeric",
    mode = "character",
    population = "ANY",
    best_population = "ANY"
  ),
  methods = list(
    #' Initialize the SVMFeature object
    #'
    #' @param data DataFrame. The input data.
    #' @param inputs Character vector. The names of the input features.
    #' @param output Character. The name of the output variable.
    #' @param costs Numeric vector. The costs associated with the features.
    #' @param pop_size Numeric. The population size.
    #' @param num_fea Numeric. The number of features.
    #' @param n_iter Numeric. The number of iterations.
    #' @return An initialized SVMFeature object.
    initialize = function(data, inputs, output, costs, pop_size, num_fea, n_iter = 10, max_time = 300, mode = "iters") {
      .self$data <- data
      .self$inputs <- inputs
      .self$output <- output
      .self$costs <- costs
      .self$pop_size <- pop_size
      .self$num_fea <- num_fea
      .self$n_iter <- n_iter
      .self$max_time <- max_time
      .self$mode <- mode

      # Normalización de los data
      scaler <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (rng[1] == rng[2]) {
          return(rep(0.5, length(x)))  # Si todos los valores son iguales, asignar 0.5 a todos
        } else {
          return((x - rng[1]) / (rng[2] - rng[1]))
        }
      }

      norm_data <- as.data.frame(lapply(.self$data[, inputs], scaler))
      norm_data[output] <- .self$data[[output]]

      # Depuramos los data cargados
      cat("Normalized data:\n")
      print(head(norm_data))

      .self$population <- SVMFeature::Population(data = norm_data, costs = .self$costs, pop_size = .self$pop_size, inputs = .self$inputs, output = .self$output, num_features = .self$num_fea)
      .self$best_population <- NULL
    },

    #' Run the SVMFeature optimization process
    #'
    #' @return The final population after the optimization process.
    run = function() {
      .self$population <- generate_initial_population(.self$population)
      .self$population <- fnds(.self$population)
      .self$update_df_solutions()

      .self$best_population <- .self$population

      init_time <- Sys.time()

      run_iteration <- function() {
        .self$population <- new_population(.self$population)
        .self$population <- fnds(.self$population)

        reduced_population <- reduce_population(.self$population)
        reduced_population <- fnds(reduced_population)

        .self$population$solution_list <- reduced_population$solution_list
        .self$update_df_solutions()
        rm(reduced_population)

        # Comparar y mantener la mejor población basada en el número de soluciones en el frente 1
        current_solutions <- nrow(filter(.self$population$df_solutions, FRONT == 1))
        best_solutions <- nrow(filter(.self$best_population$df_solutions, FRONT == 1))

        if (current_solutions > best_solutions) {
          .self$best_population <- .self$population
        }
      }

      if (.self$mode == "iters") {
        iter <- 0
        while (iter < .self$n_iter) {
          cat(sprintf("ITERATION %d: %f\n", iter, as.numeric(difftime(Sys.time(), init_time, units = "secs"))))
          run_iteration()
          iter <- iter + 1
        }
      } else if (.self$mode == "time") {
        while (as.numeric(difftime(Sys.time(), init_time, units = "secs")) < .self$max_time) {
          cat(sprintf("TIME: %f\n", as.numeric(difftime(Sys.time(), init_time, units = "secs"))))
          run_iteration()
        }
      } else {
        stop("Invalid execution mode. Usage: 'iters' or 'time'.")
      }

      print_population(.self$best_population)

      if (is.null(.self$best_population$df_solutions)) {
        stop("Error: df_solutions is NULL.")
      }

      front_population1 <- filter(.self$best_population$df_solutions, FRONT == 1)

      if (nrow(front_population1) == 0) {
        cat("There is no solution in front 1 to plot.\n")
      } else {
        front_population1$DIST <- as.numeric(front_population1$DIST)
        front_population1$EPS <- as.numeric(front_population1$EPS)

        plot <- ggplot(front_population1, aes(x = DIST, y = EPS)) +
          geom_point(color = 'blue') +
          labs(x = 'Distance', y = 'Epsilon', title = 'Solutions in Front 1') +
          theme_minimal() +
          theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray"))

        print(plot)
      }
    },

    update_df_solutions = function() {
      df_solutions <- do.call(rbind, lapply(.self$population$solution_list, function(solution) {
        solution_dict <- to_dict(solution)
        solution_dict <- lapply(solution_dict, function(item) if (is.list(item)) unlist(item) else item)
        solution_df <- as.data.frame(t(solution_dict), stringsAsFactors = FALSE)
        return(solution_df)
      }))
      df_solutions$FRONT <- as.numeric(as.character(df_solutions$FRONT))
      df_solutions <- df_solutions[order(df_solutions$FRONT), ]
      .self$population$df_solutions <- df_solutions

      # Imprimir número de soluciones por frente
      fronts_count <- df_solutions %>% group_by(FRONT) %>% summarize(count = n())
      cat("Number of solutions by front:\n")
      print(fronts_count)
    }
  )
)
