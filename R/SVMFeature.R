utils::globalVariables(c("FRONT"))
utils::globalVariables(c("DIST"))
utils::globalVariables(c("EPS"))

#' SVMFeature Class (S3 Implementation)
#'
#' Implements a feature selection algorithm using Support Vector Machines (SVM).
#'
#' @param data DataFrame. The input data.
#' @param inputs Character vector. The names of the input features.
#' @param output Character. The name of the output variable.
#' @param costs Numeric vector. The costs associated with the features.
#' @param pop_size Numeric. The population size.
#' @param num_fea Numeric. The number of features.
#' @param n_iter Numeric. The number of iterations (default: 10).
#' @param max_time Numeric. The maximum execution time (default: 300 seconds).
#' @param mode Character. "iters" for iteration-based execution, "time" for time-based execution.
#' @return An S3 object of class "SVMFeature".
#' @export
#'
#' @importFrom utils head
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize n
SVMFeature <- function(data, inputs, output, costs, pop_size, num_fea,
                       n_iter = 10, max_time = 300, mode = "iters") {

  if (!is.data.frame(data)) stop('Error: "data" must be a data frame')
  if (!mode %in% c("iters", "time")) stop('Error: "mode" must be either "iters" or "time"')
  if (!is.character(inputs) || !all(inputs %in% colnames(data))) stop('Error: "inputs" must be a character vector of existing column names in "data"')
  if (!is.character(output) || length(output) != 1 || !(output %in% colnames(data))) stop('Error: "output" must be a single column name in "data"')
  if (!is.numeric(max_time) || max_time <= 0) stop('Error: "max_time" must be a positive number')
  if (!is.numeric(n_iter) || n_iter <= 0) stop('Error: "n_iter" must be a positive number')
  if (!is.numeric(pop_size) || pop_size <= 0) stop('Error: "pop_size" must be a positive number')
  if (!is.numeric(num_fea) || num_fea <= 0) stop('Error: "num_fea" must be a positive number')

  # Normalization function
  scaler <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) return(rep(0.5, length(x)))
    return((x - rng[1]) / (rng[2] - rng[1]))
  }

  # Normalize input data
  norm_data <- as.data.frame(lapply(data[, inputs], scaler))
  norm_data[[output]] <- data[[output]]

  # Print normalized data
  cat("Normalized data:\n")
  print(head(norm_data))

  # Create the S3 object (a list with attributes)
  object <- list(
    data = data,
    inputs = inputs,
    output = output,
    costs = costs,
    pop_size = pop_size,
    num_fea = num_fea,
    n_iter = n_iter,
    max_time = max_time,
    mode = mode,
    population = Population(data = norm_data, costs = costs, pop_size = pop_size,
                            inputs = inputs, output = output, num_features = num_fea),
    best_population = NULL
  )

  class(object) <- "SVMFeature"
  return(object)
}

#' Run the SVMFeature optimization process
#' @param object An object of class "SVMFeature".
#' @return The final population after the optimization process.
#' @export
run.SVMFeature <- function(object) {

  object$population <- generate_initial_population(object$population)
  object$population <- fnds(object$population)
  object <- update_df_solutions.SVMFeature(object)

  object$best_population <- object$population
  init_time <- Sys.time()

  run_iteration <- function() {
    object$population <- new_population(object$population)
    object$population <- fnds(object$population)

    reduced_population <- reduce_population(object$population)
    reduced_population <- fnds(reduced_population)

    object$population$solution_list <- reduced_population$solution_list
    object <- update_df_solutions.SVMFeature(object)

    # Compare and keep the best population based on FRONT == 1
    current_solutions <- nrow(dplyr::filter(object$population$df_solutions, FRONT == 1))
    best_solutions <- nrow(dplyr::filter(object$best_population$df_solutions, FRONT == 1))

    if (current_solutions > best_solutions) {
      object$best_population <- object$population
    }
  }

  if (object$mode == "iters") {
    for (iter in seq_len(object$n_iter)) {
      cat(sprintf("ITERATION %d: %f seconds\n", iter, as.numeric(difftime(Sys.time(), init_time, units = "secs"))))
      run_iteration()
    }
  } else if (object$mode == "time") {
    while (as.numeric(difftime(Sys.time(), init_time, units = "secs")) < object$max_time) {
      cat(sprintf("TIME: %f seconds\n", as.numeric(difftime(Sys.time(), init_time, units = "secs"))))
      run_iteration()
    }
  } else {
    stop("Invalid execution mode. Use 'iters' or 'time'.")
  }

  print_population(object$best_population)

  # Plot the Pareto Front if solutions exist
  if (!is.null(object$best_population$df_solutions)) {
    front_population1 <- dplyr::filter(object$best_population$df_solutions, FRONT == 1)

    if (nrow(front_population1) > 0) {
      front_population1$DIST <- as.numeric(front_population1$DIST)
      front_population1$EPS <- as.numeric(front_population1$EPS)

      plot <- ggplot2::ggplot(front_population1, ggplot2::aes(x = DIST, y = EPS)) +
        ggplot2::geom_point(color = 'blue') +
        ggplot2::labs(x = 'Distance', y = 'Epsilon', title = 'Solutions in Front 1') +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid', colour = "gray"))

      print(plot)
    } else {
      cat("There is no solution in front 1 to plot.\n")
    }
  }

  return(object)
}

#' Update df_solutions for an SVMFeature object
#' @param object An object of class "SVMFeature".
#' @return The updated SVMFeature object.
#' @export
update_df_solutions.SVMFeature <- function(object) {

  df_solutions <- do.call(rbind, lapply(object$population$solution_list, function(solution) {
    solution_dict <- to_dict(solution)
    solution_dict <- lapply(solution_dict, function(item) if (is.list(item)) unlist(item) else item)
    as.data.frame(t(solution_dict), stringsAsFactors = FALSE)
  }))

  df_solutions$FRONT <- as.numeric(as.character(df_solutions$FRONT))
  df_solutions <- df_solutions[order(df_solutions$FRONT), ]
  object$population$df_solutions <- df_solutions

  # Print the number of solutions per front
  fronts_count <- df_solutions %>% dplyr::group_by(FRONT) %>% dplyr::summarize(count = dplyr::n())
  cat("Number of solutions by front:\n")
  print(fronts_count)

  return(object)
}
