#' @title Definition of the Solution Class
#'
#' @description
#'   Defines the `Solution` class, which represents an individual solution.
#'
#' @param num Solution number
#' @param data Data set
#' @param costs Cost vector
#' @param inputs Names of input variables
#' @param output Name of the output variable
#' @param num_features Number of features
#'
#' @return Object of the `Solution` class
#'
#' @export
Solution <- function(num, data, costs, inputs, output, num_features) {
  # Definir variables de instancia
  solution <- list(
    num = num,
    data = data,
    costs = costs,
    inputs = inputs,
    output = output,
    num_dim = length(inputs),
    num_features = num_features,
    data_sol = NULL,
    features = list(),
    vectors = list(),
    plane_coord = list(),
    plane_term_b = list(),
    objective = numeric(2),  # Index 1 para distancia, Index 2 para epsilon
    sol_dom_by = 0,  # Inicializar contador de dominación
    dominates_list = list(),  # Lista de soluciones que esta solución domina
    list_dominated_by = list(),  # Lista de soluciones que dominan a esta solución
    front = -1,
    successful_evaluation = FALSE,
    crowding_distance = 0
  )
  return(solution)
}



#' @title Method to convert the solution to a dictionary
#'
#' @description
#'   Method to convert the solution to a dictionary
#'
#' @param solution Solution class object
#'
#' @return List with the information of the solution
#'
#' @export
to_dict <- function(solution) {
  return(list(
    #SOL = solution$num,
    VECTORS = solution$vectors,
    #PLANO_COOR = solution$plane_coord,
    #PLANO_COOR_B = solution$plane_term_b,
    FEATURES = solution$features,
    DIST = solution$objective[1],  # Distancia
    EPS = solution$objective[2],   # Epsilon
    # COSTE = solution$obj_coste,  # Descomentar si obj_coste es utilizado
    # MC+ = solution$obj_mal_clasificados$mc_pos,  # Descomentar si se utiliza
    # MC- = solution$obj_mal_clasificados$mc_neg,  # Descomentar si se utiliza
    DOMINATES_TO = solution$dominates_list,
    DOMINATED_BY = solution$list_dominated_by,
    SOL_DOM_BY = solution$sol_dom_by,
    FRONT = solution$front,
    CROW_DIST = solution$crowding_distance
  ))
}



#' @title Method to Obtain a Random Class Vector
#'
#' @description
#'   Method to obtain a random class vector from the dataset.
#'
#' @param solution Solution class object containing the dataset
#' @param clazz Class for which to obtain a vector
#'
#' @return Index of the random class vector
#' @export
get_class_vector <- function(solution, clazz) {
  indices <- which(solution$data[[solution$output]] == clazz)
  if (length(indices) > 0) {
    return(sample(indices, 1))
  } else {
    return(NA)  # NA es un valor numérico en R, que indica 'no disponible'
  }
}



#' @title Method to Generate a Random Solution
#'
#' @description
#'   Method to generate a random solution by selecting random features and class vectors.
#'
#' @param solution Solution class object to modify
#'
#' @importFrom stats runif
#'
#' @return Solution class object with a randomly generated solution
#'
#' @export
generate_random_solution <- function(solution) {
  # Seleccionar dos puntos de la clase A y B
  classes <- sort(unique(solution$data[[solution$output]]))
  solution$vectors <- vector("list", length(classes))

  for (i in seq_along(classes)) {
    solution$vectors[[i]] <- get_class_vector(solution, classes[i])
  }

  # Seleccionar aleatoriamente p características y p coordenadas entre [-1, 1]
  if (solution$num_features == solution$num_dim) {
    solution$features <- solution$inputs
  } else {
    possible_features <- solution$inputs[which(solution$inputs %in% names(solution$data))] # Se verifica que las entradas existan como columnas
    solution$features <- sample(possible_features, solution$num_features)
  }

  # Calcular las coordenadas del plano
  solution$plane_coord <- runif(solution$num_features, -1, 1) # Coordenadas aleatorias para cada característica

  # Ajustar dataframe de solución
  indices <- unlist(solution$vectors)
  if (length(indices) == 0 || length(solution$features) == 0) {
    return(solution) # Retorna la solución sin modificar si no hay índices válidos o características seleccionadas
  }

  solution$data_sol <- solution$data[indices, solution$features, drop = FALSE]

  return(solution)
}



#' @title Method to Construct Planes
#'
#' @description
#'   Method to construct planes by calculating the coefficients based on features.
#'
#' @param solution Solution class object containing data and features to construct planes
#'
#' @return Solution class object with constructed planes
#'
#' @export
construct_planes <- function(solution) {
  # Preparar el dataframe data_sol seleccionando solo las filas y columnas necesarias
  valid_indices <- unlist(solution$vectors)
  valid_indices <- valid_indices[valid_indices <= nrow(solution$data) & valid_indices > 0]

  if (length(valid_indices) == 0) {
    stop("No valid index in 'vectors'.")
  }

  # Filtrar las características para asegurarse de que están presentes en los datos
  valid_features <- solution$features[solution$features %in% names(solution$data)]
  if (length(valid_features) == 0) {
    stop("No valid feature for selection.")
  }

  # Subseleccionar el dataframe
  solution$data_sol <- solution$data[valid_indices, valid_features, drop = FALSE]

  # Inicializar el vector de términos del plano
  solution$plane_term_b <- numeric()

  # Calcular los términos independientes de la ecuación del plano
  for (i in seq_len(nrow(solution$data_sol))) {
    independent_term <- 0
    for (j in seq_along(valid_features)) {
      independent_term <- independent_term + solution$data_sol[i, j] * solution$plane_coord[j]
    }
    # Añadir el valor negativo de 'independent_term' al vector de términos del plano
    solution$plane_term_b <- c(solution$plane_term_b, -independent_term)
  }

  # Añadir el término del plano intermedio
  if(length(solution$plane_term_b) > 0) {
    middle <- mean(solution$plane_term_b)
    solution$plane_term_b <- c(solution$plane_term_b, middle)
  } else {
    stop("No plane terms have been calculated, so the mean can not be calculated.")
  }

  return(solution)
}



#' @title Method to Calculate the Distance Objective
#'
#' @description
#'   This method calculates the distance objective based on the plane equation.
#'
#' @param solution Solution class object containing the plane equation
#'
#' @return Solution class object with the distance objective calculated
#'
#' @export
calculate_distance_objective <- function(solution) {
  # Inicializar el objective de distancia a cero
  solution$objective[[1]] <- 0
  denominator <- sum(solution$plane_coord^2) ^ 0.5
  distance <- abs(solution$plane_term_b[[1]] - solution$plane_term_b[[2]])

  # Calcular la distancia normalizada si el denominador no es cero
  if (denominator != 0) {
    solution$objective[[1]] <- distance / denominator
  } else {
    solution$objective[[1]] <- -1  # Valor de error o indicativo de no calculable
  }
  return(solution)
}



#' @title Method to Calculate the Epsilon Objective
#'
#' @description
#'   This method calculates the epsilon objective by evaluating misclassified points.
#'
#' @param solution Solution class object containing data and features for classification
#'
#' @return Solution class object with the epsilon objective calculated
#'
#' @export
calculate_epsilon_objective <- function(solution) {
  # Inicializar contadores de clasificaciones malas
  mc_pos <- 0
  mc_neg <- 0
  solution$objective[[2]] <- 0

  # Iterar sobre todos los datos para calcular la suma de errores
  for (i in 1:nrow(solution$data)) {
    clazz <- solution$data[i, solution$output]
    distance <- 0

    # Calcular la distancia lineal usando los coeficientes
    for (index in seq_along(solution$features)) {
      fea <- solution$features[[index]]
      distance <- distance + solution$plane_coord[index] * solution$data[i, fea]
    }

    # Ajustar por el término independiente y contar errores
    if (clazz == -1) {
      distance <- distance + solution$plane_term_b[[1]]
      if (distance < 0) {
        mc_neg <- mc_neg + 1
        solution$objective[[2]] <- solution$objective[[2]] + abs(distance)
      }
    } else {
      distance <- distance + solution$plane_term_b[[2]]
      if (distance > 0) {
        mc_pos <- mc_pos + 1
        solution$objective[[2]] <- solution$objective[[2]] + abs(distance)
      }
    }
  }

  # Opcionalmente, podría almacenar los conteos de clasificación incorrecta si fuera necesario
  solution$mc_pos <- mc_pos
  solution$mc_neg <- mc_neg

  return(solution)
}



# calcular_objetivo_costes <- function(solucion) {
#   solucion$obj_coste <- 0
#   for (feature in solucion$features) {
#     solucion$obj_coste <- solucion$obj_coste + solucion$costs[feature]
#   }
#   return(solucion)
# }



#' @title Evaluate Solution
#'
#' @description
#'   This function evaluates a solution by calculating its objectives.
#'
#' @param solution Solution class object to be evaluated
#'
#' @return Solution class object with updated evaluation status and objectives
#'
#' @export
evaluate_solution <- function(solution) {
  # Verifica si la suma de plane_coord es cero
  if (sum(solution$plane_coord) == 0) {
    cat("The solution cannot be evaluated, there are no coordinates.\n")
    solution$successful_evaluation <- FALSE
    return(solution)
  } else {
    # Procedimientos para evaluar la solución si hay coordenadas válidas
    solution <- construct_planes(solution)
    solution <- calculate_distance_objective(solution)
    solution <- calculate_epsilon_objective(solution)
    # solution <- calcular_objetivo_costes(solution) si es necesario

    # Asumir que las evaluaciones fueron exitosas si se alcanza este punto
    solution$successful_evaluation <- TRUE
    return(solution)
  }
}



#' @title Method to Determine Dominance
#'
#' @description
#'   This function determines if one solution dominates another based on objectives.
#'
#' @param solution1 First solution object
#' @param solution2 Second solution object
#'
#' @return TRUE if solucion1 dominates solucion2, otherwise FALSE
#'
#' @export
dominate <- function(solution1, solution2) {
  # Comprobar si alguno de los objetivos es NA o no numérico antes de comparar
  if (is.na(solution1$objective[[1]]) || !is.numeric(solution1$objective[[1]]) ||
      is.na(solution2$objective[[1]]) || !is.numeric(solution2$objective[[1]]) ||
      is.na(solution1$objective[[2]]) || !is.numeric(solution1$objective[[2]]) ||
      is.na(solution2$objective[[2]]) || !is.numeric(solution2$objective[[2]])) {
    # Retorna FALSE o maneja el caso de forma que se ajuste a tu lógica de aplicación
    return(FALSE)
  }

  # Procede con la comparación si todos los valores son adecuados
  return(solution1$objective[[1]] >= solution2$objective[[1]] &&
           solution1$objective[[2]] <= solution2$objective[[2]])
}



#' @title Method to Determine Dominance with Three States
#'
#' @description
#'   This function determines the dominance state between two solutions.
#'
#' @param solution1 First solution object
#' @param solution2 Second solution object
#'
#' @return 0 if the solutions are equal, 1 if solucion1 dominates solucion2, 2 if solucion2 dominates solucion1
#'
#' @export
dominate2 <- function(solution1, solution2) {
  # Comprobar si los objetivos son iguales
  if (solution1$objective[[1]] == solution2$objective[[1]] && solution1$objective[[2]] == solution2$objective[[2]]) {
    return(0)  # No hay dominancia, son iguales
  } else if (solution1$objective[[1]] >= solution2$objective[[1]] && solution1$objective[[2]] <= solution2$objective[[2]]) {
    return(1)  # solucion1 domina a solucion2
  } else {
    return(2)  # solucion2 domina a solucion1
  }
}



#' @title Method to Compare Solutions
#'
#' @description
#'   This function compares two solutions based on their front and dominance information.
#'
#' @param solution1 First solution object
#' @param solution2 Second solution object
#'
#' @return The solution object with better dominance or front
#'
#' @export
compare_solutions <- function(solution1, solution2) {
  if (solution1$front == solution2$front) {
    if (length(solution1$list_dominated_by) == length(solution2$list_dominated_by)) {
      coin_flip <- sample(0:1, 1)
      if (coin_flip == 0) {
        return(solution1)
      } else {
        return(solution2)
      }
    } else if (length(solution1$list_dominated_by) < length(solution2$list_dominated_by)) {
      return(solution1)
    } else {
      return(solution2)
    }
  } else {
    if (solution1$front > solution2$front) {
      return(solution2)
    } else {
      return(solution1)
    }
  }
}



#' @title Method to Mutate Vectors
#'
#' @description
#'   This function mutates the vectors in the solution by assigning new class vectors.
#'
#' @param solution Solution class object to mutate
#'
#' @return Solution class object with updated vectors
#'
#' @export
mutate_vectors <- function(solution) {
  # Obtener las clases únicas y ordenarlas
  classes <- sort(unique(solution$data[[solution$output]]))

  # Iterar sobre cada clase
  for (i in seq_along(classes)) {
    # Llamar a get_class_vector para la clase actual
    vector <- get_class_vector(solution, classes[i])

    # Verificar que el nuevo vector sea diferente al actual y si no, buscar otro
    while (identical(vector, solution$vectors[[i]])) {
      vector <- get_class_vector(solution, classes[i])
    }

    # Actualizar el vector en la solución
    solution$vectors[[i]] <- vector
  }

  return(solution)
}
