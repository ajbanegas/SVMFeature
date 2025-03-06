#' @title Definition of the Population Class
#'
#' @description
#'   This function defines the `Population` class with its attributes and methods.
#'
#' @param data Dataset
#' @param costs Associated costs
#' @param pop_size Population size
#' @param inputs Input variables
#' @param output Output variable
#' @param num_features Number of features
#' @param num_obj Number of objective (default: 2)
#' @param clones Number of clones (default: 0)
#' @param p_mutation Mutation probability (default: 0.7)
#' @param p_mut_ind Individual mutation probability (default: 0.4)
#' @param p_mut_fea Feature mutation probability (default: 0.4)
#' @param p_mut_coord Mutation coordinates (default: 0.2)
#' @param mut_coord Mutation coordinates (default: 0)
#'
#' @return Object of the `Population` class
#'
#' @export
Population <- function(data, costs, pop_size, inputs, output, num_features, num_obj = 2, clones = 0,
                      p_mutation = 0.7, p_mut_ind = 0.4, p_mut_fea = 0.4, p_mut_coord = 0.2, mut_coord = 0) {
  population <- list(
    num_features = num_features,
    inputs = inputs,
    output = output,
    num_dim = length(inputs),
    data = data,
    costs = costs,
    pop_size = pop_size,
    clones = clones,
    solution_list = list(),
    num_fronts = NULL,
    num_sol_fron = NULL,
    df_solutions = NULL,
    fronts = list(),
    p_mutation = p_mutation,
    p_mut_ind = p_mut_ind,
    p_mut_fea = p_mut_fea,
    p_mut_coord = p_mut_coord,
    mut_coord = mut_coord,
    num_obj = num_obj
  )
  return(population)
}



#' @title Generate initial population
#'
#' @description
#'   This function generates the initial population.
#'
#' @param population Population class object
#'
#' @return Population class object updated with the generated initial population
#'
#' @export
generate_initial_population <- function(population) {
  cat(sprintf("Creating initial population of size %d\n", population$pop_size))
  num_sol <- 0

  while (length(population$solution_list) < population$pop_size) {
    # Crear y generar solución aleatoria
    sol <- Solution(num_sol, population$data, population$costs, population$inputs, population$output, population$num_features)
    sol <- generate_random_solution(sol)
    sol <- evaluate_solution(sol)

    # Comprobar si la solución ha sido evaluada exitosamente
    if (sol$successful_evaluation) {
      # Comprobación adicional para clones si necesario
      if (population$clones == 1 || check_clones(population, sol) == 0) {
        # Añadir la solución evaluada y actualizada a la población
        population$solution_list[[num_sol + 1]] <- sol
        num_sol <- num_sol + 1
      } else {
        cat("Detected and skipped clone.\n")
      }
    } else {
      cat("Disregarded invalid solution.\n")
    }
  }
  return(population)
}



#' @title Print population
#'
#' @description
#'   This function prints the population.
#'
#' @param population Population class object
#'
#' @return Prints the population in the console
#'
#' @export
print_population <- function(population) {
  # Convertir lista de soluciones a data frame
  df_solutions <- do.call(rbind, lapply(population$solution_list, function(solution) {
    solution_dict <- to_dict(solution)
    # Asegurarse de que todos los valores sean vectores y no listas
    solution_dict <- lapply(solution_dict, function(item) if(is.list(item)) unlist(item) else item)
    # Convertir a dataframe
    solution_df <- as.data.frame(t(solution_dict), stringsAsFactors = FALSE)
    return(solution_df)
  }))

  # Convertir las columnas necesarias a numérico
  df_solutions$FRONT <- as.numeric(as.character(df_solutions$FRONT))

  # Verificar si la columna FRONT es numérica ahora
  if(!is.numeric(df_solutions$FRONT)) {
    stop("FRONT column could not be converted to numeric.")
  }

  # Ordenar el dataframe por la columna 'FRONT'
  df_solutions <- df_solutions[order(df_solutions$FRONT), ]

  print(df_solutions)
}



#' @title Check for Clones
#'
#' @description
#'   This function checks if a solution already exists in the solution list to avoid duplicates.
#'
#' @param population Population
#' @param solucion_eva Solution object to be evaluated for duplication
#'
#' @return 1 if a clone is found, otherwise 0
#'
#' @export
check_clones <- function(population, solution_eva) {
  clone <- 0
  for (sol in population$solution_list) {
    # Comparar todas las características relevantes
    if (all(sol$features == solution_eva$features) &&
        abs(sol$objective[[1]] - solution_eva$objective[[1]]) < 1e-5 &&
        abs(sol$objective[[2]] - solution_eva$objective[[2]]) < 1e-5 &&
        all(abs(sol$plane_coord - solution_eva$plane_coord) < 1e-5) &&
        identical(sol$vectors, solution_eva$vectors)) {
      clone <- 1
      break
    } else {
      clone <- 0
    }
  }
  return(clone)
}



#' @title Fast Non-Dominated Sorting (fnds)
#'
#' @description
#'   This function performs the Fast Non-Dominated Sorting on a given population.
#'
#' @param population Population class object
#'
#' @return Updated Population class object with sorted fronts
#'
#' @export
fnds <- function(population) {
  pop_size <- length(population$solution_list)
  fronts <- list()
  front <- list(solutions = numeric())
  # Inicializar soluciones
  for (i in 1:pop_size) {
    population$solution_list[[i]]$sol_dom_by <- 0
    population$solution_list[[i]]$dominates_list <- list()
    population$solution_list[[i]]$list_dominated_by <- list()
    population$solution_list[[i]]$front <- -1
  }

  # Calcular dominancia
  for (n in 1:(pop_size - 1)) {
    for (s in (n + 1):pop_size) {
      if (dominate(population$solution_list[[n]], population$solution_list[[s]])) {
        population$solution_list[[n]]$dominates_list <- unique(c(population$solution_list[[n]]$dominates_list, s))
        population$solution_list[[s]]$list_dominated_by <- unique(c(population$solution_list[[s]]$list_dominated_by, n))
        population$solution_list[[s]]$sol_dom_by <- population$solution_list[[s]]$sol_dom_by + 1
      } else if (dominate(population$solution_list[[s]], population$solution_list[[n]])) {
        population$solution_list[[s]]$dominates_list <- unique(c(population$solution_list[[s]]$dominates_list, n))
        population$solution_list[[n]]$list_dominated_by <- unique(c(population$solution_list[[n]]$list_dominated_by, s))
        population$solution_list[[n]]$sol_dom_by <- population$solution_list[[n]]$sol_dom_by + 1
      }
    }
  }
  # Asignar frente 1 a las soluciones no dominadas
  for (i in 1:pop_size) {
    if (population$solution_list[[i]]$sol_dom_by == 0) {
      front$solutions <- c(front$solutions, i)
      population$solution_list[[i]]$front <- 1
    }
  }
  fronts[[1]] <- front

  # Crear frentes subsiguientes
  r <- 1
  while (length(fronts[[r]]$solutions) > 0) {
    new_front <- numeric()
    for (sol_num in fronts[[r]]$solutions) {
      dominated <- population$solution_list[[sol_num]]$dominates_list
      for (i in dominated) {
        population$solution_list[[i]]$sol_dom_by <- population$solution_list[[i]]$sol_dom_by - 1
        if (population$solution_list[[i]]$sol_dom_by == 0) {
          new_front <- c(new_front, i)
          population$solution_list[[i]]$front <- r + 1
        }
      }
    }
    r <- r + 1
    front <- list(solutions = new_front)
    fronts[[r]] <- front
  }

  population$fronts <- fronts
  return(population)
}



#' @title Tournament Selection for Parent
#'
#' @description
#'   This function selects a parent solution from the population using tournament selection.
#'
#' @param population Population class object
#'
#' @return Selected parent solution
#'
#' @export
tournament_select_parent <- function(population) {
  # Seleccionar dos candidatos al azar
  num1 <- sample(seq_len(length(population$solution_list)), 1)
  num2 <- sample(seq_len(length(population$solution_list)), 1)

  # Asegurarse de que no sean el mismo candidato
  while (num1 == num2) {
    num2 <- sample(seq_len(length(population$solution_list)), 1)
  }

  # Extraer los candidatos directamente
  parent1 <- population$solution_list[[num1]]
  parent2 <- population$solution_list[[num2]]

  # Comparar los dos padres y devolver el "mejor"
  return(compare_solutions(parent1, parent2))
}



#' @title Crossover Solutions
#'
#' @description
#'   This function performs crossover between two parent solutions to create offspring.
#'
#' @param population Population class object
#' @param parent First parent solution
#' @param mother Second parent solution
#' @param num Unique identifier for the new solutions
#'
#' @importFrom stats runif
#'
#' @return List of two new offspring solutions
#'
#' @export
crossover_solutions <- function(population, parent, mother, num) {
  # Crear 4 hijos (soluciones)
  child <- Solution(num, population$data, population$costs, population$inputs, population$output, population$num_features)
  child$features <- parent$features
  child$plane_coord <- parent$plane_coord
  child$vectors <- list(parent$vectors[[1]], mother$vectors[[2]])

  child1 <- Solution(num, population$data, population$costs, population$inputs, population$output, population$num_features)
  child1$features <- mother$features
  child1$plane_coord <- mother$plane_coord
  child1$vectors <- list(parent$vectors[[1]], mother$vectors[[2]])

  child2 <- Solution(num, population$data, population$costs, population$inputs, population$output, population$num_features)
  child2$features <- parent$features
  child2$plane_coord <- parent$plane_coord
  child2$vectors <- list(mother$vectors[[1]], parent$vectors[[2]])

  child3 <- Solution(num, population$data, population$costs, population$inputs, population$output, population$num_features)
  child3$features <- mother$features
  child3$plane_coord <- mother$plane_coord
  child3$vectors <- list(mother$vectors[[1]], parent$vectors[[2]])

  # Evaluar las soluciones
  child <- evaluate_solution(child)
  child1 <- evaluate_solution(child1)
  child2 <- evaluate_solution(child2)
  child3 <- evaluate_solution(child3)

  # Comparar soluciones y determinar el ganador y el perdedor
  domina <- dominate2(child, child1)
  losing_child <- child
  if (domina == 1) {
    losing_child <- child1
  } else if (domina == 0 && runif(1) < 0.5) {
    losing_child <- child1
  }

  dominates <- dominate2(child2, child3)
  winning_child <- child2
  if (dominates == 2) {
    winning_child <- child3
  } else if (dominates == 0 && runif(1) < 0.5) {
    winning_child <- child3
  }

  # Actualizar el perdedor con los atributos del ganador
  losing_child$features <- winning_child$features
  losing_child$plane_coord <- winning_child$plane_coord
  losing_child$vectors[[1]] <- winning_child$vectors[[1]]
  losing_child$vectors[[2]] <- winning_child$vectors[[2]]

  losing_child$num <- num
  winning_child$num <- num + 1

  return(list(losing_child, winning_child))
}



#' @title Create New Population
#'
#' @description
#'   This function creates a new population through crossover.
#'
#' @param population Population class object
#'
#' @importFrom stats runif
#'
#' @return Updated Population class object with the new population generated
#'
#' @export
new_population <- function(population) {
  cat(sprintf("Creating new population by crossover from %d to %d....\n", population$pop_size, 2 * population$pop_size))
  i <- population$pop_size

  while (length(population$solution_list) < 2 * population$pop_size) {
    parent <- tournament_select_parent(population)
    mother <- tournament_select_parent(population)
    while (parent$num == mother$num) {
      mother <- tournament_select_parent(population)
    }

    offspring <- crossover_solutions(population, parent, mother, i)
    for (child in offspring) {
      if (length(population$solution_list) >= 2 * population$pop_size) {
        break  # Evitar agregar más hijos si se alcanza el tamaño máximo
      }

      if (runif(1) < population$p_mutation) {
        child <- mutate_solution(child, population)
      }

      child <- evaluate_solution(child)
      while (!child$successful_evaluation) {
        child <- generate_random_solution(child)
        child <- evaluate_solution(child)
      }

      while (check_clones(population, child) == 1) {
        child <- mutate_solution(child, population)
        child <- evaluate_solution(child)
      }

      population$solution_list[[length(population$solution_list) + 1]] <- child
      i <- i + 1
    }
  }

  return(population)
}



#' @title Mutate Solution
#'
#' @description
#'   This function mutates a solution based on mutation probabilities and coordinates.
#'
#' @param solution Solution object to mutate
#' @param population Population class object that contains mutation parameters
#'
#' @importFrom stats runif
#'
#' @return Mutated solution object
#'
#' @export
mutate_solution <- function(solution, population) {
  random <- runif(1)

  # Mutación de vectores
  if (random < population$p_mut_ind) {
    solution <- mutate_vectors(solution)
  }

  # Mutación de características
  if (population$num_dim > population$num_features) {
    for (i in seq_along(solution$features)) {
      rand_val <- runif(1)
      if (rand_val < population$p_mut_fea) {
        possible_features <- setdiff(population$inputs, solution$features)
        solution$features[i] <- sample(possible_features, 1)
        solution$plane_coord[i] <- runif(1, -1, 1)
      }
    }
  }

  # Mutación de planos
  for (i in seq_along(solution$plane_coord)) {
    rand_val <- runif(1)
    if (rand_val < population$p_mut_coord) {
      if (population$mut_coord == 0) {
        percentage <- runif(1, 0, 0.25)
        rand_direction <- runif(1)
        if (rand_direction < 0.5) {  # incrementa
          solution$plane_coord[i] <- solution$plane_coord[i] * (1 + percentage)
          if (solution$plane_coord[i] > 1) {
            solution$plane_coord[i] <- solution$plane_coord[i] / 10
          }
        } else {  # decrementa
          solution$plane_coord[i] <- solution$plane_coord[i] * (1 - percentage)
          if (solution$plane_coord[i] < -1) {
            solution$plane_coord[i] <- solution$plane_coord[i] / 10
          }
        }
      } else {
        solution$plane_coord[i] <- runif(1, -1, 1)
      }
    }
  }

  return(solution)
}



#' @title Reduce Population Size
#'
#' @description
#'   This function reduces the population size based on fronts and crowding distance.
#'
#' @param population Population class object to reduce
#'
#' @return Reduced Population class object
#'
#' @export
reduce_population <- function(population) {
  cat("Size reduction from 2N to N....\n")

  # Crear nueva población vacía
  new_population <- list(
    data = population$data,
    costs = population$costs,
    pop_size = population$pop_size,
    inputs = population$inputs,
    output = population$output,
    num_features = population$num_features,
    solution_list = list()
  )

  # Contador de frentes
  front_counter <- 1

  # Añadir soluciones de los frentes a la nueva población hasta completar el tamaño
  while (length(new_population$solution_list) + length(population$fronts[[front_counter]]$solutions) <= population$pop_size) {
    front_solutions <- lapply(population$solution_list, function(sol) {
      if (sol$front == front_counter) return(sol) else return(NULL)
    })
    front_solutions <- front_solutions[!sapply(front_solutions, is.null)]
    new_population$solution_list <- c(new_population$solution_list, front_solutions)
    cat(sprintf("FRONT %d: %d solutions\n", front_counter, length(new_population$solution_list)))
    front_counter <- front_counter + 1
  }

  # Si no se ha completado la población, seleccionar las mejores soluciones restantes usando crowding distance
  if (length(new_population$solution_list) < population$pop_size && front_counter <= length(population$fronts)) {
    front_solutions <- crowding_distance(population, front_counter)
    remaining <- population$pop_size - length(new_population$solution_list)
    new_population$solution_list <- c(new_population$solution_list, front_solutions[1:remaining])
  }

  # Renombrar las soluciones
  for (i in seq_along(new_population$solution_list)) {
    new_population$solution_list[[i]]$num <- i - 1
  }

  return(new_population)
}



#' @title Crowding Distance Calculation
#'
#' @description
#'   This function calculates the crowding distance for solutions in a specified front.
#'
#' @param population Population class object
#' @param num_front The front number to calculate the crowding distance for
#'
#' @return List of solutions with updated crowding distances
#'
#' @export
crowding_distance <- function(population, num_front) {
  # Extraer soluciones del frente especificado
  front_solutions <- lapply(population$solution_list, function(sol) {
    if (sol$front == num_front) return(sol) else return(NULL)
  })
  front_solutions <- front_solutions[!sapply(front_solutions, is.null)]

  num_sol_front <- length(front_solutions)

  # Inicializar la distancia de hacinamiento a cero para todas las soluciones
  for (sol in front_solutions) {
    sol$crowding_distance <- 0
  }

  # Procesar cada objective si hay suficientes soluciones para comparar
  if (num_sol_front > 2) {
    for (obj in 1:population$num_obj) {
      # Ordenar las soluciones basándose en el objective actual
      front_solutions <- front_solutions[order(sapply(front_solutions, function(x) x$objective[[obj]]))]

      # Asignar la distancia de hacinamiento infinita a las soluciones extremas
      front_solutions[[1]]$crowding_distance <- Inf
      front_solutions[[num_sol_front]]$crowding_distance <- Inf

      max_obj <- front_solutions[[num_sol_front]]$objective[[obj]]
      min_obj <- front_solutions[[1]]$objective[[obj]]

      if (max_obj > min_obj) {
        for (j in 2:(num_sol_front - 1)) {
          increment <- (front_solutions[[j + 1]]$objective[[obj]] - front_solutions[[j - 1]]$objective[[obj]]) / (max_obj - min_obj)
          front_solutions[[j]]$crowding_distance <- front_solutions[[j]]$crowding_distance + increment
        }
      }
    }
  }

  # Ordenar las soluciones por la distancia de hacinamiento de mayor a menor
  front_solutions <- front_solutions[order(sapply(front_solutions, function(x) x$crowding_distance), decreasing = TRUE)]

  return(front_solutions)
}
