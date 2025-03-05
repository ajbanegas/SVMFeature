#' @title Definition of the Poblacion Class
#'
#' @description
#'   This function defines the `Poblacion` class with its attributes and methods.
#'
#' @param data Dataset
#' @param costes Associated costs
#' @param tam_pob Population size
#' @param inputs Input variables
#' @param output Output variable
#' @param num_features Number of features
#' @param num_obj Number of objective (default: 2)
#' @param clones Number of clones (default: 0)
#' @param p_mutacion Mutation probability (default: 0.7)
#' @param p_mut_ind Individual mutation probability (default: 0.4)
#' @param p_mut_fea Feature mutation probability (default: 0.4)
#' @param p_mut_coord Mutation coordinates (default: 0.2)
#' @param mut_coord Mutation coordinates (default: 0)
#'
#' @return Object of the `Poblacion` class
#'
#' @export
Poblacion <- function(data, costes, tam_pob, inputs, output, num_features, num_obj = 2, clones = 0,
                      p_mutacion = 0.7, p_mut_ind = 0.4, p_mut_fea = 0.4, p_mut_coord = 0.2, mut_coord = 0) {
  poblacion <- list(
    num_features = num_features,
    inputs = inputs,
    output = output,
    num_dim = length(inputs),
    data = data,
    costes = costes,
    tam_pob = tam_pob,
    clones = clones,
    lista_soluciones = list(),
    num_fronts = NULL,
    num_sol_fron = NULL,
    df_soluciones = NULL,
    fronts = list(),
    p_mutacion = p_mutacion,
    p_mut_ind = p_mut_ind,
    p_mut_fea = p_mut_fea,
    p_mut_coord = p_mut_coord,
    mut_coord = mut_coord,
    num_obj = num_obj
  )
  return(poblacion)
}



#' @title Generate initial population
#'
#' @description
#'   This function generates the initial population.
#'
#' @param poblacion Population class object
#'
#' @return Population class object updated with the generated initial population
generar_poblacion_inicial <- function(poblacion) {
  cat(sprintf("Generando poblacion inicial de tamaño %d\n", poblacion$tam_pob))
  num_sol <- 0

  while (length(poblacion$lista_soluciones) < poblacion$tam_pob) {
    # Crear y generar solución aleatoria
    sol <- Solucion(num_sol, poblacion$data, poblacion$costes, poblacion$inputs, poblacion$output, poblacion$num_features)
    sol <- generar_solucion_aleatoria(sol)
    sol <- evaluar_solucion(sol)

    # Comprobar si la solución ha sido evaluada exitosamente
    if (sol$evaluacion_exitosa) {
      # Comprobación adicional para clones si necesario
      if (poblacion$clones == 1 || comprobar_clones(poblacion, sol) == 0) {
        # Añadir la solución evaluada y actualizada a la población
        poblacion$lista_soluciones[[num_sol + 1]] <- sol
        num_sol <- num_sol + 1
      } else {
        cat("Clon detectado y omitido.\n")
      }
    } else {
      cat("Solución no válida generada y descartada.\n")
    }
  }
  return(poblacion)
}



#' @title Print population
#'
#' @description
#'   This function prints the population.
#'
#' @param poblacion Population class object
#'
#' @return Prints the population in the console
imprimir_poblacion <- function(poblacion) {
  # Convertir lista de soluciones a data frame
  df_soluciones <- do.call(rbind, lapply(poblacion$lista_soluciones, function(solucion) {
    solucion_dict <- to_dict(solucion)
    # Asegurarse de que todos los valores sean vectores y no listas
    solucion_dict <- lapply(solucion_dict, function(item) if(is.list(item)) unlist(item) else item)
    # Convertir a dataframe
    solucion_df <- as.data.frame(t(solucion_dict), stringsAsFactors = FALSE)
    return(solucion_df)
  }))

  # Convertir las columnas necesarias a numérico
  df_soluciones$FRONT <- as.numeric(as.character(df_soluciones$FRONT))

  # Verificar si la columna FRONT es numérica ahora
  if(!is.numeric(df_soluciones$FRONT)) {
    stop("La columna FRONT no se pudo convertir a numérica.")
  }

  # Ordenar el dataframe por la columna 'FRONT'
  df_soluciones <- df_soluciones[order(df_soluciones$FRONT), ]

  print(df_soluciones)
}



#' @title Check for Clones
#'
#' @description
#'   This function checks if a solution already exists in the solution list to avoid duplicates.
#'
#' @param poblacion Population
#' @param solucion_eva Solution object to be evaluated for duplication
#'
#' @return 1 if a clone is found, otherwise 0
comprobar_clones <- function(poblacion, solucion_eva) {
  clon <- 0
  for (sol in poblacion$lista_soluciones) {
    # Comparar todas las características relevantes
    if (all(sol$features == solucion_eva$features) &&
        abs(sol$objetivo[[1]] - solucion_eva$objetivo[[1]]) < 1e-5 &&
        abs(sol$objetivo[[2]] - solucion_eva$objetivo[[2]]) < 1e-5 &&
        all(abs(sol$plano_coord - solucion_eva$plano_coord) < 1e-5) &&
        identical(sol$vectors, solucion_eva$vectors)) {
      clon <- 1
      break
    } else {
      clon <- 0
    }
  }
  return(clon)
}



#' @title Fast Non-Dominated Sorting (fnds)
#'
#' @description
#'   This function performs the Fast Non-Dominated Sorting on a given population.
#'
#' @param poblacion Population class object
#'
#' @return Updated Population class object with sorted fronts
fnds <- function(poblacion) {
  tam_pob <- length(poblacion$lista_soluciones)
  fronts <- list()
  front <- list(soluciones = numeric())
  # Inicializar soluciones
  for (i in 1:tam_pob) {
    poblacion$lista_soluciones[[i]]$sol_dom_por <- 0
    poblacion$lista_soluciones[[i]]$lista_domina_a <- list()
    poblacion$lista_soluciones[[i]]$lista_dominado_por <- list()
    poblacion$lista_soluciones[[i]]$front <- -1
  }

  # Calcular dominancia
  for (n in 1:(tam_pob - 1)) {
    for (s in (n + 1):tam_pob) {
      if (dominar(poblacion$lista_soluciones[[n]], poblacion$lista_soluciones[[s]])) {
        poblacion$lista_soluciones[[n]]$lista_domina_a <- unique(c(poblacion$lista_soluciones[[n]]$lista_domina_a, s))
        poblacion$lista_soluciones[[s]]$lista_dominado_por <- unique(c(poblacion$lista_soluciones[[s]]$lista_dominado_por, n))
        poblacion$lista_soluciones[[s]]$sol_dom_por <- poblacion$lista_soluciones[[s]]$sol_dom_por + 1
      } else if (dominar(poblacion$lista_soluciones[[s]], poblacion$lista_soluciones[[n]])) {
        poblacion$lista_soluciones[[s]]$lista_domina_a <- unique(c(poblacion$lista_soluciones[[s]]$lista_domina_a, n))
        poblacion$lista_soluciones[[n]]$lista_dominado_por <- unique(c(poblacion$lista_soluciones[[n]]$lista_dominado_por, s))
        poblacion$lista_soluciones[[n]]$sol_dom_por <- poblacion$lista_soluciones[[n]]$sol_dom_por + 1
      }
    }
  }
  # Asignar frente 1 a las soluciones no dominadas
  for (i in 1:tam_pob) {
    if (poblacion$lista_soluciones[[i]]$sol_dom_por == 0) {
      front$soluciones <- c(front$soluciones, i)
      poblacion$lista_soluciones[[i]]$front <- 1
    }
  }
  fronts[[1]] <- front

  # Crear frentes subsiguientes
  r <- 1
  while (length(fronts[[r]]$soluciones) > 0) {
    nuevo_frente <- numeric()
    for (sol_num in fronts[[r]]$soluciones) {
      dominados <- poblacion$lista_soluciones[[sol_num]]$lista_domina_a
      for (i in dominados) {
        poblacion$lista_soluciones[[i]]$sol_dom_por <- poblacion$lista_soluciones[[i]]$sol_dom_por - 1
        if (poblacion$lista_soluciones[[i]]$sol_dom_por == 0) {
          nuevo_frente <- c(nuevo_frente, i)
          poblacion$lista_soluciones[[i]]$front <- r + 1
        }
      }
    }
    r <- r + 1
    front <- list(soluciones = nuevo_frente)
    fronts[[r]] <- front
  }

  poblacion$fronts <- fronts
  return(poblacion)
}



#' @title Tournament Selection for Parent
#'
#' @description
#'   This function selects a parent solution from the population using tournament selection.
#'
#' @param poblacion Population class object
#'
#' @return Selected parent solution
torneo_elegir_padre <- function(poblacion) {
  # Seleccionar dos candidatos al azar
  num1 <- sample(seq_len(length(poblacion$lista_soluciones)), 1)
  num2 <- sample(seq_len(length(poblacion$lista_soluciones)), 1)

  # Asegurarse de que no sean el mismo candidato
  while (num1 == num2) {
    num2 <- sample(seq_len(length(poblacion$lista_soluciones)), 1)
  }

  # Extraer los candidatos directamente
  padre1 <- poblacion$lista_soluciones[[num1]]
  padre2 <- poblacion$lista_soluciones[[num2]]

  # Comparar los dos padres y devolver el "mejor"
  return(comparar_soluciones(padre1, padre2))
}



#' @title Crossover Solutions
#'
#' @description
#'   This function performs crossover between two parent solutions to create offspring.
#'
#' @param poblacion Population class object
#' @param padre First parent solution
#' @param madre Second parent solution
#' @param num Unique identifier for the new solutions
#'
#' @return List of two new offspring solutions
cruzar_soluciones <- function(poblacion, padre, madre, num) {
  # Crear 4 hijos (soluciones)
  hijo <- Solucion(num, poblacion$data, poblacion$costes, poblacion$inputs, poblacion$output, poblacion$num_features)
  hijo$features <- padre$features
  hijo$plano_coord <- padre$plano_coord
  hijo$vectors <- list(padre$vectors[[1]], madre$vectors[[2]])

  hijo1 <- Solucion(num, poblacion$data, poblacion$costes, poblacion$inputs, poblacion$output, poblacion$num_features)
  hijo1$features <- madre$features
  hijo1$plano_coord <- madre$plano_coord
  hijo1$vectors <- list(padre$vectors[[1]], madre$vectors[[2]])

  hijo2 <- Solucion(num, poblacion$data, poblacion$costes, poblacion$inputs, poblacion$output, poblacion$num_features)
  hijo2$features <- padre$features
  hijo2$plano_coord <- padre$plano_coord
  hijo2$vectors <- list(madre$vectors[[1]], padre$vectors[[2]])

  hijo3 <- Solucion(num, poblacion$data, poblacion$costes, poblacion$inputs, poblacion$output, poblacion$num_features)
  hijo3$features <- madre$features
  hijo3$plano_coord <- madre$plano_coord
  hijo3$vectors <- list(madre$vectors[[1]], padre$vectors[[2]])

  # Evaluar las soluciones
  hijo <- evaluar_solucion(hijo)
  hijo1 <- evaluar_solucion(hijo1)
  hijo2 <- evaluar_solucion(hijo2)
  hijo3 <- evaluar_solucion(hijo3)

  # Comparar soluciones y determinar el ganador y el perdedor
  domina <- dominar2(hijo, hijo1)
  hijo_perdedor <- hijo
  if (domina == 1) {
    hijo_perdedor <- hijo1
  } else if (domina == 0 && runif(1) < 0.5) {
    hijo_perdedor <- hijo1
  }

  domina <- dominar2(hijo2, hijo3)
  hijo_ganador <- hijo2
  if (domina == 2) {
    hijo_ganador <- hijo3
  } else if (domina == 0 && runif(1) < 0.5) {
    hijo_ganador <- hijo3
  }

  # Actualizar el perdedor con los atributos del ganador
  hijo_perdedor$features <- hijo_ganador$features
  hijo_perdedor$plano_coord <- hijo_ganador$plano_coord
  hijo_perdedor$vectors[[1]] <- hijo_ganador$vectors[[1]]
  hijo_perdedor$vectors[[2]] <- hijo_ganador$vectors[[2]]

  hijo_perdedor$num <- num
  hijo_ganador$num <- num + 1

  return(list(hijo_perdedor, hijo_ganador))
}



#' @title Create New Population
#'
#' @description
#'   This function creates a new population through crossover.
#'
#' @param poblacion Population class object
#'
#' @return Updated Population class object with the new population generated
nueva_poblacion <- function(poblacion) {
  cat(sprintf("Creando nueva población por cruce desde %d hasta %d....\n", poblacion$tam_pob, 2 * poblacion$tam_pob))
  i <- poblacion$tam_pob

  while (length(poblacion$lista_soluciones) < 2 * poblacion$tam_pob) {
    padre <- torneo_elegir_padre(poblacion)
    madre <- torneo_elegir_padre(poblacion)
    while (padre$num == madre$num) {
      madre <- torneo_elegir_padre(poblacion)
    }

    hijos <- cruzar_soluciones(poblacion, padre, madre, i)
    for (hijo in hijos) {
      if (length(poblacion$lista_soluciones) >= 2 * poblacion$tam_pob) {
        break  # Evitar agregar más hijos si se alcanza el tamaño máximo
      }

      if (runif(1) < poblacion$p_mutacion) {
        hijo <- mutar_solucion(hijo, poblacion)
      }

      hijo <- evaluar_solucion(hijo)
      while (!hijo$evaluacion_exitosa) {
        hijo <- generar_solucion_aleatoria(hijo)
        hijo <- evaluar_solucion(hijo)
      }

      while (comprobar_clones(poblacion, hijo) == 1) {
        hijo <- mutar_solucion(hijo, poblacion)
        hijo <- evaluar_solucion(hijo)
      }

      poblacion$lista_soluciones[[length(poblacion$lista_soluciones) + 1]] <- hijo
      i <- i + 1
    }
  }

  return(poblacion)
}



#' @title Mutate Solution
#'
#' @description
#'   This function mutates a solution based on mutation probabilities and coordinates.
#'
#' @param solucion Solution object to mutate
#' @param poblacion Population class object that contains mutation parameters
#'
#' @return Mutated solution object
mutar_solucion <- function(solucion, poblacion) {
  random <- runif(1)

  # Mutación de vectores
  if (random < poblacion$p_mut_ind) {
    solucion <- mutar_vectores(solucion)
  }

  # Mutación de características
  if (poblacion$num_dim > poblacion$num_features) {
    for (i in seq_along(solucion$features)) {
      aleatorio <- runif(1)
      if (aleatorio < poblacion$p_mut_fea) {
        features_posibles <- setdiff(poblacion$inputs, solucion$features)
        solucion$features[i] <- sample(features_posibles, 1)
        solucion$plano_coord[i] <- runif(1, -1, 1)
      }
    }
  }

  # Mutación de planos
  for (i in seq_along(solucion$plano_coord)) {
    aleatorio <- runif(1)
    if (aleatorio < poblacion$p_mut_coord) {
      if (poblacion$mut_coord == 0) {
        porcentaje <- runif(1, 0, 0.25)
        aleatorio2 <- runif(1)
        if (aleatorio2 < 0.5) {  # incrementa
          solucion$plano_coord[i] <- solucion$plano_coord[i] * (1 + porcentaje)
          if (solucion$plano_coord[i] > 1) {
            solucion$plano_coord[i] <- solucion$plano_coord[i] / 10
          }
        } else {  # decrementa
          solucion$plano_coord[i] <- solucion$plano_coord[i] * (1 - porcentaje)
          if (solucion$plano_coord[i] < -1) {
            solucion$plano_coord[i] <- solucion$plano_coord[i] / 10
          }
        }
      } else {
        solucion$plano_coord[i] <- runif(1, -1, 1)
      }
    }
  }

  return(solucion)
}



#' @title Reduce Population Size
#'
#' @description
#'   This function reduces the population size based on fronts and crowding distance.
#'
#' @param poblacion Population class object to reduce
#'
#' @return Reduced Population class object
reducir_poblacion <- function(poblacion) {
  cat("Reduciendo de tamaño 2N a tamaño N....\n")

  # Crear nueva población vacía
  nueva_poblacion <- list(
    data = poblacion$data,
    costes = poblacion$costes,
    tam_pob = poblacion$tam_pob,
    inputs = poblacion$inputs,
    output = poblacion$output,
    num_features = poblacion$num_features,
    lista_soluciones = list()
  )

  # Contador de frentes
  cont_front <- 1

  # Añadir soluciones de los frentes a la nueva población hasta completar el tamaño
  while (length(nueva_poblacion$lista_soluciones) + length(poblacion$fronts[[cont_front]]$soluciones) <= poblacion$tam_pob) {
    soluciones_front <- lapply(poblacion$lista_soluciones, function(sol) {
      if (sol$front == cont_front) return(sol) else return(NULL)
    })
    soluciones_front <- soluciones_front[!sapply(soluciones_front, is.null)]
    nueva_poblacion$lista_soluciones <- c(nueva_poblacion$lista_soluciones, soluciones_front)
    cat(sprintf("FRONT %d: %d soluciones\n", cont_front, length(nueva_poblacion$lista_soluciones)))
    cont_front <- cont_front + 1
  }

  # Si no se ha completado la población, seleccionar las mejores soluciones restantes usando crowding distance
  if (length(nueva_poblacion$lista_soluciones) < poblacion$tam_pob && cont_front <= length(poblacion$fronts)) {
    soluciones_front <- crowding_distance(poblacion, cont_front)
    faltantes <- poblacion$tam_pob - length(nueva_poblacion$lista_soluciones)
    nueva_poblacion$lista_soluciones <- c(nueva_poblacion$lista_soluciones, soluciones_front[1:faltantes])
  }

  # Renombrar las soluciones
  for (i in seq_along(nueva_poblacion$lista_soluciones)) {
    nueva_poblacion$lista_soluciones[[i]]$num <- i - 1
  }

  return(nueva_poblacion)
}



#' @title Crowding Distance Calculation
#'
#' @description
#'   This function calculates the crowding distance for solutions in a specified front.
#'
#' @param poblacion Population class object
#' @param num_front The front number to calculate the crowding distance for
#'
#' @return List of solutions with updated crowding distances
crowding_distance <- function(poblacion, num_front) {
  # Extraer soluciones del frente especificado
  soluciones_front <- lapply(poblacion$lista_soluciones, function(sol) {
    if (sol$front == num_front) return(sol) else return(NULL)
  })
  soluciones_front <- soluciones_front[!sapply(soluciones_front, is.null)]

  num_sol_front <- length(soluciones_front)

  # Inicializar la distancia de hacinamiento a cero para todas las soluciones
  for (sol in soluciones_front) {
    sol$crowding_distance <- 0
  }

  # Procesar cada objetivo si hay suficientes soluciones para comparar
  if (num_sol_front > 2) {
    for (obj in 1:poblacion$num_obj) {
      # Ordenar las soluciones basándose en el objetivo actual
      soluciones_front <- soluciones_front[order(sapply(soluciones_front, function(x) x$objetivo[[obj]]))]

      # Asignar la distancia de hacinamiento infinita a las soluciones extremas
      soluciones_front[[1]]$crowding_distance <- Inf
      soluciones_front[[num_sol_front]]$crowding_distance <- Inf

      max_obj <- soluciones_front[[num_sol_front]]$objetivo[[obj]]
      min_obj <- soluciones_front[[1]]$objetivo[[obj]]

      if (max_obj > min_obj) {
        for (j in 2:(num_sol_front - 1)) {
          incremento <- (soluciones_front[[j + 1]]$objetivo[[obj]] - soluciones_front[[j - 1]]$objetivo[[obj]]) / (max_obj - min_obj)
          soluciones_front[[j]]$crowding_distance <- soluciones_front[[j]]$crowding_distance + incremento
        }
      }
    }
  }

  # Ordenar las soluciones por la distancia de hacinamiento de mayor a menor
  soluciones_front <- soluciones_front[order(sapply(soluciones_front, function(x) x$crowding_distance), decreasing = TRUE)]

  return(soluciones_front)
}
