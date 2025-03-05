#' @title Definition of the Solucion Class
#'
#' @description
#'   Defines the `Solucion` class, which represents an individual solution.
#'
#' @param num Solution number
#' @param data Data set
#' @param costes Cost vector
#' @param inputs Names of input variables
#' @param output Name of the output variable
#' @param num_features Number of features
#'
#' @return Object of the `Solucion` class
#'
#' @export
Solucion <- function(num, data, costes, inputs, output, num_features) {
  # Definir variables de instancia
  solucion <- list(
    num = num,
    data = data,
    costes = costes,
    inputs = inputs,
    output = output,
    num_dim = length(inputs),
    num_features = num_features,
    data_sol = NULL,
    features = list(),
    vectors = list(),
    plano_coord = list(),
    plano_termino_b = list(),
    objetivo = numeric(2),  # Index 1 para distancia, Index 2 para epsilon
    sol_dom_por = 0,  # Inicializar contador de dominación
    lista_domina_a = list(),  # Lista de soluciones que esta solución domina
    lista_dominado_por = list(),  # Lista de soluciones que dominan a esta solución
    front = -1,
    evaluacion_exitosa = FALSE,
    crowding_distance = 0
  )
  return(solucion)
}



#' @title Method to convert the solution to a dictionary
#'
#' @description
#'   Method to convert the solution to a dictionary
#'
#' @param solucion Solucion class object
#'
#' @return List with the information of the solution
to_dict <- function(solucion) {
  return(list(
    #SOL = solucion$num,
    VECTORS = solucion$vectors,
    #PLANO_COOR = solucion$plano_coord,
    #PLANO_COOR_B = solucion$plano_termino_b,
    FEATURES = solucion$features,
    DIST = solucion$objetivo[1],  # Distancia
    EPS = solucion$objetivo[2],   # Epsilon
    # COSTE = solucion$obj_coste,  # Descomentar si obj_coste es utilizado
    # MC+ = solucion$obj_mal_clasificados$mc_pos,  # Descomentar si se utiliza
    # MC- = solucion$obj_mal_clasificados$mc_neg,  # Descomentar si se utiliza
    DOMINA_A = solucion$lista_domina_a,
    DOMINADO_POR = solucion$lista_dominado_por,
    SOL_DOM_POR = solucion$sol_dom_por,
    FRONT = solucion$front,
    CROW_DIST = solucion$crowding_distance
  ))
}



#' @title Method to Obtain a Random Class Vector
#'
#' @description
#'   Method to obtain a random class vector from the dataset.
#'
#' @param solucion Solucion class object containing the dataset
#' @param clase Class for which to obtain a vector
#'
#' @return Index of the random class vector
obtener_vector_clase <- function(solucion, clase) {
  indices <- which(solucion$data[[solucion$output]] == clase)
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
#' @param solucion Solucion class object to modify
#'
#' @return Solucion class object with a randomly generated solution
generar_solucion_aleatoria <- function(solucion) {
  # Seleccionar dos puntos de la clase A y B
  clases <- sort(unique(solucion$data[[solucion$output]]))
  solucion$vectors <- vector("list", length(clases))

  for (i in seq_along(clases)) {
    solucion$vectors[[i]] <- obtener_vector_clase(solucion, clases[i])
  }

  # Seleccionar aleatoriamente p características y p coordenadas entre [-1, 1]
  if (solucion$num_features == solucion$num_dim) {
    solucion$features <- solucion$inputs
  } else {
    fea_posibles <- solucion$inputs[which(solucion$inputs %in% names(solucion$data))] # Se verifica que las entradas existan como columnas
    solucion$features <- sample(fea_posibles, solucion$num_features)
  }

  # Calcular las coordenadas del plano
  solucion$plano_coord <- runif(solucion$num_features, -1, 1) # Coordenadas aleatorias para cada característica

  # Ajustar dataframe de solución
  indices <- unlist(solucion$vectors)
  if (length(indices) == 0 || length(solucion$features) == 0) {
    return(solucion) # Retorna la solución sin modificar si no hay índices válidos o características seleccionadas
  }

  solucion$data_sol <- solucion$data[indices, solucion$features, drop = FALSE]

  return(solucion)
}



#' @title Method to Construct Planes
#'
#' @description
#'   Method to construct planes by calculating the coefficients based on features.
#'
#' @param solucion Solucion class object containing data and features to construct planes
#'
#' @return Solucion class object with constructed planes
construir_planos <- function(solucion) {
  # Preparar el dataframe data_sol seleccionando solo las filas y columnas necesarias
  indices_validos <- unlist(solucion$vectors)
  indices_validos <- indices_validos[indices_validos <= nrow(solucion$data) & indices_validos > 0]

  if (length(indices_validos) == 0) {
    stop("No hay índices válidos en 'vectors'.")
  }

  # Filtrar las características para asegurarse de que están presentes en los datos
  features_validos <- solucion$features[solucion$features %in% names(solucion$data)]
  if (length(features_validos) == 0) {
    stop("No hay características válidas para la selección.")
  }

  # Subseleccionar el dataframe
  solucion$data_sol <- solucion$data[indices_validos, features_validos, drop = FALSE]

  # Inicializar el vector de términos del plano
  solucion$plano_termino_b <- numeric()

  # Calcular los términos independientes de la ecuación del plano
  for (i in seq_len(nrow(solucion$data_sol))) {
    indep <- 0
    for (j in seq_along(features_validos)) {
      indep <- indep + solucion$data_sol[i, j] * solucion$plano_coord[j]
    }
    # Añadir el valor negativo de 'indep' al vector de términos del plano
    solucion$plano_termino_b <- c(solucion$plano_termino_b, -indep)
  }

  # Añadir el término del plano intermedio
  if(length(solucion$plano_termino_b) > 0) {
    medio <- mean(solucion$plano_termino_b)
    solucion$plano_termino_b <- c(solucion$plano_termino_b, medio)
  } else {
    stop("No se han calculado términos del plano, por lo que no se puede calcular el medio.")
  }

  return(solucion)
}



#' @title Method to Calculate the Distance Objective
#'
#' @description
#'   This method calculates the distance objective based on the plane equation.
#'
#' @param solucion Solucion class object containing the plane equation
#'
#' @return Solucion class object with the distance objective calculated
calcular_objetivo_distancia <- function(solucion) {
  # Inicializar el objetivo de distancia a cero
  solucion$objetivo[[1]] <- 0
  denominador <- sum(solucion$plano_coord^2) ^ 0.5
  distancia <- abs(solucion$plano_termino_b[[1]] - solucion$plano_termino_b[[2]])

  # Calcular la distancia normalizada si el denominador no es cero
  if (denominador != 0) {
    solucion$objetivo[[1]] <- distancia / denominador
  } else {
    solucion$objetivo[[1]] <- -1  # Valor de error o indicativo de no calculable
  }
  return(solucion)
}



#' @title Method to Calculate the Epsilon Objective
#'
#' @description
#'   This method calculates the epsilon objective by evaluating misclassified points.
#'
#' @param solucion Solucion class object containing data and features for classification
#'
#' @return Solucion class object with the epsilon objective calculated
calcular_objetivo_epsilon <- function(solucion) {
  # Inicializar contadores de clasificaciones malas
  mc_pos <- 0
  mc_neg <- 0
  solucion$objetivo[[2]] <- 0

  # Iterar sobre todos los datos para calcular la suma de errores
  for (i in 1:nrow(solucion$data)) {
    clase <- solucion$data[i, solucion$output]
    distancia <- 0

    # Calcular la distancia lineal usando los coeficientes
    for (indice_coor in seq_along(solucion$features)) {
      fea <- solucion$features[[indice_coor]]
      distancia <- distancia + solucion$plano_coord[indice_coor] * solucion$data[i, fea]
    }

    # Ajustar por el término independiente y contar errores
    if (clase == -1) {
      distancia <- distancia + solucion$plano_termino_b[[1]]
      if (distancia < 0) {
        mc_neg <- mc_neg + 1
        solucion$objetivo[[2]] <- solucion$objetivo[[2]] + abs(distancia)
      }
    } else {
      distancia <- distancia + solucion$plano_termino_b[[2]]
      if (distancia > 0) {
        mc_pos <- mc_pos + 1
        solucion$objetivo[[2]] <- solucion$objetivo[[2]] + abs(distancia)
      }
    }
  }

  # Opcionalmente, podría almacenar los conteos de clasificación incorrecta si fuera necesario
  solucion$mc_pos <- mc_pos
  solucion$mc_neg <- mc_neg

  return(solucion)
}



# calcular_objetivo_costes <- function(solucion) {
#   solucion$obj_coste <- 0
#   for (feature in solucion$features) {
#     solucion$obj_coste <- solucion$obj_coste + solucion$costes[feature]
#   }
#   return(solucion)
# }



#' @title Evaluate Solution
#'
#' @description
#'   This function evaluates a solution by calculating its objectives.
#'
#' @param solucion Solucion class object to be evaluated
#'
#' @return Solucion class object with updated evaluation status and objectives
evaluar_solucion <- function(solucion) {
  # Verifica si la suma de plano_coord es cero
  if (sum(solucion$plano_coord) == 0) {
    cat("No se puede evaluar la solucion, no hay coordenadas.\n")
    solucion$evaluacion_exitosa <- FALSE
    return(solucion)
  } else {
    # Procedimientos para evaluar la solución si hay coordenadas válidas
    solucion <- construir_planos(solucion)
    solucion <- calcular_objetivo_distancia(solucion)
    solucion <- calcular_objetivo_epsilon(solucion)
    # solucion <- calcular_objetivo_costes(solucion) si es necesario

    # Asumir que las evaluaciones fueron exitosas si se alcanza este punto
    solucion$evaluacion_exitosa <- TRUE
    return(solucion)
  }
}



#' @title Method to Determine Dominance
#'
#' @description
#'   This function determines if one solution dominates another based on objectives.
#'
#' @param solucion1 First solution object
#' @param solucion2 Second solution object
#'
#' @return TRUE if solucion1 dominates solucion2, otherwise FALSE
dominar <- function(solucion1, solucion2) {
  # Comprobar si alguno de los objetivos es NA o no numérico antes de comparar
  if (is.na(solucion1$objetivo[[1]]) || !is.numeric(solucion1$objetivo[[1]]) ||
      is.na(solucion2$objetivo[[1]]) || !is.numeric(solucion2$objetivo[[1]]) ||
      is.na(solucion1$objetivo[[2]]) || !is.numeric(solucion1$objetivo[[2]]) ||
      is.na(solucion2$objetivo[[2]]) || !is.numeric(solucion2$objetivo[[2]])) {
    # Retorna FALSE o maneja el caso de forma que se ajuste a tu lógica de aplicación
    return(FALSE)
  }

  # Procede con la comparación si todos los valores son adecuados
  return(solucion1$objetivo[[1]] >= solucion2$objetivo[[1]] &&
           solucion1$objetivo[[2]] <= solucion2$objetivo[[2]])
}



#' @title Method to Determine Dominance with Three States
#'
#' @description
#'   This function determines the dominance state between two solutions.
#'
#' @param solucion1 First solution object
#' @param solucion2 Second solution object
#'
#' @return 0 if the solutions are equal, 1 if solucion1 dominates solucion2, 2 if solucion2 dominates solucion1
dominar2 <- function(solucion1, solucion2) {
  # Comprobar si los objetivos son iguales
  if (solucion1$objetivo[[1]] == solucion2$objetivo[[1]] && solucion1$objetivo[[2]] == solucion2$objetivo[[2]]) {
    return(0)  # No hay dominancia, son iguales
  } else if (solucion1$objetivo[[1]] >= solucion2$objetivo[[1]] && solucion1$objetivo[[2]] <= solucion2$objetivo[[2]]) {
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
#' @param solucion1 First solution object
#' @param solucion2 Second solution object
#'
#' @return The solution object with better dominance or front
comparar_soluciones <- function(solucion1, solucion2) {
  if (solucion1$front == solucion2$front) {
    if (length(solucion1$lista_dominado_por) == length(solucion2$lista_dominado_por)) {
      moneda <- sample(0:1, 1)
      if (moneda == 0) {
        return(solucion1)
      } else {
        return(solucion2)
      }
    } else if (length(solucion1$lista_dominado_por) < length(solucion2$lista_dominado_por)) {
      return(solucion1)
    } else {
      return(solucion2)
    }
  } else {
    if (solucion1$front > solucion2$front) {
      return(solucion2)
    } else {
      return(solucion1)
    }
  }
}



#' @title Method to Mutate Vectors
#'
#' @description
#'   This function mutates the vectors in the solution by assigning new class vectors.
#'
#' @param solucion Solucion class object to mutate
#'
#' @return Solucion class object with updated vectors
mutar_vectores <- function(solucion) {
  # Obtener las clases únicas y ordenarlas
  clases <- sort(unique(solucion$data[[solucion$output]]))

  # Iterar sobre cada clase
  for (i in seq_along(clases)) {
    # Llamar a obtener_vector_clase para la clase actual
    vector <- obtener_vector_clase(solucion, clases[i])

    # Verificar que el nuevo vector sea diferente al actual y si no, buscar otro
    while (identical(vector, solucion$vectors[[i]])) {
      vector <- obtener_vector_clase(solucion, clases[i])
    }

    # Actualizar el vector en la solución
    solucion$vectors[[i]] <- vector
  }

  return(solucion)
}
