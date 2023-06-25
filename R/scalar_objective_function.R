#' @param alloc_vector vector containing the allocation information for
#' each cluster of observations. Must be a vector of
#' length `nrow(class_counts)`. Each position __i__ should contain
#' an integer between 1 and __nsplits__ indicating the split number to
#' which the __i__-th cluster is allocated.
#' @param class_counts data frame with one row per cluster.
#' Must contain a column Cluster` (with the cluster number/ID) plus
#' additional columns with counts of occurrences for each unique class
#' and columns with the class proportions in each cluster.
#' @param obj_functions character vector with the names of the objective
#' functions to include.
#' @param weights vector of weights for function aggregation. Must be of
#' the same length as `obj_functions`. Weights must be non-negative. If a
#' scalar value is passed then the weights are taken as equal for all
#' objectives.
#' @param cluster_sizes a vector with the total number of observations in
#'  each cluster. Calculated as `rowSums(class_counts[, -1])`.
#' @param class_balance A vector with the proportions of each class in
#' the data. Calculated as `colSums(class_counts[, -1]) / sum(class_counts[, -1])`.
#'

scalar_objective_function <- function(alloc_vector,
                                      class_counts,
                                      obj_functions,
                                      weights = 1,
                                      cluster_sizes = NULL,
                                      class_balance = NULL){

  # =======================================================================
  # Sanity checks and initial definitions



  if(is.null(class_balance)){
    class_balance <- colSums(class_counts[, -1]) / sum(class_counts[, -1])
  }

  if(is.null(cluster_sizes)){
    cluster_sizes <- rowSums(class_counts[, -1])
  }

  if(length(weights == 1)) {
    weights <- rep(1 / length(obj_functions),
                   times = length((obj_functions)))
  }

  weights <- weights / sum(weights)


  # =======================================================================

  f <- numeric(length(obj_functions))

  for (i in seq_along(f)){
    f[i] <- do.call(obj_functions[i],
                    args = list(alloc_vector = alloc_vector,
                                class_counts = class_counts,
                                cluster_sizes = cluster_sizes,
                                class_balance = class_balance))
  }

  y = sum(weights * f) # <--- change to Aug Chebyshev


  return(list(f.vector = f,
              f.scalar = y))

}
