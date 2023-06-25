#' @param alloc_vector vector containing the allocation information for
#' each cluster of observations. Must be a vector of
#' length `nrow(class_counts)`. Each position __i__ should contain
#' an integer between 1 and __nsplits__ indicating the split number to
#' which the __i__-th cluster is allocated.
#' @param class_counts data frame with one row per cluster.
#' Must contain a column Cluster` (with the cluster number/ID) plus
#' additional columns with counts of occurrences for each unique class
#' and columns with the class proportions in each cluster, plus a
#' column with the total cluster size.
#' @param obj_functions character vector with the names of the objective
#' functions to include.
#' @param weights vector of weights for function aggregation. Must be of
#' the same length as `obj_functions`. Weights must be non-negative. If a
#' scalar value is passed then the weights are taken as equal for all
#' objectives.
#' @param class_balance A vector with the proportions of each class in
#' the data. Calculated internally if `NULL`.
#'

scalar_objective_function <- function(alloc_vector,
                                      class_counts,
                                      target_proportions,
                                      obj_functions,
                                      weights = 1,
                                      class_balance = NULL){

  # =======================================================================
  # Sanity checks and initial definitions


  idx <- grep("Class.", names(class_counts))
  if(is.null(class_balance)){
    class_balance <- colSums(class_counts[, idx]) / sum(class_counts[, idx])
  }

  if(!("Size" %in% names(class_counts))){
    class_counts$Size <- rowSums(class_counts[, grep("Class\\.", names(class_counts))])
  }

  if(length(weights == 1)) {
    weights <- rep(1 / length(obj_functions),
                   times = length((obj_functions)))
  }

  weights <- weights / sum(weights)

  target_proportions <- target_proportions / sum(target_proportions)



  # =======================================================================

  f <- numeric(length(obj_functions))

  for (i in seq_along(f)){
    f[i] <- do.call(obj_functions[i],
                    args = list(alloc_vector = alloc_vector,
                                class_counts = class_counts,
                                class_balance = class_balance,
                                target_proportions = target_proportions))
  }

  y = sum(weights * f) # <--- change to Aug Chebyshev


  return(list(f.vector = f,
              f.scalar = y))

}
