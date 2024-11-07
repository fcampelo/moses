#' @export
calculate_scalarised_objective <- function(X, C, delta, w, rho = 1e-4){

  Fvals <- calculate_objective_functions(X, C, delta)

  return(max(w * Fvals$Primary) + rho * sum(w * Fvals$Primary))
}
