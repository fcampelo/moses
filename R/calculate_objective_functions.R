#' @export
calculate_objective_functions <- function(X, C, delta){

  Dsize <- calculate_size_deviation(X, C, delta)
  Dbal  <- calculate_balance_deviation(X, C)
  Dhomo <- calculate_split_homogeneity(X)

  # Size deviation objective functions
  f1  <- max(Dsize)
  f1. <- mean(Dsize)

  # Balance deviation objective functions
  f2  <- max(apply(Dbal, 1, max))
  f2. <- mean(rowMeans(Dbal))

  # Split homogeneity objective functions
  f3  <- max(Dhomo)
  f3. <- mean(Dhomo)

  return(data.frame(Objective = c("Size", "Balance", "Homogeneity"),
                    Primary   = c(f1, f2, f3),
                    Secondary = c(f1., f2., f3.)))
}
