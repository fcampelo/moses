calc_objective_functions <- function(X, C, delta){

  Dsize <- calc_size_deviation(X, C, delta)
  Dbal  <- calc_balance_deviation(X, C)
  Dhomo <- calc_split_homogeneity(X, C)

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
