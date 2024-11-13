#' Calculate the scalarised objective function for data splitting
#'
#' This function calculates the objective function for the data splitting
#' optimisation. See section **Scalarisation details** for a full description.
#'
#' @section Scalarisation details:
#' This function calculates a scalarised objective function that is minimised
#' as part of the definition of the optimal data splits for data
#' under a clustering/grouping structure, with potential data imbalance and
#' with an arbitrary number of classes.
#' The three objective functions to be minimised are:
#' \itemize{
#'    \item f1 (*Size deviation*), which measures how much the resulting splits
#'    (given a tentative allocation) deviate from the desired proportion.
#'    \item f2 (*Balance deviation*), which measures how much the class balance of
#'    the splits deviates from the overall class balance of the full dataset.
#'    \item f3 (*Homogeneity score*), which quantifies how homogeneous (using the
#'    number of different clusters/groups in each split as a measure of
#'    diversity) the splits are.
#' }
#'
#'
#' These three objective functions are aggregated using the augmented
#' Tchebycheff scalarisation function,
#'
#'\deqn{
#' f_{AT} = \max_\ell \left[ w_\ell f_\ell(x)\right ]  + \rho\sum_\ell w_\ell f_\ell(x)
#' }
#'
#' where \eqn{w} represent the weights provided by the user, and \eqn{\rho} is a small
#' positive value.
#'
#' Depending on the application, and particularly on the relative size of different
#' groups, users may find it useful to reduce or even zero the value of w[3], i.e.,
#' disregard the homogeneity objective, particularly if the number of splits is relatively
#' large (e.g., when using this function to define cross-validation folds).
#'
#' @param X matrix of binary allocation variables. Each position `x_{ki}`
#' indicates the allocation or not of the i-th group to the k-th split. Note
#' that each group must be allocated to a single split, i.e.,
#' colSums(X) must be a vector of ones.
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position `c_{ij}` contains
#' the number of examples of the j-th Class contained in the i-th group.
#' @param delta vector of desired split proportions (must add up to one). Splits
#' are always created in decreasing order of (desired) size, so it is useful
#' (but not mandatory) to pass a vector delta that is already sorted in
#' decreasing order, to prevent later confusion.
#' @param w vector of weights for function aggregation. Must be of
#' length 3 and add to 1. All weights must be non-negative. If a
#' scalar value is passed then the weights are taken as equal for all
#' objectives. Weights are scaled proportionally if they don't add to 1.
#' @param rho small positive value for augmented Tchebycheff scalarisation.
#' @param which version of the objectives should be evaluated? Accepts
#' "primary", "secondary", "both" (vector with primary and secondary values), or
#' "raw" (matrix of all individual primary and secondary objective functions).
#'
#' @seealso `calc_size_deviation`, `calc_balance_deviation`, `calc_split_homogeneity`
#'
#' @export
#'
#' @examples
#'
#' library(moses)
#'
#' # Example 1
#' C <- matrix(c(10, 10, 1, 5, 1, 4, 5, 15, 0),
#'             nrow = 3, byrow = TRUE)
#'
#' delta <- c(.75, .25)
#'
#' X <- matrix(c(0, 0, 1, 1, 1, 0),
#'            nrow = 2, byrow = TRUE)
#'
#' w = c(.4, .4, .2)
#'
#' calc_scalarised_objective(X, C, delta, w)
#'
#' \dontrun{
#' # Example 2
#' fpath1 <- system.file("diamond", "bfv_proteins.fa", package="moses")
#' fpath2 <- system.file("diamond", "bfv_peptides.rds", package="moses")
#'
#' # Calculate clusters using sequence data
#' mycl <- extract_clusters_cdhit(seqfile = fpath1, diss_threshold = 0.2)
#'
#' # Load data frame with classes
#' df <- readRDS(fpath2)
#'
#' # Consolidate class counts
#' C <- consolidate_class_counts(mycl$clusters, df)
#'
#' # Desired allocation proportions
#' delta <- c(.6, .2, .2)
#'
#' w = c(.5, .4, .1)
#'
#' # Generate a random allocation matrix
#' X = matrix(0,
#'            nrow = length(delta),
#'            ncol = nrow(C))
#' for (i in 1:ncol(X)) X[sample.int(nrow(X), 1), i] <- 1
#'
#' calc_scalarised_objective(X, C, delta, w, which = "both")
#' }

calc_scalarised_objective <- function(X, C, delta, w, rho = 1e-4, which = "primary"){

  Fvals <- calc_objective_functions(X, C, delta)
  pri   <- max(w * Fvals$Primary) + rho * sum(w * Fvals$Primary)
  sec   <- max(w * Fvals$Secondary) + rho * sum(w * Fvals$Secondary)

  if (which == "primary") toRet <- pri
  if (which == "secondary") toRet <- sec
  if (which == "both") toRet <- c(pri, sec)
  if (which == "raw") toRet <- Fvals

  return(toRet)
}
