#' Calculate split homogeneity for all splits given an allocation matrix.
#'
#' @param X matrix of binary allocation variables. Each position `x_{ki}`
#' indicates the allocation or not of the i-th group to the k-th split. Note
#' that each group must be allocated to a single split, i.e.,
#' colSums(X) must be a vector of ones.
#'
#' @return @return A vector containing the normalized homogeneity scores of
#' each split.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(moses)
#'
#' # Build example with 3 splits, 6 groups
#'
#' X = matrix(0, nrow = 3, ncol = 6)
#' X[1, 1] <- 1
#' X[2, 2:3] <- 1
#' X[3, 4:6] <- 1
#'
#' calculate_split_homogeneity(X)
#' }
#'

calculate_split_homogeneity <- function(X, C){
  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(X),
                          all(as.vector(X) %in% c(0, 1)))

  # =======================================================================

  Nk <- rowSums(X)
  return(1 / pmax(1, Nk))
}
