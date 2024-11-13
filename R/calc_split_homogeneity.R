#' Calculate split homogeneity for all splits given an allocation matrix.
#'
#' @param X matrix of binary allocation variables. Each position \eqn{x_{ki}}
#' indicates the allocation or not of the i-th group to the k-th split. Note
#' that each group must be allocated to a single split, i.e.,
#' colSums(X) must be a vector of ones.
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position \eqn{c_{ij}} contains
#' the number of examples of the j-th Class contained in the i-th group.
#'
#' @return A vector containing the normalized homogeneity scores of
#' each split.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#'
#' library(moses)
#'
#' C <- matrix(c(10, 10, 1, 5, 1, 4, 5, 15, 0),
#'             nrow = 3, byrow = TRUE)
#'
#' X <- matrix(c(0, 0, 1, 1, 1, 0),
#'            nrow = 2, byrow = TRUE)
#'
#' calc_split_homogeneity(X, C)
#'

calc_split_homogeneity <- function(X, C){
  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(X),
                          all(as.vector(X) %in% c(0, 1)))

  # =======================================================================

  Nk <- rowSums(X)
  G <- nrow(C)
  ifelse(Nk == 0, 0, 1 - Nk / G)
}
