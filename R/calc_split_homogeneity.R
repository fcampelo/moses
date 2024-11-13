#' Calculate split homogeneity for all splits given an allocation matrix.
#'
#' @param X matrix of binary allocation variables. Each position `x_{ki}`
#' indicates the allocation or not of the i-th group to the k-th split. Note
#' that each group must be allocated to a single split, i.e.,
#' colSums(X) must be a vector of ones.
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position `c_{ij}` contains
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
#' X <- matrix(c(1, 0, 0, 0, 0, 0,
#'               0, 1, 1, 0, 0, 0,
#'               0, 0, 0, 1, 1, 1),
#'            nrow = 3, byrow = TRUE)
#'
#' calc_split_homogeneity(X)
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
