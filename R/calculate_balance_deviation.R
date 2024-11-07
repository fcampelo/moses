#' Calculate balance deviation for all splits given an allocation matrix.
#'
#' @param X matrix of binary allocation variables. Each position `x_{ki}`
#' indicates the allocation or not of the i-th group to the k-th split. Note
#' that each group must be allocated to a single split, i.e.,
#' colSums(X) must be a vector of ones.
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position `c_{ij}` contains
#' the number of examples of the j-th Class contained in the i-th group.
#'
#' @return A matrix containing the normalized balance deviation scores of each
#' class on each split (splits as rows, classes as columns).
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
#' # Build example with 6 groups, 3 classes, 2 splits
#'
#' C <- matrix(c(1, 0, 0,
#'               3, 5, 0,
#'               4, 5, 6,
#'               4, 5, 6,
#'               4, 5, 6,
#'               4, 5, 6),
#'               ncol = 3, byrow = TRUE)
#'
#' X = matrix(0,
#'            nrow = 2,
#'            ncol = nrow(C))
#' X[1, 1:2] <- 1
#' X[2, 3:6] <- 1
#'
#' calculate_balance_deviation(X, C)
#' }
#'

calculate_balance_deviation <- function(X, C){
  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(C),
                          is.numeric(C),
                          all(C >= 0),
                          is.matrix(X),
                          all(as.vector(X) %in% c(0, 1)),
                          nrow(C) == ncol(X))

  # =======================================================================

  pi     <- colSums(C) / sum(C)
  M      <- X %*% C
  Mtilde <- M[order(rowSums(M), decreasing = TRUE), ]

  mkdot <- rowSums(Mtilde)

  pimkdot <- expand.grid(mkdot, pi)
  pimkdot <- matrix(pimkdot[, 1] * pimkdot[, 2],
                    nrow = nrow(Mtilde),
                    byrow = FALSE)
  npimkdot <- expand.grid(mkdot, 1-pi)
  npimkdot <- matrix(npimkdot[, 1] * npimkdot[, 2],
                     nrow = nrow(Mtilde),
                     byrow = FALSE)

  mask <- (Mtilde < pimkdot)*1

  fkj  <- mask * (1 - Mtilde / (pimkdot + 1e-9)) +
    (1 - mask) * (Mtilde - pimkdot) / (npimkdot + 1e-9)

  return(fkj)

}
