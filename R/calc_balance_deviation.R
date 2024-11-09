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
#'
#' library(moses)
#'
#' # Example 1
#' C <- matrix(c(10, 10, 1, 5, 1, 4, 5, 15, 0),
#'             nrow = 3, byrow = TRUE)
#'
#' X <- matrix(c(0, 0, 1, 1, 1, 0),
#'            nrow = 2, byrow = TRUE)
#'
#' calc_balance_deviation(X, C)
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
#' # Generate a random allocation matrix
#' X = matrix(0,
#'            nrow = length(delta),
#'            ncol = nrow(C))
#' for (i in 1:ncol(X)) X[sample.int(nrow(X), 1), i] <- 1
#'
#' calc_balance_deviation(X, C)
#' }
#'

calc_balance_deviation <- function(X, C){
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

  pitilde <- Mtilde / matrix(rowSums(Mtilde), nrow = nrow(M), ncol = ncol(M), byrow = FALSE)
  return(abs(matrix(pi, nrow = nrow(M), ncol = ncol(M), byrow = TRUE) - pitilde))


  # mkdot <- rowSums(Mtilde) pimkdot <- expand.grid(mkdot, pi) pimkdot <-
  # matrix(pimkdot[, 1] * pimkdot[, 2], nrow = nrow(Mtilde), byrow = FALSE)
  # npimkdot <- expand.grid(mkdot, sum(pi)-pi) npimkdot <- matrix(npimkdot[, 1]
  # * npimkdot[, 2], nrow = nrow(Mtilde), byrow = FALSE)
  #
  # mask <- (Mtilde < pimkdot)*1
  #
  # fkj  <- mask * (1 - (Mtilde + 1e-9) / (pimkdot + 1e-9)) + (1 - mask) *
  # (Mtilde - pimkdot + 1e-9) / (npimkdot + 1e-9)
  #
  # return(fkj)

}
