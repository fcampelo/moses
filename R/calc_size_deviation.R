#' Calculate size deviation for all splits given an allocation matrix
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
#'
#' @return A vector containing the normalized size deviation scores of each split.
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
#' delta <- c(.75, .25)
#'
#' X <- matrix(c(0, 0, 1, 1, 1, 0),
#'            nrow = 2, byrow = TRUE)
#'
#' calc_size_deviation(X, C, delta)
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
#' # Generate a random allocation matrix
#' X = matrix(0,
#'            nrow = length(delta),
#'            ncol = nrow(C))
#' for (i in 1:ncol(X)) X[sample.int(nrow(X), 1), i] <- 1
#'
#' calc_size_deviation(X, C, delta)
#' }
#'

calc_size_deviation <- function(X, C, delta){
  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(C),
                          is.numeric(C),
                          all(C >= 0),
                          is.matrix(X),
                          all(as.vector(X) %in% c(0, 1)),
                          is.vector(delta),
                          is.numeric(delta),
                          sum(delta) == 1,
                          nrow(C) == ncol(X),
                          nrow(X) == length(delta))

  # Force delta in decreasing order
  delta <- sort(delta, decreasing = TRUE)
  # =======================================================================

  M  <- X %*% C
  Mtilde <- M[order(rowSums(M), decreasing = TRUE), ]

  Delta_s <- abs(delta - rowSums(Mtilde) / sum(Mtilde + 1e-9))

  return(Delta_s)

}
