#' Generate an initial group allocation solution based on a greedy heuristic.
#'
#' Generates an initial allocation matrix.
#'
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position \eqn{c_{ij}} contains
#' the number of examples of the j-th Class contained in the i-th group.
#' @param delta vector of desired split proportions (must add up to one). It
#' is useful (but not mandatory) to use a vector delta that is sorted
#' in decreasing order, to prevent later confusion.
#' @param w vector of weights for function aggregation. Must be of
#' length 3 and add to 1. All weights must be non-negative. If a
#' scalar value is passed then the weights are taken as equal for all
#' objectives. Weights are scaled proportionally if they don't add to 1.
#' @param rho small positive value for augmented Tchebycheff scalarisation.
#'
#' @return Matrix of binary allocation variables. Each position \eqn{x_{ki}}
#' indicates the allocation or not of the i-th group to the k-th split.
#' **NOTE**: The allocation matrix is always returned with splits in decreasing
#' order of size.
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
#' \dontrun{
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
#' # Objective weights
#' w = c(.5, .4, .1)
#'
#' X <- constructive_heuristic(C, delta, w)
#'
#' # Check allocation:
#' M <- X %*% C
#'
#' data.frame(Desired.prop = sort(delta, decreasing = TRUE),
#'            Actual.prop  = rowSums(M) / sum(M),
#'            Groups.per.split = rowSums(X))
#'
#' # Balance of data in each split
#' Dbal <- cbind(data.frame(Overall = colSums(C) / sum(C)),
#'               apply(M, 1, function(z) z/sum(z)))
#' names(Dbal)[-1] <- paste0("Split.", names(Dbal)[-1])
#' Dbal
#' }


constructive_heuristic <- function(C, delta, w, rho = 1e-4){

  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(C),
                          is.numeric(C),
                          all(C >= 0),
                          is.vector(delta),
                          is.numeric(delta),
                          sum(delta) == 1)

  # Force delta in decreasing order
  delta <- sort(delta, decreasing = TRUE)

  # Number of possible allocations:
  nstates <- length(delta) ^ nrow(C)
  # =======================================================================

  X <- matrix(0, nrow = length(delta), ncol = nrow(C))

  # index from largest to smallest group
  idx <- order(rowSums(C), decreasing = TRUE)

  for (j in seq_along(idx)){
    i    <- idx[j]
    icum <- idx[1:j]
    Ctmp <- C[icum, , drop = FALSE]
    Xtmp <- X[, icum, drop = FALSE]

    trials <- lapply(1:nrow(Xtmp),
                     function(k, Y){Y[k, ncol(Y)] <- 1; Y},
                     Y = Xtmp)

    f <- lapply(trials,
                calc_scalarised_objective,
                C = Ctmp, delta = delta, w = w, rho = rho,
                which = "both") %>%
      do.call(what = rbind)

    lb <- which(f[, 1] == min(f[, 1]))
    if(length(lb) > 1) lb <- which(f[, 1] == min(f[, 1]) & f[, 2] == min(f[, 2]))
    if(length(lb) > 1) lb <- sample(lb, 1)

    X[, icum] <- trials[[lb]]
  }

  # Return allocation matrix in decreasing order of split size
  M <- X %*% C
  idx <- order(rowSums(M), decreasing = TRUE)
  X <- X[idx, ]

  return(X)
}
