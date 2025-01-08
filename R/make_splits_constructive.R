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
#' length 3 and add to 1. All weights must be non-negative.
#' @param X0 Initial (partial) allocation. Must be a binary numerical matrix
#' containing only zeroes and ones. Must have `nrow(C)` columns and `length(delta)`
#' rows. Each position \eqn{x_{ki}}
#' indicates the allocation (or not) of the i-th group to the k-th split.
#' Each group can only be allocated to (at most) one split (i.e.,
#' `colSums(X0)` \eqn{\leq 1}). The constructive heuristic will not change the
#' pre-allocations provided in X0, it will only allocate the remaining groups
#' to the splits.
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
#' X <- make_splits_constructive(C, delta, w)
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

make_splits_constructive <- function(C, delta, w = c(.5, .4, .1), X0 = NULL, rho = 1e-4){

  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(C),
                          is.numeric(C),
                          all(C >= 0),
                          is.vector(delta),
                          is.numeric(delta),
                          sum(delta) == 1,
                          is.numeric(w), length(w) == 3,
                          all(w >= 0), sum(w) == 1,
                          is.null(X0) | is.matrix(X0),
                          is.numeric(rho), length(rho) == 1, rho > 0)


  if(!is.null(X0)){
    assertthat::assert_that(nrow(X0) == length(delta),
                            ncol(X0) == nrow(C),
                            all(X0 %in% c(0, 1)))
    X <- X0
  } else {
    X <- matrix(0, nrow = length(delta), ncol = nrow(C))
  }

  # Force delta in decreasing order
  delta <- sort(delta, decreasing = TRUE)

  # Number of possible allocations:
  # nstates <- length(delta) ^ nrow(C)
  # =======================================================================


  # index from largest to smallest group
  idx <- order(rowSums(C), decreasing = TRUE)

  # Put pre-allocated groups at the top of the list
  iii <- which(colSums(X) == 1)

  if(length(iii) == nrow(C)){
    warning("X0 already provides a full allocation of groups. Skipping constructive_heuristic()")
    return(X0)
  }

  if(length(iii) > 0) idx <- c(iii, idx[-which(idx %in% iii)])

  for (j in (1+length(iii)):length(idx)){
    #i    <- idx[j]
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

    # Allocation criteria:

    # Smallest value of primary objective function
    lb <- which(f[, 1] == min(f[, 1]))

    # Tie breaker 1: allocate to empty splits
    if(length(lb) > 1) {
      ff <- rowSums(Xtmp)[lb]
      if(any(ff == 0)) lb <- lb[which(ff == 0)]
    }

    # Tie breaker 2: smallest value of secondary objective function
    if(length(lb) > 1) {
      ff <- f[lb, ]
      lb <- lb[which(ff[, 1] == min(ff[, 1]) & ff[, 2] == min(ff[, 2]))]
    }

    # Tie breaker 3: random allocation
    if(length(lb) > 1) lb <- sample(lb, 1)

    X[, icum] <- trials[[lb]]
  }

  # Return allocation matrix in decreasing order of split size
  M <- X %*% C
  idx <- order(rowSums(M), decreasing = TRUE)
  X <- X[idx, ]

  return(X)
}
