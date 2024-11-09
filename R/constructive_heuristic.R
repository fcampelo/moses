#' Generate an initial group allocation solution based on a greedy heuristic.
#'
#' Generates an initial allocation matrix.
#'
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position `c_{ij}` contains
#' the number of examples of the j-th Class contained in the i-th group.
#' @param delta vector of desired split proportions (must add up to one). It
#' is useful (but not mandatory) to use a vector delta that is sorted
#' in decreasing order, to prevent later confusion.
#' @param w vector of weights for function aggregation. Must be of
#' length 3 and add to 1. All weights must be non-negative. If a
#' scalar value is passed then the weights are taken as equal for all
#' objectives. Weights are scaled proportionally if they don't add to 1.
#' @param rho small positive value for augmented Tchebycheff scalarisation.
#' @param cl a cluster object, created by package `parallel` or package `snow`.
#' If `NULL` the function will use the registered default cluster (or none).
#'
#' @return Matrix of binary allocation variables. Each position `x_{ki}`
#' indicates the allocation or not of the i-th group to the k-th split.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'

constructive_heuristic <- function(C, delta, w, rho = 1e-4, cl = NULL){

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

}

# constructive_heuristic_old <- function(class_counts,
#                                    target_proportions,
#                                    weights = 1,
#                                    rho     = 1e-2,
#                                    class_balance = NULL){
#
#   class_counts <- class_counts[order(class_counts$Size,
#                                      decreasing = TRUE), ]
#   class_counts$Split <- NA
#
#   alloc <- as.data.frame(matrix(data = NA,
#                                 nrow = nrow(class_counts),
#                                 ncol = length(target_proportions),
#                                 byrow = TRUE))
#
#   for (i in seq_along( class_counts$Cluster)){
#     alloc[i, ] <- 1:length(target_proportions)
#     f0 <- lapply(alloc,
#                  scalar_objective_function,
#                  class_counts = class_counts,
#                  target_proportions = target_proportions,
#                  weights = weights,
#                  rho     = rho,
#                  class_balance = class_balance)
#
#     # Extract best allocation
#     idx <- which(sapply(f0, function(x)x$f.AT) == min(sapply(f0, function(x)x$f.AT)))
#
#     if (length(idx) > 1){
#       idx2 <- which(sapply(f0[idx], function(x)x$f.AT2) == min(sapply(f0[idx], function(x)x$f.AT2)))
#
#       if (length(idx2) > 1) idx2 <- sample(idx2, 1)
#
#       idx <- idx[idx2]
#     }
#
#     alloc[i, ] <- idx
#   }
# }
