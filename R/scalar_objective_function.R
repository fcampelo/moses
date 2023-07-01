#' Calculate the scalarised objective function for data splitting
#'
#' This function calculates the objective function for the data splitting
#' optimisation. See section **Scalarisation details** for a full description.
#'
#' @section Scalarisation details:
#' This function calculates a scalarised objective function that is minimised
#' as part of the definition of the optimal data splits for sequence data
#' under a clustering/grouping structure, with potential data imbalance and
#' with an arbitrary number of classes.
#' The three objective functions (all of which are set to be minimised, i.e.,
#' smaller = better) are:
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
#' These three objective functions are aggregated using the augmented
#' Tchebycheff scalarisation function,
#'
#' $$
#' f_{AT} = \max_\ell\left{w_\ell f_\ell(x)\right} + \rho\sum_{\ell}w_\ell f_\ell(x)
#' $$
#'
#' where $w$ represent the weights provided by the user, and $\rho$ is a small
#' positive value.
#'
#' @param alloc_vector vector containing the allocation information
#' for each cluster of observations. Must be a vector of
#' length `nrow(class_counts)`. Each position __i__ should contain
#' an integer between 1 and __nsplits__ indicating the split number to
#' which the __i__-th cluster is allocated.
#' @param class_counts data frame with one row per cluster.
#' Must contain a column Cluster` (with the cluster number/ID) plus
#' additional columns with counts of occurrences for each unique class
#' and columns with the class proportions in each cluster, plus a
#' column with the total cluster size.
#' @param weights vector of weights for function aggregation. Must be of
#' length 3 and add to 1. All weights must be non-negative. If a
#' scalar value is passed then the weights are taken as equal for all
#' objectives. Weights are scaled proportionally if they don't add to 1.
#' @param rho small positive value for augmented Tchebycheff scalarisation.
#' @param class_balance A data frame with the proportions of each class in
#' the data. Calculated internally if `NULL`.
#'

scalar_objective_function <- function(alloc_vector,
                                      class_counts,
                                      target_proportions,
                                      weights = 1,
                                      rho     = 1e-6,
                                      class_balance = NULL){

  # =======================================================================
  # Sanity checks and initial definitions


  idx <- grep("Class.", names(class_counts))
  if(is.null(class_balance)){
    class_balance <- as.data.frame(t(colSums(class_counts[, idx]) / sum(class_counts[, idx])))
  }

  if(!("Size" %in% names(class_counts))){
    class_counts$Size <- rowSums(class_counts[, grep("Class\\.", names(class_counts))])
  }

  weights <- weights / sum(weights)

  target_proportions <- target_proportions / sum(target_proportions)



  # =======================================================================

  splits <- class_counts %>%
    dplyr::mutate(Split  = alloc_vector) %>%
    dplyr::group_by(Split) %>%
    dplyr::summarise(C.. = sum(class_counts$Size),
                     cardSk = sum(Size),
                     normSk = cardSk / C..,
                     Nk     = dplyr::n(),
                     dplyr::across(dplyr::starts_with("Class"),
                                   sum))

  CCount  <- dplyr::select(splits, dplyr::starts_with("Class"))
  CTarget <- splits$cardSk * cbind(splits[, 1], class_balance)[, -1]
  flag <- CCount < CTarget

  fj  <- flag * (1 - CCount/CTarget) +
    (!flag) * ((Ccount - CTarget) / (splits$cardSk - CTarget))


  Delta.Sk     <- abs(splits$normSk - sort(target_proportions))
  Delta.Bk.max <- apply(fj, 1, max)
  Delta.Bk.avg <- rowMeans(fj)
  Hk           <- ifelse(splits$Nk == 0, 1, 1 / splits$Nk)

  obj.funs <- data.frame(prime  = c(max(Delta.Sk),
                                    max(Delta.Bk.max),
                                    max(Hk)),
                         second = c(mean(Delta.Sk),
                                    mean(Delta.Bk.avg),
                                    median(Hk)))


  y1 <- max(weights * obj.funs$prime) + rho * sum((weights * obj.funs$prime))
  y2 <- max(weights * obj.funs$second) + rho * sum((weights * obj.funs$second))

  return(list(obj.funs = obj.funs,
              f.AT     = y1,
              f.AT2    = y2))

}
