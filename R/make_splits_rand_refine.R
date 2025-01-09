#' Create a randomised data split.
#'
#' Generates a randomised allocation matrix that attempts to fulfill the
#' objectives of (i) adherence to the desired split proportions, and (ii)
#' homogeneous class balance across splits.
#'
#' @param C matrix of class counts per group. Commonly calculated using
#' [consolidate_class_counts()]. Each position \eqn{c_{ij}} contains
#' the number of examples of the j-th Class contained in the i-th group.
#' @param delta vector of desired split proportions (must add up to one). It
#' is useful (but not mandatory) to use a vector delta that is sorted
#' in decreasing order, to prevent later confusion.
#' @param w vector of weights for function aggregation. Must be of
#' length 2 and add to 1. All weights must be non-negative. `w[1]` represents the
#' relative importance of the split sizes and `w[2]` that of maintaining a
#' generally similar class balance in the splits.
#' @param maxiter integer; maximum number or refinement iterations.
#' @param X0 Initial (partial) allocation. Must be a binary numerical matrix
#' containing only zeroes and ones. Must have `nrow(C)` columns and `length(delta)`
#' rows. Each position \eqn{x_{ki}}
#' indicates the allocation (or not) of the i-th group to the k-th split.
#' Each group can only be allocated to (at most) one split (i.e.,
#' `colSums(X0)` \eqn{\leq 1}). Allocations provided in X0 are maintained in the
#' final solution returned.
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
#' # Desired allocation proportions (train-validation-test)
#' delta <- c(.6, .2, .2)
#'
#' X <- make_splits_rand_refine(C, delta)
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
#'
#' # Desired allocation proportions (5-fold CV folds + test)
#' delta <- c(rep(.16, 5), .2)
#' X <- make_splits_rand_refine(C, delta, w = c(.8, .2))
#' }


make_splits_rand_refine <- function(C, delta, X0 = NULL, w = c(2/3, 1/3), maxiter = 20){

  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.matrix(C),
                          is.numeric(C),
                          all(C >= 0),
                          is.vector(delta),
                          is.numeric(delta),
                          sum(delta) == 1,
                          is.null(X0) | is.matrix(X0))


  if(!is.null(X0)){
    assertthat::assert_that(nrow(X0) == length(delta),
                            ncol(X0) == nrow(C),
                            all(X0 %in% c(0, 1)))
    X <- X0
  } else {
    X <- matrix(0, nrow = length(delta), ncol = nrow(C))
  }

  vargroups <- which(colSums(X) == 0)

  # Force delta in decreasing order
  delta <- sort(delta, decreasing = TRUE)

  # =======================================================================

  # Initialise allocation list
  alloc <- apply(X, 2, \(x) {
    ii <- which(x == 1)
    if(length(ii) == 0) ii <- NA
    ii})

  # Complete with random allocation
  alloc[is.na(alloc)] <- sample(1:length(delta), length(vargroups), replace = TRUE)

  # Build summary matrix of groups with allocation info
  Cdf <- as.data.frame(C) %>%
    dplyr::mutate(Size    = rowSums(C),
                  ID      = 1:nrow(C),
                  movable = (1:nrow(C)) %in% vargroups,
                  Split   = alloc)


  # Extract candidate splitting info and evaluation
  spleva <- spliteval(Cdf, delta, w)

  # Split with largest deviation from desired size/balance
  Erridx <- order(spleva$Error, decreasing = TRUE)

  # Is the worst-case split too big or too small?
  direction <- sign(spleva$Size.dev[Erridx[1]])

  OF.track  <- spleva$Error.total
  for(iter in 1:maxiter){
    toBR <- FALSE
    for (i in 1:(length(Erridx) - 1)){
      for (j in (i+1):length(Erridx)){
        Cdf.try <- move_give(Cdf,
                             Erridx[i], Erridx[j],
                             direction, delta, w)


        if(Cdf.try$OF < spleva$Error.total){
          Cdf <- Cdf.try$Cdf

          # Extract candidate splitting info and evaluation
          spleva <- spliteval(Cdf, delta, w)

          # Split with largest deviation from desired size/balance
          Erridx <- order(spleva$Error, decreasing = TRUE)

          # Track progress
          OF.track  <- c(OF.track, spleva$Error.total)

          toBR <- TRUE
          break
        }
      }

      if(toBR) break
    }

    if((i == length(Erridx) - 1 && j == length(Erridx))) {
      break
    }
  }

  for (i in 1:nrow(X)){
    X[i, which(Cdf$Split == i)] <- 1
  }

  attr(X, "opt.history") <- OF.track
  attr(X, "split.eval")  <- spleva
  return(X)
}

calcsplitsummary <- function(Cdf, delta){
  summ <- Cdf %>%
    dplyr::group_by(dplyr::across(dplyr::all_of("Split"))) %>%
    dplyr::summarise(dplyr::across(dplyr::starts_with("class"), sum), .groups = "drop") %>%
    dplyr::mutate(Size = rowSums(dplyr::across(dplyr::starts_with("class"))))

  # correct cases with counts of zero
  idx <- which(!((1:length(delta)) %in% summ$Split))

  if(length(idx) > 0){
    for (i in seq_along(idx)){
      summ <- rbind(summ, summ[1, ])
      summ$Split[nrow(summ)] <- idx[i]
      summ[nrow(summ), -1] <- 0
    }
  }

  return(summ[order(summ$Split), ])

}

spliteval <- function(Cdf, delta, w){
  splsum <- calcsplitsummary(Cdf, delta)

  sizedev    <- splsum$Size/(sum(splsum$Size)+1e-9) - delta
  trgt.prop  <- matrix(colSums(splsum[, -c(1, ncol(splsum))]) / (sum(splsum$Size)+1e-9),
                       nrow = nrow(splsum), ncol = ncol(splsum) - 2, byrow = TRUE)
  act.prop   <- splsum[, -c(1, ncol(splsum))] / matrix(rowSums(splsum[, -c(1, ncol(splsum))]) + 1e-9,
                                                       nrow = nrow(splsum), ncol = ncol(trgt.prop), byrow = FALSE)
  propdev    <- act.prop - trgt.prop

  Objs       <- data.frame(size = abs(sizedev),
                           Balance = apply(abs(propdev), 1, max))
  Error      <- w[1] * Objs$size + w[2] * Objs$Balance
  ErrDriver  <- ifelse(w[1] * Objs$size > w[2] * Objs$Balance, "Size", "Balance")

  return(list(Size.dev = sizedev,
              Prop.dev = propdev,
              Objectives = Objs,
              Error = Error,
              Error.total = sum(Error),
              Error.driver = ErrDriver))
}


move_give <- function(Cdf.try, s1, s2, direction, delta, w){

  if(direction == 1) {
    isource <- which(Cdf.try$Split == s1 & Cdf.try$movable)
    target  <- s2
  } else {
    isource <- which(Cdf.try$Split == s2 & Cdf.try$movable)
    target  <- s1
  }

  trials <- sapply(X = isource,
                   function(ii, target, Y, delta, w){
                     Y$Split[ii] <- target
                     spliteval(Y, delta, w)$Error.total
                   },
                   target = target, Y = Cdf.try, delta = delta, w = w)

  bestmove <- isource[which.min(trials)]
  bestOF   <- min(trials)

  Cdf.try$Split[bestmove] <- target
  return(list(Cdf = Cdf.try, OF = bestOF))
}
