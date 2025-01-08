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
#' # Desired allocation proportions
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
#' }


make_splits_rand_refine <- function(C, delta, X0 = NULL, w = c(2/3, 1/3), maxiter = 100){

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


        if(Cdf.try$OF <= spleva$Error.total){
          Cdf <- Cdf.try$Cdf
          spleva <- spliteval(Cdf, delta, w)
          # Extract candidate splitting info and evaluation
          spleva <- spliteval(Cdf, delta, w)

          # Split with largest deviation from desired size/balance
          Erridx <- order(spleva$Error, decreasing = TRUE)

          toBR <- TRUE
          break
        }
      }

      if(toBR) break
    }

    OF.track  <- c(OF.track, spleva$Error.total)

    if((i == length(Erridx) - 1 && j == length(Erridx))) {
      break
    }
  }

  for (i in 1:nrow(X)){
    X[i, which(Cdf$Split == i)] <- 1
  }

  return(X)
}

calcsplitsummary <- function(Cdf){
  Cdf %>%
    dplyr::group_by(Split) %>%
    dplyr::summarise(across(starts_with("class"), sum), .groups = "drop") %>%
    dplyr::mutate(Size = rowSums(across(starts_with("class"))))
}

spliteval <- function(Cdf, delta, w){
  splsum <- calcsplitsummary(Cdf)

  sizedev    <- splsum$Size/sum(splsum$Size) - delta
  trgt.prop  <- matrix(colSums(splsum[, -c(1, ncol(splsum))]) / sum(splsum$Size),
                       nrow = nrow(splsum), ncol = ncol(splsum) - 2, byrow = TRUE)
  act.prop   <- splsum[, -c(1, ncol(splsum))] / matrix(splsum$Size, nrow = nrow(splsum), ncol = length(trgt.prop), byrow = FALSE)
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
