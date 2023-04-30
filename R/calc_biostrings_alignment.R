#' Calculate dissimilarity scores based on Smith-Waterman or Needleman-Wunsch
#'
#' This function invokes [Biostrings::pairwiseAlignment()] to calculate the
#' normalised dissimilarity scores between all pairs of sequences in a
#' data frame.
#'
#' @param X data frame with two fields, `IDs` (with sequence ids) and `SEQs`
#' (containing strings with the sequences to be aligned)
#' @param ncpus number of cores to use
#' @param seqtype type of sequence being aligned. Accepts "aa", "dna" or "rna".
#' @param par.list list object with parameters to be passed to
#' [Biostrings::pairwiseAlignment()]. See `Details`.
#' @param vrb logical flag: should progress be printed to console?
#'
#' @return list object with two elements: `scores` (matrix of
#' local or global alignment scores) and `diss_matrix` (dissimilarity matrix)

calc_biostrings_alignment <- function(X,
                                      ncpus = 1,
                                      seqtype = c("aa","dna","rna"),
                                      par.list = list(),
                                      vrb = TRUE){

  # ========================================================================
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(X),
                          all(c("IDs", "SEQs") %in% names(X)),
                          assertthat::is.count(ncpus),
                          is.character(seqtype), length(seqtype) == 1,
                          seqtype %in% c("aa","dna","rna"),
                          is.list(par.list),
                          is.logical(vrb), length(vrb) == 1)

  # Check/Set default par.list attributes
  if(!("type" %in% names(par.list)))
    par.list$type <- "local"
  if(!("substitutionMatrix" %in% names(par.list)))
    par.list$substitutionMatrix <- "BLOSUM62"
  if(!("gapOpening" %in% names(par.list)))
    par.list$gapOpening <- 10
  if(!("gapExtension" %in% names(par.list)))
    par.list$gapExtension <- 4

  # ========================================================================

  if(seqtype == "aa"){
    mymsg("Calculating similarities", vrb)

    utils::data(list    = par.list$substitutionMatrix,
                package = "Biostrings")

    toxp <- list(substitution_matrix = par.list$substitutionMatrix)

    scores <- mypblapply(X   = seq_along(X$SEQs),
                         FUN = myalign,
                         ncpus = ncpus,
                         toexport = toxp,
                         SEQs = X$SEQs,
                         SM = par.list$substitutionMatrix,
                         type = par.list$type) %>%
      dplyr::bind_cols() %>%
      as.matrix()
  }

  # Build denominator matrix: D_{ij} = min(scores_{i,i}, scores{j,j})
  denom <- matrix(pmin(rep(diag(scores), times = nrow(scores)),
                       rep(diag(scores), each = nrow(scores))),
                  nrow  = nrow(scores), byrow = FALSE)

  # Calculate normalized dissimilarity
  rownames(scores) <- colnames(scores) <- X$IDs
  diss_matrix <- 1 - scores / denom

  return(list(scores = scores,
              diss_matrix = diss_matrix))

}
