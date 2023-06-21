#' Sequence dissimilarity scores based on Smith-Waterman or Needleman-Wunsch
#'
#' This function invokes [Biostrings::pairwiseAlignment()] to calculate the
#' normalised dissimilarity scores between all pairs of sequences in a
#' data frame.
#'
#' @section Details:
#' Parameter `par.list` should contain the list of parameters that the user
#' wants to pass down to [Biostrings::pairwiseAlignment()]. By default,
#' this list will be populated with:
#'
#' \itemize{
#'  \item `type = "local"`
#'  \item `substitutionMatrix = "BLOSUM62"`
#'  \item `gapOpening = 10`
#'  \item `gapExtension = 4`
#' }
#'
#' The user can overwrite these values and add any other parameter accepted by
#' [Biostrings::pairwiseAlignment()] _except_ ` scoreOnly` (which is always
#'  set to `TRUE`).
#'
#' @section Dissimilarity value:
#' The dissimilarity values returned by this function are calculated as:
#'
#' diss(a, b) = 1 - score(a, b) / min(score(a, a), score(b, b)),
#'
#' that is, based on the alignment score between a and b normalised by the
#' smallest self-alignment score of sequences a and b.
#'
#' @param X data frame with two fields, `IDs` (with sequence ids) and `SEQs`
#' (containing strings with the sequences to be aligned). Ignored if a file
#' path is provided in `seqfile`.
#' @param seqfile FASTA file containing the sequences. If `NULL` then sequences
#' must be provided as a data frame in `X`.
#' @param ncpus number of cores to use.
#' @param par.list list object with parameters to be passed to
#' [Biostrings::pairwiseAlignment()]. See `Details`.
#' `NULL` is translated to the default values.
#' @param vrb logical flag: should progress be printed to console?
#'
#' @return list object with two elements: `scores` (matrix of
#' local or global alignment scores) and `diss_matrix` (dissimilarity matrix)
#'
#' @importFrom dplyr %>%
#'
#' @export

calc_biostrings_alignment <- function(X = NULL,
                                      seqfile = NULL,
                                      ncpus = 1,
                                      par.list = list(type = "local",
                                                      substitutionMatrix = "BLOSUM62",
                                                      gapOpening = 10,
                                                      gapExtension = 4),
                                      vrb = TRUE){

  # ========================================================================
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(X) || is.null(X),
                          is.null(X) || all(c("IDs", "SEQs") %in% names(X)),
                          is.null(seqfile) || is.character(seqfile),
                          is.null(seqfile) || file.exists(seqfile),
                          assertthat::is.count(ncpus),
                          is.null(par.list) || is.list(par.list),
                          is.logical(vrb), length(vrb) == 1)

  if(is.null(par.list)){
    par.list <- list(type = "local",
                     substitutionMatrix = "BLOSUM62",
                     gapOpening = 10,
                     gapExtension = 4)
  }

  # Check/Set default par.list attributes
  if(!("type" %in% names(par.list)))
    par.list$type <- "local"
  if(!("substitutionMatrix" %in% names(par.list)))
    par.list$substitutionMatrix <- "BLOSUM62"
  if(!("gapOpening" %in% names(par.list)))
    par.list$gapOpening <- 10
  if(!("gapExtension" %in% names(par.list)))
    par.list$gapExtension <- 4

  if(!is.null(seqfile)){
    X <- seqinr::read.fasta(seqfile, as.string = TRUE)
    X <- data.frame(IDs = attributes(X)$name,
                    SEQs = toupper(as.character(X)))
  }

  # ========================================================================

  mymsg("Calculating similarities", vrb)

  utils::data(list    = par.list$substitutionMatrix,
              package = "Biostrings")

  toxp <- list(substitution_matrix = par.list$substitutionMatrix)

  scores <- mypblapply(X   = seq_along(X$SEQs),
                       FUN = myalign,
                       ncpus = ncpus,
                       toexport = toxp,
                       vrb = vrb,
                       SEQs = X$SEQs,
                       pars = par.list) %>%
    dplyr::bind_rows() %>%
    t() %>%
    as.matrix()


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
