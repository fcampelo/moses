#' Extract protein sequence clusters using CDHIT
#'
#' This function extract protein sequence clusters using CDHIT
#' (<https://doi.org/10.1093/bioinformatics/btl158>), based on the R
#' function [CellaRepertorium::cdhit()]. Please refer to that function for
#' details.
#'
#' @param X data frame with two fields, `IDs` (with sequence ids) and `SEQs`
#' (containing strings with the sequences to be aligned). Ignored if a file
#' path is provided in `seqfile`.
#' @param seqfile FASTA file containing the sequences. If `NULL` then sequences
#' must be provided as a data frame in `X`.
#' @param diss_threshold the desired dissimilarity cutting threshold for
#' attributing clusters. Passed down to [CellaRepertorium::cdhit()] as
#' parameter `identity = 1 - diss_treshold`.
#' @param vrb logical flag: should progress be printed to console? Passed down
#' to [CellaRepertorium::cdhit()] as parameter `showProgress = vrb`.
#' @param ncpus number of cores to use. Passed down
#' to [CellaRepertorium::cdhit()] as parameter `T = ncpus`.
#' @param par.list further parameters to be passed down to
#' [CellaRepertorium::cdhit()]. Please check that function for details.
#'
#' @return list object with two elements: `cl.df` (A data frame returned by
#' CDHIT) and `clusters` (a data frame
#' with cluster memberships at the dissimilarity level `diss_threshold`).
#'
#'
#' @export
#'
extract_cdhit_clusters <- function(X = NULL,
                                   seqfile = NULL,
                                   diss_threshold = .3,
                                   vrb = TRUE,
                                   ncpus = 1,
                                   par.list = list(kmerSize = NULL,
                                                   min_length = 8,
                                                   s = 1,
                                                   G = 1)){

  # =======================================================================
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(X) || is.null(X),
                          is.null(X) || all(c("IDs", "SEQs") %in% names(X)),
                          is.null(seqfile) || is.character(seqfile),
                          is.null(seqfile) || file.exists(seqfile),
                          assertthat::is.count(ncpus),
                          is.null(par.list) || is.list(par.list),
                          is.logical(vrb), length(vrb) == 1,
                          is.numeric(diss_threshold),
                          length(diss_threshold) > 0,
                          all(diss_threshold >= 0), all(diss_threshold <= 1))

  # =======================================================================
  mymsg("Extracting CDHIT clusters", vrb)

  if(is.null(seqfile)){
    seqs <- X$SEQs
    names(seqs) <- X$IDs
    seqs <- Biostrings::AAStringSet(x = seqs, use.names = TRUE)
  } else {
    seqs <- Biostrings::readAAStringSet(seqfile)
  }

  par.list$T <- ncpus
  par.list$showProgress <- vrb
  par.list$identity <- 1 - diss_threshold
  par.list$seqs <- seqs
  par.list$only_index <- FALSE

  clusters <- do.call(CellaRepertorium::cdhit, args = par.list)

  return(cl.df    = clusters,
         clusters = data.frame(ID      = clusters$query_name,
                               Cluster = clusters$cluster_idx))

}
