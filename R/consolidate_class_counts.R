#' Consolidate ID, cluster and class count data
#'
#' This function consolidates the output of [extract_clusters()] or
#' [extract_cdhit_clusters()] with data related to the count of class examples
#' to be associated with each entry.
#'
#' @param clusters data frame containing entry IDs (column `IDs`) and their
#' respective cluster identifiers (column `Cluster`).
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
#' @return list object with two elements: `clusters` (An object of class
#' 'hclust'. It encodes a stepwise dendrogram.) and `cl.labels` (a vector with
#' cluster memberships at the dissiimlarity level `diss_threshold`). Notice that
#' this is the output of a call to [stats::cutree()], so if a vector is passed
#' in `diss_threshold` then `cl.labels` will be a matrix.
#'
#'
#' @export
#'
consolidate_class_counts <- function(X = NULL,
                                     seqfile = NULL,
                                     diss_threshold = .3,
                                     vrb = TRUE,
                                     ncpus = 1,
                                     par.list = list(kmerSize = NULL,
                                                     min_length = 8,
                                                     s = ,
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

  return(clusters  = clusters,
         cl.labels = data.frame(IDs = clusters$query_name,
                                Cluster = clusters$cluster_idx))

}
