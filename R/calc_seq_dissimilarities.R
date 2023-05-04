#' Calculate dissimilarity scores for biological sequences
#'
#' This function calculates either a heuristic (using DIAMOND/BLAST) or exact
#' alignment scores (Smith-Waterman or Needleman-Wunsch) and uses it to
#' generate a dissimilarity matrix for all pairs of sequences passed as
#' its argument. Please check [calc_biostrings_alignment()] (Smith-Waterman
#' and Needleman-Wunsch) or [calc_diamond_alignment()] (DIAMOND/BLAST) for
#' details. Note that you may need to install DIAMOND in your system (see
#' below).
#'
#' @section Installing DIAMOND:
#' This routine may require DIAMOND to be available in the system. Please refer
#' to [https://github.com/bbuchfink/diamond/wiki/2.-Installation](https://github.com/bbuchfink/diamond/wiki/2.-Installation)
#' for installation details.
#'
#' @param X data frame with two fields, `IDs` (with sequence ids) and `SEQs`
#' (containing strings with the sequences to be aligned). Ignored if a file
#' path is provided in `seqfile`.
#' @param seqfile FASTA file containing the sequences. If `NULL` then sequences
#' must be provided as a data frame in `X`.
#' @param ncpus number of cores to use
#' @param seqtype type of sequence being aligned. Accepts "aa", "dna" or "rna".
#' @param par.list parameters to be passed to [calc_biostrings_alignment()] or
#'  [calc_diamond_alignment()]. Please check these functions for details.
#' @param vrb logical flag: should progress be printed to console? Note that
#' this does not control the echoing of DIAMOND (check
#' [calc_diamond_alignment()]) for details).
#' @param diamond.path path to the folder where DIAMOND can be executed
#' (e.g., where the executable file is located). This is also the folder where
#' the DIAMOND files will be (temporarily) saved.
#' @param cleanup logical flag: should files created by DIAMOND be deleted at
#' the end?
#' @param min.hit shortest allowed length of alignment hit for DIAMOND.
#'
#' @return list object with two elements: `scores` (matrix of
#' alignment scores) and `diss_matrix` (resulting dissimilarity matrix).
#'
#'
#' @export

# calc_seq_dissimilarities <- function(X = NULL,
#                                      seqfile = NULL,
#                                      ncpus = 1,
#                                      seqtype = c("aa","dna","rna"),
#                                      disstype = "diamond",
#                                      par.list = NULL,
#                                      diamond.path = "./",
#                                      min.hit = 8,
#                                      cleanup = TRUE,
#                                      vrb = TRUE){

  # ========================================================================
  # Sanity checks and initial definitions
#
#
#   return(list(scores = scores,
#               diss_matrix = diss_matrix))
# }
