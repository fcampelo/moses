#' Calculate dissimilarity scores based on DIAMOND-BLAST alignment
#'
#' This function invokes DIAMOND (needs to be installed externally) to
#'  calculate the
#' normalised dissimilarity scores between all pairs of sequences in a
#' data frame.
#'
#' @section Installing DIAMOND:
#' This routine requires DIAMOND to be available in the system. Please refer to
#' [https://github.com/bbuchfink/diamond/wiki/2.-Installation](https://github.com/bbuchfink/diamond/wiki/2.-Installation) for installation details.
#'
#'
#' @section Details:
#' Parameter `par.list` should contain a character vector or a list of options
#' that the user wants to pass to DIAMOND.
#' By default, this list will be populated with:
#'
#' \itemize{
#'  \item `--ultra-sensitive` (the most sensitive setting)
#'  \item `--block-size 0.5`   (to limit memory usage to about 3GB, which should enable it to run in most personal machines)
#'  \item `--matrix BLOSUM62` (DIAMOND standard)
#'  \item `--gapopen 11` (DIAMOND standard)
#'  \item `--gapextend 1` (DIAMOND standard)
#' }
#'
#' Additionally, this routine will translate the value of `ncpus`
#' to the DIAMOND option `--threads ncpus`. Similarly, setting `vrb` to
#' `TRUE` will be translated as option `--quiet` for DIAMOND.
#'
#' The user can pass any desired parameters for DIAMOND, except `--outfmt`
#' (which is always set to 6, the BLAST tabular format with
#' the default output fields defined in DIAMOND:
#' "qseqid", "sseqid", "pident", "length",  "mismatch", "gapopen",
#' "qstart", "qend", "sstart", "send", "evalue", "bitscore")
#'
#' @param X data frame with two fields, `IDs` (with sequence ids) and `SEQs`
#' (containing strings with the sequences to be aligned)
#' @param ncpus number of cores to use
#' @param seqtype type of sequence being aligned. Accepts "aa", "dna" or "rna".
#' @param par.list list object with parameters to be passed to
#' DIAMOND. See `Details`.
#' @param vrb logical flag: should progress be printed to console?
#' @param diamond.path path to the folder where DIAMOND can be executed
#' (e.g., where the executable file is located). This is also the folder where
#' the DIAMOND files will be (temporarily) saved.
#' @param cleanup logical flag: should files created by DIAMOND be deleted at
#' the end?
#'
#' @return list object with two elements: `scores` (matrix of
#' local or global alignment scores) and `diss_matrix` (dissimilarity matrix)

calc_diamond_alignment <- function(X,
                                   ncpus = 1,
                                   seqtype = c("aa","dna","rna"),
                                   par.list = c("--ultra-sensitive",
                                                "--matrix BLOSUM62",
                                                "--gapopen 11",
                                                "--gapextend 1",
                                                "--block-size 0.5"),
                                   diamond.path = "./",
                                   cleanup = TRUE,
                                   vrb = TRUE){

  # ========================================================================
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(X),
                          all(c("IDs", "SEQs") %in% names(X)),
                          assertthat::is.count(ncpus),
                          is.character(seqtype), length(seqtype) == 1,
                          seqtype %in% c("aa","dna","rna"),
                          is.list(par.list) || is.character(par.list),
                          is.character(diamond.path), length(diamond.path) == 1,
                          dir.exists(diamond.path),
                          is.logical(cleanup), length(cleanup) == 1,
                          is.logical(vrb), length(vrb) == 1)

  # # ========================================================================
  diamond.path <- gsub("\\/$", "", diamond.path)
  torm <- paste0(diamond.path,
                 c("/seqs.fa",                  # sequences file
                   "/proteins-reference"))      # database file



  # Save sequences as a FASTA file
  seqinr::write.fasta(sequences = as.list(X$SEQs),
                      names     = X$IDs,
                      file.out  = torm[1])

  # Build DIAMOND makedb string
  mdbst <- paste0(diamond.path, "/diamond makedb --in ", torm[1],
                  " -d ", torm[2])

  system("diamond makedb --in ../data/proteins.fa -d ../data/diamond/proteins-reference")

  callst <- paste0(diamond.path, "/diamond ",
                   ifelse(seqtype == "aa", "blastp", "blastx"),
                   )


    #
    # if(seqtype == "aa"){
    #   mymsg("Calculating similarities", vrb)
    #
    #   utils::data(list    = par.list$substitutionMatrix,
    #               package = "Biostrings")
    #
    #   toxp <- list(substitution_matrix = par.list$substitutionMatrix)
    #
  #   scores <- mypblapply(X   = seq_along(X$SEQs),
  #                        FUN = myalign,
  #                        ncpus = ncpus,
  #                        toexport = toxp,
  #                        vrb = vrb,
  #                        SEQs = X$SEQs,
  #                        pars = par.list) %>%
  #     dplyr::bind_rows() %>%
  #     t() %>%
  #     as.matrix()
  # }
  #
  # # Build denominator matrix: D_{ij} = min(scores_{i,i}, scores{j,j})
  # denom <- matrix(pmin(rep(diag(scores), times = nrow(scores)),
  #                      rep(diag(scores), each = nrow(scores))),
  #                 nrow  = nrow(scores), byrow = FALSE)
  #
  # # Calculate normalized dissimilarity
  # rownames(scores) <- colnames(scores) <- X$IDs
  # diss_matrix <- 1 - scores / denom
  #
  # return(list(scores = scores,
  #             diss_matrix = diss_matrix))

}
