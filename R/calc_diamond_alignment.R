#' Sequence dissimilarity scores based on DIAMOND-BLAST alignment
#'
#' This function invokes DIAMOND (needs to be installed externally) to
#' calculate the normalised dissimilarity scores between all pairs of
#' sequences in a data frame.
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
#'  \item `--sensitive` (more sensitive setting)
#'  \item `--block-size 0.5`   (to limit memory usage to about 3GB, which should enable it to run in most personal machines)
#'  \item `--matrix BLOSUM62` (DIAMOND standard)
#'  \item `--gapopen 11` (DIAMOND standard)
#'  \item `--gapextend 1` (DIAMOND standard)
#' }
#'
#' Additionally, this routine will translate the value of `ncpus`
#' to the DIAMOND option `--threads ncpus`.
#'
#' The user can pass any desired parameters for DIAMOND, except
#' `--threads` (which is passed separately as `ncpus`) and `--outfmt`
#' (which is always set to 6, the BLAST tabular format with
#' the default output fields defined in DIAMOND:
#' "qseqid", "sseqid", "pident", "length",  "mismatch", "gapopen",
#' "qstart", "qend", "sstart", "send", "evalue", "bitscore")
#'
#' @section Dissimilarity value:
#' The dissimilarity values returned by this function are calculated as:
#'
#' diss(a, b) = 1 - max(pident(a, b)) / 100,
#'
#' that is, they are based on the largest percent identity (as returned by
#' DIAMOND) identified between two sequences "a" and "b". Note that only hits
#' longer than `min.hit` are considered.
#'
#'
#' @param X data frame with two fields, `IDs` (with sequence ids) and `SEQs`
#' (containing strings with the sequences to be aligned). Ignored if a file
#' path is provided in `seqfile`.
#' @param seqfile FASTA file containing the sequences. If `NULL` then sequences
#' must be provided as a data frame in `X`.
#' @param ncpus number of cores to use
#' @param par.list list object with parameters to be passed to
#' DIAMOND. See `Details`. A `NULL` par.list is translated to the default values.
#' @param vrb logical flag: should progress be printed to console? Note that
#' this does not control the echoing of DIAMOND - for that, add option
#' `"--quiet"` to `par.list`.
#' @param diamond.path path to the folder where DIAMOND can be executed
#' (e.g., where the executable file is located). This is also the folder where
#' the DIAMOND files will be (temporarily) saved.
#' @param cleanup logical flag: should files created by DIAMOND be deleted at
#' the end?
#' @param min.hit shortest allowed length of alignment hit.
#'
#' @return list object with two elements: `scores` (matrix of
#' local or global alignment scores) and `diss_matrix` (dissimilarity matrix)
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export

calc_diamond_alignment <- function(X = NULL,
                                   seqfile = NULL,
                                   ncpus = 1,
                                   par.list = c("--sensitive",
                                                "--matrix BLOSUM62",
                                                "--gapopen 11",
                                                "--gapextend 1",
                                                "--block-size 0.5",
                                                "--quiet"),
                                   diamond.path = "./",
                                   min.hit = 8,
                                   cleanup = TRUE,
                                   vrb = TRUE){

  # ========================================================================
  # Sanity checks and initial definitions
  if(is.list(par.list)) unlist(par.list)

  assertthat::assert_that(is.data.frame(X) || is.null(X),
                          is.null(X) || all(c("IDs", "SEQs") %in% names(X)),
                          is.null(seqfile) || is.character(seqfile),
                          is.null(seqfile) || file.exists(seqfile),
                          assertthat::is.count(ncpus),
                          is.null(par.list) || is.character(par.list),
                          is.character(diamond.path), length(diamond.path) == 1,
                          dir.exists(diamond.path),
                          is.logical(cleanup), length(cleanup) == 1,
                          is.logical(vrb), length(vrb) == 1,
                          assertthat::is.count(min.hit))

  if(is.null(par.list)){
    par.list <- c("--sensitive", "--matrix BLOSUM62", "--gapopen 11",
                 "--gapextend 1", "--block-size 0.5", "--quiet")
  }

  par.list <- c(par.list, paste0("--threads ", ncpus))

  diamond.path <- gsub("\\/$", "", diamond.path)

  # List of files that will be created (for later deletion if needed)
  torm <- paste0(diamond.path,
                 c("/seqs.fa",                  # sequences file
                   "/proteins-reference.dmnd",  # database file
                   "/protein-matches.tsv"))

  # ========================================================================
  # Run DIAMOND

  # Save sequences as a FASTA file
  if(is.null(seqfile)){
    seqinr::write.fasta(sequences = as.list(X$SEQs),
                        names     = X$IDs,
                        file.out  = torm[1])
  } else {
    file.copy(from = seqfile, to = torm[1], overwrite = TRUE)
  }


  # Call makedb
  mymsg("Building DIAMOND database", vrb)
  mdbst <- paste0(diamond.path, "/diamond makedb --in ", torm[1],
                  " -d ", torm[2],
                  ifelse(any(grepl("--quiet", par.list)), " --quiet ", ""),
                  "--threads ", ncpus)

  system(mdbst)

  # Call aligner
  mymsg("Running DIAMOND aligner", vrb)
  callst <- paste0(diamond.path, "/diamond blastp ",
                   "-q ", torm[1], " -d ", torm[2], " -o ", torm[3], " ",
                   paste(par.list, collapse = " "))

  system(callst)

  # ========================================================================
  mymsg("Calculating dissimilarity matrix", vrb)
  if(!file.exists(torm[3])){
    # No alignments found.
    scores <- as.data.frame(matrix(ncol=12, nrow=1))
    names(scores) <- c("qseqid", "sseqid", "pident", "length",
                       "mismatch", "gapopen", "qstart", "qend",
                       "sstart", "send", "evalue", "bitscore")
    scores <- scores[-1, ]
    protIDs <- names(seqinr::read.fasta(torm[1], as.string = TRUE))
    diss_matrix <- matrix(1, ncol = length(protIDs), nrow = length(protIDs))
    diag(diss_matrix) <- 0
    rownames(diss_matrix) <- colnames(diss_matrix) <- protIDs

  } else {
    scores <- utils::read.csv(torm[3], sep = "\t",
                       header = FALSE,
                       stringsAsFactors = FALSE)
    names(scores) <- c("qseqid", "sseqid", "pident", "length",
                       "mismatch", "gapopen", "qstart", "qend",
                       "sstart", "send", "evalue", "bitscore")

    diss_matrix <- scores %>%
      dplyr::filter(.data$length >= 8) %>%
      dplyr::group_by(.data$qseqid, .data$sseqid) %>%
      dplyr::arrange(dplyr::desc(.data$pident)) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::first),
                       .groups = "drop") %>%
      dplyr::mutate(diss = 1 - .data$pident / 100) %>%
      dplyr::select(.data$qseqid, .data$sseqid, .data$diss) %>%
      tidyr::pivot_wider(names_from = .data$sseqid,
                         values_from = .data$diss,
                         values_fill = NA) %>%
      dplyr::arrange(.data$qseqid) %>%
      as.data.frame()

    rownames(diss_matrix) <- diss_matrix$qseqid

    # Make the dissimilarity scores matrix
    diss_matrix <- diss_matrix %>%
      dplyr::select(order(colnames(diss_matrix)),
                    -.data$qseqid)

    pnames <- rownames(diss_matrix)

    protIDs <- names(seqinr::read.fasta(torm[1], as.string = TRUE))
    missing <- protIDs[which(!(protIDs %in% rownames(diss_matrix)))]

    diss_matrix[(nrow(diss_matrix) + 1):(nrow(diss_matrix) + length(missing)), ] <- 1
    diss_matrix[, (ncol(diss_matrix) + 1):(ncol(diss_matrix) + length(missing))] <- 1

    diss_matrix <- as.matrix(diss_matrix)
    diag(diss_matrix) <- 0

    colnames(diss_matrix) <- c(pnames, missing)
    rownames(diss_matrix) <- c(pnames, missing)
    diss_matrix[is.na(diss_matrix)] <- 1
  }

  if(cleanup){
    mymsg("Removing DIAMOND files.", vrb)
    tmp <- file.remove(torm)
  }

  return(list(scores = scores,
              diss_matrix = diss_matrix))
}
