% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_seq_dissimilarities.R
\name{calc_seq_dissimilarities}
\alias{calc_seq_dissimilarities}
\title{Calculate dissimilarity scores for biological sequences}
\usage{
calc_seq_dissimilarities(
  X = NULL,
  seqfile = NULL,
  ncpus = 1,
  aligner = "diamond",
  par.list = NULL,
  diamond.path = "./",
  min.hit = 8,
  cleanup = TRUE,
  vrb = TRUE
)
}
\arguments{
\item{X}{data frame with two fields, `ID` (with sequence ids) and `SEQ`
(containing strings with the sequences to be aligned). Ignored if a file
path is provided in `seqfile`.}

\item{seqfile}{FASTA file containing the sequences. If `NULL` then sequences
must be provided as a data frame in `X`.}

\item{ncpus}{number of cores to use.}

\item{aligner}{aligner to use. Accepts `"SW"` (Smith-Waterman local
alignment), `"NW"` (Needleman-Wunsch global alignment) or `"diamond"`
(DIAMOND BLASTp aligner).}

\item{par.list}{parameters to be passed to [calc_biostrings_alignment()]
(if `aligner = "SW"` or `aligner = "NW"`) or to [calc_diamond_alignment()]
(if `aligner = "diamond"`). Please check these functions for details.
`NULL` defaults to the standard parameter values of each function.}

\item{diamond.path}{path to the folder where DIAMOND can be executed
(e.g., where the executable file is located). This is also the folder where
the DIAMOND files will be (temporarily) saved.}

\item{min.hit}{shortest allowed length of alignment hit for DIAMOND.}

\item{cleanup}{logical flag: should files created by DIAMOND be deleted at
the end?}

\item{vrb}{logical flag: should progress be printed to console? Note that
this does not control the echoing of DIAMOND (check
[calc_diamond_alignment()]) for details).}
}
\value{
list object with two elements: `scores` (matrix of
alignment scores) and `diss_matrix` (resulting dissimilarity matrix).
}
\description{
This function calculates either a heuristic (using DIAMOND/BLAST) or exact
alignment scores (Smith-Waterman or Needleman-Wunsch) and uses it to
generate a dissimilarity matrix for all pairs of sequences passed as
its argument. Please check [calc_biostrings_alignment()] (Smith-Waterman
and Needleman-Wunsch) or [calc_diamond_alignment()] (DIAMOND/BLAST) for
details. Note that you may need to install DIAMOND in your system (see
below).
}
\section{Installing DIAMOND}{

This routine may require DIAMOND to be available in the system. Please refer
to [https://github.com/bbuchfink/diamond/wiki/2.-Installation](https://github.com/bbuchfink/diamond/wiki/2.-Installation)
for installation details.
}

\examples{
\dontrun{
x <- calc_seq_dissimilarities(seqfile = system.file("diamond", "bfv_proteins.fa", package="moses"),
                              aligner = "SW")
}

}
