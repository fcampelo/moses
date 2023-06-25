#' Consolidate ID, cluster and class count data
#'
#' This function consolidates the output of [extract_clusters()] or
#' [extract_cdhit_clusters()] with data related to the count of class examples
#' to be associated with each entry.
#'
#' @param clusters data frame containing entry IDs (column `ID`) and their
#' respective cluster identifiers (column `Cluster`).
#' @param class_counts data frame containing entry IDs (column `ID` - may
#' contain repeated values, in which case they are aggregated using the sum of
#' each class label) and either:
#' (i) a single column named `Class` with the class label associated with the
#' respective ID, or (ii) two or more columns (with any arbitrary names
#' **except** `Class`. `Class.A`, `Class.B` etc. are permitted) each containing
#' the counts of different classes associated with the respective ID.
#' Note that in the second case there should be no other columns except `ID` and
#' the class columns. See "Examples".

#' @return A list object containing:
#' \itemize{
#'    \item A data frame with one row per cluster. Contains a column
#' `   Cluster` (with the cluster number/ID) and additional columns with
#'     counts of occurrences for each unique class (with
#'     names taken either from the unique class values in
#' `   class_counts$Class` or from the column names of `class_counts`),
#'     plus columns with the class proportions in each cluster and a
#'     column with the total cluster size.
#'    \item A vector with the proportions of each class in the data.
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(moses)
#'
#' # using any fasta file.
#' x <- calc_seq_dissimilarities(seqfile = "diamond/bfv_proteins.fa",
#'                               aligner = "SW")
#' cl <- extract_clusters(x$diss_matrix)
#'
#' epits <- readRDS("diamond/bfv_peptides.rds")
#'
#' cc <- consolidate_class_counts(cl$clusters, epits)
#' }
#'
#'
consolidate_class_counts <- function(clusters, class_counts){

  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.data.frame(clusters),
                          all(c("ID", "Cluster") %in% names(clusters)),
                          is.data.frame(class_counts),
                          "ID" %in% names(class_counts))


  # =======================================================================

  X <- dplyr::left_join(clusters, class_counts, by = "ID")

  if("Class" %in% names(class_counts)){
    X <- X %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("Cluster", "Class")))) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = "Class", values_from = "Count") %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), .fn = ~ifelse(is.na(.x), 0, .x))) %>%
      dplyr::rename_with(.cols = !dplyr::starts_with("Cluster"),
                         .fn = ~paste0("Class.", .x))
  } else {
    X <- X %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("Cluster")))) %>%
      dplyr::summarise(across(everything(), ~sum(.x))) %>%
      dplyr::rename_with(.cols = !dplyr::starts_with("Cluster"),
                         .fn = ~paste0("Class.", .x))

  }

  class.balance <- colSums(X[, grep("Class", names(X))]) / sum(X[, grep("Class", names(X))])

  tmp <- X[, -1] / rowSums(X[ ,-1])
  names(tmp) <- gsub("Class", "Prop", names(tmp))
  X <- cbind(X, tmp)

  X$Size <- rowSums(X[, grep("Class\\.", names(X))])

  return(list(class_counts  = X,
              class_balance = class.balance))
}
