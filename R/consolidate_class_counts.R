#' Consolidate ID, group and class count data
#'
#' This function consolidates the output of [extract_clusters()] or
#' [extract_clusters_cdhit()] with data related to the count of class examples
#' to be associated with each entry.
#'
#' @param groups data frame containing entry IDs (column `ID`) and their
#' respective group identifiers (column `Cluster`).
#' @param class_counts data frame containing entry IDs (column `ID` - may
#' contain repeated values, in which case they are aggregated using the sum of
#' each class label) and either:
#' (i) a single column named `Class` with the class label associated with the
#' respective ID, or (ii) two or more columns (with any arbitrary names
#' **except** `Class`. `Class.A`, `Class.B` etc. are permitted) each containing
#' the counts of different classes associated with the respective ID.
#' Note that in the second case there should be no other columns except `ID` and
#' the class columns. See "Examples".

#' @return A matrix with one row per group, and the counts of each unique class
#'    (with names taken either from the unique class values in
#'    ` class_counts$Class` or from the column names of `class_counts`) as
#'    columns.
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
#' }
#'

consolidate_class_counts <- function(groups, class_counts){

  # =======================================================================
  # Sanity checks and initial definitions

  assertthat::assert_that(is.data.frame(groups),
                          all(c("ID", "Cluster") %in% names(groups)),
                          is.data.frame(class_counts),
                          "ID" %in% names(class_counts))


  # =======================================================================

  X <- dplyr::left_join(groups, class_counts, by = "ID")

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
      dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(.x))) %>%
      dplyr::rename_with(.cols = !dplyr::starts_with("Cluster"),
                         .fn = ~paste0("Class.", .x))

  }

  class.balance <- colSums(X[, grep("Class", names(X))]) / sum(X[, grep("Class", names(X))])

  # tmp <- X[, -1] / rowSums(X[ ,-1])
  # names(tmp) <- gsub("Class", "Prop", names(tmp))
  # X <- cbind(X, tmp)

  return(as.matrix(X[, -1]))
}

