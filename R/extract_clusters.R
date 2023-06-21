#' Extract clusters based on dissimilarity matrix
#'
#' This function calculates a hierarchical clustering of entries based on
#' a dissimilarity matrix (generated, e.g., using [calc_seq_dissimilarities()]).
#' It invokes [fastcluster::hclust()] internally for higher efficiency.
#'
#' @param diss_matrix a dissimilarity matrix.
#' @param diss_threshold the desired dissimilarity cutting threshold for
#' attributing clusters.
#' @param linkage the linkage function to use. This is passed down to
#' [fastcluster::hclust()] as parameter `method`. Please check the
#' documentation of [fastcluster::hclust()] for details.
#' @param vrb logical flag: should progress be printed to console?
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
extract_clusters <- function(diss_matrix,
                             diss_threshold = .3,
                             linkage = "single",
                             vrb = TRUE){

  # =======================================================================
  # Sanity checks and initial definitions
  assertthat::assert_that(is.matrix(diss_matrix),
                          nrow(diss_matrix) == ncol(diss_matrix),
                          all(rownames(diss_matrix) == colnames(diss_matrix)),
                          is.numeric(diss_matrix),
                          is.numeric(diss_threshold),
                          length(diss_threshold) > 0,
                          all(diss_threshold >= 0), all(diss_threshold <= 1),
                          is.character(linkage), length(linkage) == 1,
                          is.logical(vrb), length(vrb) == 1)

  # =======================================================================
  mymsg("Extracting hierarchical clusters", vrb)

  clusters  <- do.call(fastcluster::hclust,
                       args = list(d = stats::as.dist(diss_matrix),
                                   method = linkage))

  mymsg(paste0("Cutting clusters at dissimilarity threshold ",
               diss_threshold), vrb)
  cluster.labels <- stats::cutree(clusters, h = diss_threshold)


  return(list(clusters = clusters,
              cl.labels = data.frame(ID = names(cluster.labels),
                                     Cluster = unname(cluster.labels))))

}
