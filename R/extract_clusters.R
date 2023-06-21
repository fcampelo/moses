
extract_clusters <- function(diss_matrix,
                             diss_threshold = .3,
                             linkage = "single",
                             vrb = TRUE){

  # =======================================================================
  # Sanity checks and initial definitions
  assertthat::assert_that(is.matrix(diss_matrix),
                          nrow(diss_matrix) == ncol(diss_matrix),
                          rownames(diss_matrix) == colnames(diss_matrix),
                          is.numeric(diss_matrix),
                          is.numeric(diss_threshold),
                          length(diss_threshold) == 1,
                          diss_threshold >= 0, diss_threshold <= 1,
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
              cl.labels = cluster.labels))

}
