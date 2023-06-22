#'
#' @importFrom dplyr %>%
#'
objective_size_deviation <- function(alloc_vector, target_alloc, class_counts, ...){

  my_alloc <- class_counts %>%
    dplyr::bind_cols(data.frame(Split = alloc_vector)) %>%
    dplyr::group_by(Split) %>%
    dplyr::summarise(dplyr::across(!dplyr::starts_with("Cluster"),
                                   ~sum(.x)))

    my_alloc <- data.frame(Split = my_alloc$Split,
                           Sk    = rowSums(my_alloc[, -1]),
                           Target = target_alloc * sum(my_alloc[, -1]))





}
