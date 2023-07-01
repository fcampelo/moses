#'
#' @importFrom dplyr %>%
#'
size_deviation <- function(alloc_vector,
                           class_counts,
                           target_proportions,
                           ...){

  Sk <- class_counts %>%
    dplyr::mutate(Alloc = alloc_vector) %>%
    dplyr::group_by(Alloc) %>%
    dplyr::summarise(Prop = sum(Size) / sum(class_counts$Size)) %>%
    dplyr::arrange(Prop) %>%
    dplyr::mutate(Target = sort(target_proportions),
                  Delta.s = abs(Target - Prop))

  return(list(Sk.max = max(Sk$Delta.s),
              Sk.avg = mean(Sk$Delta.s)))

}


balance_deviation <- function(alloc_vector,
                              class_counts,
                              target_proportions,
                              class_balance,
                              ...){

  Dk <- class_counts %>%
    dplyr::mutate(Alloc = alloc_vector) %>%
    dplyr::group_by(Alloc) %>%
    dplyr::summarise(Size = sum(Size),
                     dplyr::across(dplyr::starts_with("Class"),
                                   .names = "Prop.{col}",
                                   ~sum(.x)/Size)) %>%
    dplyr::arrange(Size) %>%
    dplyr::mutate(Target = sort(target_proportions * sum(class_counts$Size)),
                  Delta.size = Size - Target) %>%
    dplyr::rowwise()



                                return(list(Sk.max = max(Sk$Delta.s),
                                            Sk.avg = mean(Sk$Delta.s)))

}
