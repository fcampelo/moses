

scalar_objective_function <- function(alloc_vector,
                                      class_counts,
                                      target_proportions,
                                      weights = 1,
                                      rho     = 1e-2,
                                      class_balance = NULL){

  # =======================================================================
  # Sanity checks and initial definitions


  idx <- grep("Class.", names(class_counts))
  if(is.null(class_balance)){
    class_balance <- as.data.frame(t(colSums(class_counts[, idx]) / sum(class_counts[, idx])))
  }

  if(!("Size" %in% names(class_counts))){
    class_counts$Size <- rowSums(class_counts[, grep("Class\\.", names(class_counts))])
  }

  weights <- weights / sum(weights)

  target_proportions <- target_proportions / sum(target_proportions)



  # =======================================================================

  splits <- class_counts %>%
    dplyr::mutate(Split  = alloc_vector) %>%
    dplyr::filter(!is.na(Split))

  if (nrow(splits) > 0){

    splits <- splits %>%
      dplyr::group_by(Split) %>%
      dplyr::summarise(C.. = sum(class_counts$Size[!is.na(alloc_vector)]),
                       cardSk = sum(Size),
                       normSk = cardSk / C..,
                       Nk     = dplyr::n(),
                       dplyr::across(dplyr::starts_with("Class"),
                                     sum))

    CCount  <- dplyr::select(splits, dplyr::starts_with("Class"))
    CTarget <- splits$cardSk * cbind(splits[, 1], class_balance)[, -1]
    flag <- CCount < CTarget

    fj  <- flag * (1 - CCount/CTarget) +
      (!flag) * ((CCount - CTarget) / (splits$cardSk - CTarget))

    normSk <- 0 * target_proportions
    normSk[1:length(splits$normSk)] <- splits$normSk
    Delta.Sk     <- abs(normSk - sort(target_proportions))

    Delta.Bk.max <- apply(fj, 1, max)
    Delta.Bk.avg <- rowMeans(fj)

    Nk <- 0 * target_proportions
    Nk[1:length(splits$Nk)] <- splits$Nk
    Hk           <- ifelse(Nk == 0, 1, 1 / splits$Nk)

    obj.funs <- data.frame(prime  = c(max(Delta.Sk),
                                      max(Delta.Bk.max),
                                      max(Hk)),
                           second = c(mean(Delta.Sk),
                                      mean(Delta.Bk.avg),
                                      stats::median(Hk)))
    y1 <- max(weights * obj.funs$prime) +
      rho * sum((weights * obj.funs$prime))
    y2 <- max(weights * obj.funs$second) +
      rho * sum((weights * obj.funs$second))

  } else {
    obj.funs <- data.frame(prime  = c(1, 1, 1),
                           second = c(1, 1, 1))
    y1 <- 10
    y2 <- 10
  }

  return(list(obj.funs = obj.funs,
              f.AT     = y1,
              f.AT2    = y2))

}
