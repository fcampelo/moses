constructive_heuristic <- function(class_counts,
                                   target_proportions,
                                   weights = 1,
                                   rho     = 1e-2,
                                   class_balance = NULL){

  class_counts <- class_counts[order(class_counts$Size,
                                     decreasing = TRUE), ]
  class_counts$Split <- NA

  alloc <- as.data.frame(matrix(data = NA,
                                nrow = nrow(class_counts),
                                ncol = length(target_proportions),
                                byrow = TRUE))

  for (i in seq_along( class_counts$Cluster)){
    alloc[i, ] <- 1:length(target_proportions)
    f0 <- lapply(alloc,
                 scalar_objective_function,
                 class_counts = class_counts,
                 target_proportions = target_proportions,
                 weights = weights,
                 rho     = rho,
                 class_balance = class_balance)

    # Extract best allocation
    idx <- which(sapply(f0, function(x)x$f.AT) == min(sapply(f0, function(x)x$f.AT)))

    if (length(idx) > 1){
      idx2 <- which(sapply(f0[idx], function(x)x$f.AT2) == min(sapply(f0[idx], function(x)x$f.AT2)))

      if (length(idx2) > 1) idx2 <- sample(idx2, 1)

      idx <- idx[idx2]
    }

    alloc[i, ] <- idx
  }
}
