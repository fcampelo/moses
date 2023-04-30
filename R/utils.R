# Little selective messaging function
mymsg <- function(str, vrb){ if(vrb) message(str) }


# Little field check function
nullcheck <- function(x) { ifelse(is.null(x), yes = NA, no = x) }


# ================================================================
# Parallel processing functions
set_mc <- function(ncpus){
  if (ncpus > 1 && .Platform$OS.type == "windows"){
    cl <- parallel::makeCluster(ncpus, setup_timeout = 2)
  } else {
    cl <- max(1, min(ncpus, parallel::detectCores() - 1))
  }
  return(cl)
}

close_mc <- function(cl){
  # Stop cluster
  if("cluster" %in% class(cl)) parallel::stopCluster(cl)
  invisible(TRUE)
}


mypblapply <- function(X, FUN, ncpus, toexport = list(), ...){
  cl  <- set_mc(ncpus)

  if(ncpus > 1 && length(toexport) > 0 && .Platform$OS.type == "windows"){
    parallel::clusterExport(cl = cl,
                            varlist = toexport)
  }
  res <- pbapply::pblapply(cl = cl, X = X, FUN = FUN, ...)

  close_mc(cl)
  return(res)
}



# ================================================================
# Progress bar function
mypb <- function(i, max_i, t0, npos){
  nb <- max(1, ceiling(max_i / npos))
  if (i == 0){
    pbstr <- paste0("  |", paste(rep("_", npos), collapse = ""), "|")
    cat(pbstr, "0% processed. Elapsed time: 0s")
  } else if (i >= (max_i - 1)) {
    pbstr <- paste(rep(">", times = npos), collapse = "")
    td <- Sys.time() - t0
    perc_done <- 100
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  } else if (!(i %% nb)) {
    nn <- i / nb
    pbstr <- paste(rep(c("+", "_"), times = c(nn, npos - nn)),
                   collapse = "")
    td <- Sys.time() - t0
    perc_done <- round(100 * i / max_i, digits = 0)
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  }
  invisible(NULL)
}


# ================================================================
# Alignment function
myalign <- function(i, SEQs, SM, type){
  #utils::data(list = SM, package = "Biostrings")
  patt <- rep(SEQs[i], times = 1 + length(SEQs) - i)
  subj <- SEQs[i:length(SEQs)]
  vals <- Biostrings::pairwiseAlignment(pattern = patt,
                                        subject = subj,
                                        substitutionMatrix = SM,
                                        type = type,
                                        scoreOnly = TRUE)
  return(c(rep(NA, i - 1), vals))
}
