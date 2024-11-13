.onAttach <- function(...) {

  # This is just to prevent CHECK errors due to knitr not being
  # explicitly called anywhere.
  # if(FALSE) knitr::opts_chunk

  tocheck <- c("Biostrings", "CellaRepertorium")
  vers    <- c('2.60.0',     '1.10.0')

  for (i in seq_along(tocheck)){
    # Check if Biostrings is installed
    bsp <- suppressWarnings(utils::packageDescription(tocheck[i]))

    if(all(is.na(bsp))) {
      packageStartupMessage("\nDependency package ", tocheck[i],
                            " not detected.\n",
                            "Please run install_bioc_dependencies(force = TRUE)\n",
                            "before using this package.")
    } else if(utils::packageVersion(tocheck[i]) < vers[i]) {
      packageStartupMessage("\nPackage  ", tocheck[i],  "version ",
                            vers[i], " or later is required\n",
                            "Please run install_bioc_dependencies(force = TRUE)\n",
                            "before using this package.")
    }

  }
  invisible(TRUE)
}
