.onLoad <- function(libname, pkgname) {

  required_pkgs <- c(
    clusterProfiler = "4.13.0",
    igraph          = "2.0",
    textclean       = "0.9.3"
  )

  # Loop through and check versions
  for (pkg in names(required_pkgs)) {
    req_ver <- required_pkgs[[pkg]]

    # Check if installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        sprintf(
          "Package %s requires '%s' (>= %s) to be installed.",
          pkgname, pkg, req_ver
        ),
        call. = FALSE
      )
    }

    # Check version number
    if (utils::packageVersion(pkg) < req_ver) {
      stop(
        sprintf(
          "Package %s requires '%s' version >= %s. Installed version is %s.",
          pkgname,
          pkg,
          req_ver,
          as.character(utils::packageVersion(pkg))
        ),
        call. = FALSE
      )
    }
  }
}
