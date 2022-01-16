#Borrowed from Parsons, that also borrowed from another guy (and so on :) )

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage(pkgname, " is stil software in progress! Please report any bugs.")
}
