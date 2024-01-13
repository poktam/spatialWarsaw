.onAttach <- function(lib, pkg) {
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage(paste("\n","Welcome to spatialWarsaw version ",ver,sep=""))
}
