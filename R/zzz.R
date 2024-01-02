.onAttach <- function(lib, pkg) {
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage(paste("Welcome to spatialWarsaw version ",ver,sep=""))
}
