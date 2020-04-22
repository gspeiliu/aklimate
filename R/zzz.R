.onLoad <- function(libname, pkgname) {
  `%docomb%` <<- if(foreach::getDoParRegistered()) foreach::`%dopar%` else foreach::`%do%`
}
