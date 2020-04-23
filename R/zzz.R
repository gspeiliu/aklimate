.onLoad <- function(libname, pkgname) {
  ##create a composite binary operator in global environment
  `%docomb%` <<- if(foreach::getDoParRegistered()) foreach::`%dopar%` else foreach::`%do%`
}
