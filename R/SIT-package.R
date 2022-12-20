## usethis namespace: start
#' @useDynLib SIT, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @importFrom stats complete.cases pnorm runif sd
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("SIT", libpath)
}
