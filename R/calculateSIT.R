#' Compute the cross rank coefficient sit on two vectors.
#'
#' This function computes the sit coefficient between two vectors x and y.
#'
#' @aliases  sitcorcoefficient
#' @param x Vector of numeric values in the first coordinate.
#' @param y Vector of numeric values in the second coordinate.
#' @param c The number of observations in each slice.
#' @return The function returns the value of the
#' sit
#' coefficient.
#' @note Auxiliary function with no checks for NA, etc.
#' @author Yilin Zhang, Canyi Chen & Liping Zhu
#' @seealso sitcor
#' @references Zhang Y., Chen C., & Zhu L. (2021). Sliced Independence Test. Statistica Sinica. https://doi.org/10.5705/ss.202021.0203.
#' @keywords ~methods ~htest
#' @export
#' @examples
#' # Compute one of the coefficients
#' library("psychTools")
#' data(peas)
#' calculateSIT(peas$parent,peas$child)
#' calculateSIT(peas$child,peas$parent)
#' @export
calculateSIT <- function(x, y, c = 2) {
  n <- length(x)

  h <- floor(n / c)
  xr <- rank(x, ties.method = "random")
  ix <- order(xr)
  ys <- y[ix]
  r1 <- rank(ys)

  s <- blocksum(r1, c)
  a <- s / n ^ 2 / (c - 1)
  b <- sum(r1) / n / (n - 1) - sum(r1 ^ 2) / n ^ 2 / (n - 1)

  sit <- 1 - a / b
  return(sit)
}
