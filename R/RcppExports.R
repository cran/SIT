# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute the block-wise sum of a vector.
#'
#'
#' @param r An integer vector
#' @param c The number of observations in each block
#' @return The function returns the block sum of the vector.
#' @export
blocksum <- function(r, c) {
    .Call(`_SIT_blocksum`, r, c)
}

