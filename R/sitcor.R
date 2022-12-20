#' Conduct the sliced independence test.
#'
#' This function computes the sit coefficient between two vectors x and y,
#' possibly all paired coefficients for a matrix.
#'
#' @aliases  sit sitcor
#'
#' @param x Vector of numeric values in the first coordinate.
#' @param y Vector of numeric values in the second coordinate.
#' @param c The number of observations in each slice.
#' @param pvalue Whether or not to return the p-value of rejecting
#' independence, if TRUE the function also returns the standard deviation of
#' sit.
#' @param ties Do we need to handle ties? If ties=TRUE the algorithm assumes
#' that the data has ties and employs the more elaborated theory for
#' calculating s.d. and P-value. Otherwise, it uses the simpler theory. There
#' is no harm in putting ties = TRUE even if there are no ties.
#' @param method If method = "asymptotic" the function returns P-values
#' computed by the asymptotic theory (not available in the presence of ties). If method = "permutation", a permutation
#' test with nperm permutations is employed to estimate the P-value. Usually,
#' there is no need for the permutation test. The asymptotic theory is good
#' enough.
#' @param nperm In the case of a permutation test, \code{nperm} is the number
#' of permutations to do.
#' @param factor Whether to transform integers into factors, the default is to
#' leave them alone.
#' @return In the case pvalue=FALSE, function returns the value of the sit
#' coefficient, if the input is a matrix, a matrix of coefficients is returned.
#' In the case pvalue=TRUE is chosen, the function returns a list:
#' \describe{\item{sitcor}{The
#' value of the sit coefficient.}
#' \item{sd}{The standard deviation.}
#' \item{pval}{The test p-value.}
#' }
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' library("psychTools")
#' data(peas)
#' # Visualize       the peas data
#' library(ggplot2)
#' ggplot(peas,aes(parent,child)) +
#' geom_count() + scale_radius(range=c(0,5)) +
#'        xlim(c(13.5,24))+ylim(c(13.5,24))+       coord_fixed() +
#'        theme(legend.position="bottom")
#' # Compute one of the coefficients
#' sitcor(peas$parent,peas$child, c = 4, pvalue=TRUE)
#' sitcor(peas$child,peas$parent, c = 4)
#' # Compute all the coefficients
#' sitcor(peas, c = 4)
#'
#' @author Yilin Zhang, Canyi Chen & Liping Zhu
#' @references Zhang Y., Chen C., & Zhu L. (2022). Sliced Independence Test. Statistica Sinica. https://doi.org/10.5705/ss.202021.0203.
#' @keywords ~methods ~htest
#' @export
sitcor <- function(x,
                   y = NULL,
                   c = 2,
                   pvalue = FALSE,
                   ties = FALSE,
                   method = "asymptotic",
                   nperm = 199,
                   factor = FALSE) {
  # Factor variables are converted to integers here.
  if (factor == TRUE) {
    if (!is.numeric(x))
      x <- as.numeric(factor(x))
    if (!is.numeric(y))
      y <- as.numeric(factor(y))
  }
  if (is.data.frame(y))
    y <- as.matrix(y)
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x)))
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y)))
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }



  if (is.null(y)) {
    ncy <- ncx <- ncol(x)
    if (ncx == 0)
      stop("'x' is empty")
    if (pvalue == TRUE)
      stop("testing is not available for matrices")
    r <- matrix(0, nrow = ncx, ncol = ncy)
    for (i in seq_len(ncx)) {
      for (j in seq_len(i)) {
        x2 <- x[, i]
        y2 <- x[, j]
        ok <- complete.cases(x2, y2)
        x2 <- x2[ok]
        y2 <- y2[ok]
        if (any(ok))
        {
          r[i, j] <- calculateSIT(x2, y2, c = c)
          ###it's not symmetric, we have to compute both
          r[j, i] <- calculateSIT(y2, x2, c = c)
        }
        else
          NA
      }
    }
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(x)
    return(r)
  } else {
    if (ncol(as.matrix(x)) > 1 & ncol(as.matrix(y)) == 1) {
      if (pvalue == TRUE)
        stop("testing is not available for matrices and vectors pairs")
      r <- rep(NA, ncol(x))
      for (i in seq_len(ncol(x))) {
        x2 <- x[, i]
        ok <- complete.cases(cbind(x2, y))
        x2 <- x2[ok]
        y2 <- y[ok]
        if (any(ok)) {
          r[i] <- calculateSIT(x2, y2, c = c)
        } else {
          NA
        }
      }
      names(r) <- colnames(x)
      return(r)
    }
    # Two vectors case
    if (ncol(as.matrix(x)) == 1 & ncol(as.matrix(y)) == 1) {
      ok <- complete.cases(x, y)
      x <- x[ok]
      y <- y[ok]
      sit <- calculateSIT(x, y, c = c)
      n <- length(x)
    }
  }



  # Arguments checking
  stopifnot(length(x) == length(y))
  # stopifnot(h <= floor(length(x) / 2))

  sitcor <- calculateSIT(x, y, c = c)
  n <- length(x)
  # If P-value needs to be computed:
  if (pvalue) {
    if (ties == FALSE) {
      var1 <- 2 / 5 * 2 / (c - 1) / n
      return(list(
        sitcor = sitcor,
        sd = sqrt(var1),
        pval = 1 - pnorm(sqrt(n) * sitcor / sqrt(var1))
      ))
    } else {
      # Arguments checking
      if (!(method %in% c("asymptotic", "permutation")))
        stop("method for test can only be asymptotic or permutation")

      # If there are ties, and the theoretical variance is used:
      if (method == "asymptotic") {
        # stop("Asymptotic methods have not been available in the presence of ties.")
        # The following steps calculate the theoretical variance in the presence of ties:
        PI <- rank(x, ties.method = "random")
        # fr[i] is number of j s.t. y[j] <= y[i], divided by n.
        fr <- rank(y, ties.method = "max")/n
        # gr[i] is number of j s.t. y[j] >= y[i], divided by n.
        gr <- rank((- y), ties.method = "max")/n
        ord <- order(PI)
        fr <- fr[ord]

        CU <- mean(gr* (1 - gr))

        qfr <- sort(fr)
        ind <- c(1:n)
        ind2 <- 2*n - 2*ind + 1
        ai <- mean(ind2*qfr*qfr)/n
        ci <- mean(ind2*qfr)/n
        cq <- cumsum(qfr)
        m <- (cq + (n - ind)*qfr)/n
        b <- mean(m^2)
        v <- (ai - 2*b + ci^2)/(CU^2)

        var1 <- v * 2 / (c - 1) / n
        return(list(
          sitcor = sitcor,
          sd = sqrt(var1),
          pval = 1 - pnorm(sqrt(n) * sitcor / sqrt(var1))
        ))


      }
      #
      # If the permutation test is to be used for calculating P-value:
      if (method == "permutation") {
        rp <- rep(0, nperm)
        for (i in 1:nperm) {
          x1 <- runif(n, 0, 1)
          rp[i] <- calculateSIT(x1, y, c = c)
        }
        return(list(
          sitcor = sitcor,
          sd = sd(rp),
          pval = mean(rp > sitcor)
        ))
      }
    }
  }
  # If only sitcor is desired, return value for sit:
  return(sitcor)
}
