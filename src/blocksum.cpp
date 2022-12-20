#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Compute the block-wise sum of a vector.
//'
//'
//' @param r An integer vector
//' @param c The number of observations in each block
//' @return The function returns the block sum of the vector.
//' @export
// [[Rcpp::export]]
double blocksum(arma::ivec r, int c)
{
  int n = r.n_elem;
  int h = floor(n / c);
  arma::imat M = arma::zeros<arma::imat>(c, c);
  int idx_begin = 0;
  int idx_end = 0;
  double s = 0;
  for (int ih = 0; ih < h; ih++)
  {
    idx_begin = c * ih;
    idx_end = c * (ih + 1) - 1;
    M = arma::repmat(r.subvec(idx_begin, idx_end), 1, c);
    s += arma::accu(arma::abs(arma::trimatl(M) - arma::trimatl(M.t())));
  }
  idx_begin = c * h;
  idx_end = n - 1;
  if (idx_begin <= idx_end)
  {
    arma::imat tmp = arma::repmat(r.subvec(idx_begin, idx_end), 1, idx_end - idx_begin + 1);
    s += arma::accu(arma::abs(arma::trimatl(tmp) - arma::trimatl(tmp.t())));
  }

  Rcpp::checkUserInterrupt();
  return s;
}
