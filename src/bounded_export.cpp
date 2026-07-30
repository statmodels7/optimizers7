// The bound transforms, exposed so that the test suite can hold them against
// linkfunctions7. They are not part of the interface and exist for that check:
// the transforms are written out in C++ rather than called through R, so the
// only thing standing between them and silent drift is an assertion that they
// still agree with the package they were copied from.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "bounded.h"

using namespace optimizers7;

// [[Rcpp::export]]
Rcpp::List bounded_transform_cpp(Rcpp::NumericVector b, arma::vec eta) {
  Rcpp::List bl = Rcpp::List::create(b);
  std::vector<Coord> cs = make_coords(bl);
  const Coord c = cs[0];

  const arma::uword n = eta.n_elem;
  arma::vec h(n), d1(n), d2(n);
  for (arma::uword i = 0; i < n; ++i) {
    h[i]  = h_value(c, eta[i]);
    d1[i] = h_deriv1(c, eta[i]);
    d2[i] = h_deriv2(c, eta[i]);
  }
  return Rcpp::List::create(
    Rcpp::Named("h")  = as_r_vector(h),
    Rcpp::Named("d1") = as_r_vector(d1),
    Rcpp::Named("d2") = as_r_vector(d2)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector bounded_forward_cpp(Rcpp::NumericVector b, arma::vec theta) {
  Rcpp::List bl = Rcpp::List::create(b);
  std::vector<Coord> cs = make_coords(bl);
  const Coord c = cs[0];
  arma::vec e(theta.n_elem);
  for (arma::uword i = 0; i < theta.n_elem; ++i) e[i] = g_value(c, theta[i]);
  return as_r_vector(e);
}
