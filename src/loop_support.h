#ifndef OPTIMIZERS7_LOOP_SUPPORT_H
#define OPTIMIZERS7_LOOP_SUPPORT_H

// The parts of a run that are the same whatever the algorithm: asking the
// R-side stopping rule, collecting the trace, and assembling the result.
//
// Factored out when Adam arrived, because Adam does not fit the descent loop --
// it has no line search, guarantees no monotone descent, and computes its own
// stochastic gradient -- and duplicating the bookkeeping would have meant two
// places to keep the reporting honest.

#include <RcppArmadillo.h>
#include "objective.h"

namespace optimizers7 {

inline bool ask_criterion(SEXP crit_fn, SEXP criterion, int iter,
                          double f_new, double f_old,
                          const arma::vec& x_new, const arma::vec& x_old,
                          const arma::vec& g, bool have_g, bool have_old,
                          double stationarity = 0.0, bool have_stat = false) {
  // Built with the optional slots empty and filled in afterwards: a ternary
  // whose branches are a NumericVector and R_NilValue has no common type.
  Rcpp::List state = Rcpp::List::create(
    Rcpp::Named("iter")         = iter,
    Rcpp::Named("f_new")        = f_new,
    Rcpp::Named("f_old")        = R_NilValue,
    Rcpp::Named("x_new")        = as_r_vector(x_new),
    Rcpp::Named("x_old")        = R_NilValue,
    Rcpp::Named("gradient")     = R_NilValue,
    Rcpp::Named("stationarity") = R_NilValue
  );
  if (have_old) {
    state["f_old"] = f_old;
    state["x_old"] = as_r_vector(x_old);
  }
  if (have_g) state["gradient"] = as_r_vector(g);
  if (have_stat) state["stationarity"] = stationarity;

  Rcpp::LogicalVector out = Rcpp::Function(crit_fn)(criterion, state);
  return out.size() == 1 && out[0] == TRUE;
}


struct Trace {
  std::vector<int>    iter;
  std::vector<double> value, gnorm, step, stat;
  std::vector<std::string> guard;
  bool has_stat = false;   // set by the derivative-free loops only

  void push(int it, double f, double gn, double s, const std::string& gd) {
    iter.push_back(it); value.push_back(f);
    gnorm.push_back(gn); step.push_back(s); guard.push_back(gd);
    stat.push_back(NA_REAL);
  }

  // The derivative-free methods have no gradient to report and a stationarity
  // measure instead. The column appears only when one of them filled it in,
  // rather than sitting as a row of NAs in every descent trace.
  void push_stat(int it, double f, double st, double s, const std::string& gd) {
    has_stat = true;
    iter.push_back(it); value.push_back(f);
    gnorm.push_back(NA_REAL); step.push_back(s); guard.push_back(gd);
    stat.push_back(st);
  }

  // SEXP, not Rcpp::List: assigning R_NilValue to an Rcpp::List builds an empty
  // list rather than NULL, and the R side tests for NULL.
  SEXP as_r(bool keep) const {
    if (!keep) return R_NilValue;
    if (has_stat) {
      return Rcpp::DataFrame::create(
        Rcpp::Named("iteration")    = iter,
        Rcpp::Named("value")        = value,
        Rcpp::Named("stationarity") = stat,
        Rcpp::Named("step")         = step,
        Rcpp::Named("safeguard")    = guard,
        Rcpp::Named("stringsAsFactors") = false
      );
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("iteration") = iter,
      Rcpp::Named("value")     = value,
      Rcpp::Named("gnorm")     = gnorm,
      Rcpp::Named("step")      = step,
      Rcpp::Named("safeguard") = guard,
      Rcpp::Named("stringsAsFactors") = false
    );
  }
};

} // namespace optimizers7

#endif
