// The R entry point of the descent loop. The loop itself is in
// inst/include/optimizers7/descent_loop.h, where other packages reach it.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <optimizers7/descent_loop.h>
#include <optimizers7/rules.h>

using namespace optimizers7;

// [[Rcpp::export]]
Rcpp::List descent_run(Rcpp::List spec,
                       arma::vec par,
                       SEXP criterion,
                       SEXP crit_fn,
                       Rcpp::List method,
                       Rcpp::List line_search,
                       int maxit,
                       int max_eval,
                       bool verbose,
                       int refresh,
                       bool keep_trace,
                       double step,
                       Rcpp::List bounds,
                       SEXP rule_spec = R_NilValue) {

  std::unique_ptr<Objective> inner(make_objective(spec));

  // A box is removed rather than enforced: the objective is wrapped so that the
  // algorithm sees an unconstrained problem in eta and never learns a box was
  // involved.
  std::vector<Coord> coords;
  if (bounds.size() > 0) coords = make_coords(bounds);
  std::unique_ptr<Direction> dir(make_direction(method, par.n_elem));
  const LineSearchSpec ls = LineSearchSpec::from_list(line_search);
  // The rule in its compiled form where criterion_spec() gave one, and the R
  // rule asked through crit_met() otherwise.
  RStopRule r_rule(crit_fn, criterion);
  std::unique_ptr<CompiledStopRule> c_rule;
  if (rule_spec != R_NilValue) {
    c_rule.reset(new CompiledStopRule(make_compiled_rule(Rcpp::List(rule_spec))));
  }
  StopRule& rule = c_rule ? static_cast<StopRule&>(*c_rule) : static_cast<StopRule&>(r_rule);

  DescentOut r = descent_loop(*inner, par, rule, *dir, ls, maxit, max_eval,
                              verbose, refresh, keep_trace, step, coords);

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(r.par),
    Rcpp::Named("value")      = r.value,
    Rcpp::Named("gradient")   = as_r_vector(r.gradient),
    Rcpp::Named("iterations") = r.iterations,
    Rcpp::Named("converged")  = r.converged,
    Rcpp::Named("stopped_by") = r.stopped_by,
    Rcpp::Named("message")    = r.message,
    Rcpp::Named("n_value")    = r.n_value,
    Rcpp::Named("n_grad")     = r.n_grad,
    Rcpp::Named("n_hess")     = r.n_hess,
    Rcpp::Named("trace")      = r.trace.as_r(keep_trace)
  );
}
