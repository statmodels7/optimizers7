// The descent loop, written once for every method that has a direction.
//
// This replaces the per-algorithm driver the skeleton started with. Steepest
// descent, Newton, BFGS and L-BFGS all run through here; what differs is the
// Direction they carry, and nothing else. Adding a fifth is a Direction and a
// constructor.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "objective.h"
#include "line_search.h"
#include "direction.h"

using namespace optimizers7;

static bool criterion_met(SEXP crit_fn, SEXP criterion, int iter,
                          double f_new, double f_old,
                          const arma::vec& x_new, const arma::vec& x_old,
                          const arma::vec& g, bool have_g, bool have_old) {
  // Built with the optional slots empty and filled in afterwards: a ternary
  // whose branches are a NumericVector and R_NilValue has no common type.
  Rcpp::List state = Rcpp::List::create(
    Rcpp::Named("iter")     = iter,
    Rcpp::Named("f_new")    = f_new,
    Rcpp::Named("f_old")    = R_NilValue,
    Rcpp::Named("x_new")    = as_r_vector(x_new),
    Rcpp::Named("x_old")    = R_NilValue,
    Rcpp::Named("gradient") = R_NilValue
  );
  if (have_old) {
    state["f_old"] = f_old;
    state["x_old"] = as_r_vector(x_old);
  }
  if (have_g) state["gradient"] = as_r_vector(g);

  Rcpp::LogicalVector out = Rcpp::Function(crit_fn)(criterion, state);
  return out.size() == 1 && out[0] == TRUE;
}

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
                       double step) {

  std::unique_ptr<Objective> obj(make_objective(spec));
  std::unique_ptr<Direction> dir(make_direction(method, par.n_elem));
  const LineSearchSpec ls = LineSearchSpec::from_list(line_search);

  arma::vec x = par;
  double f = obj->value(x);
  if (!std::isfinite(f)) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  std::vector<int>    tr_iter;
  std::vector<double> tr_value, tr_gnorm, tr_step;
  std::vector<std::string> tr_guard;

  bool converged = false;
  std::string stopped_by = "maxit";
  std::string note = "";
  int it = 0;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective     |grad|max        step  safeguard\n";
  }

  arma::vec g = obj->grad(x);

  for (it = 1; it <= maxit; ++it) {
    // A C++ loop that ignores interrupts is a session that cannot be stopped.
    Rcpp::checkUserInterrupt();

    const double gnorm = g.is_empty() ? 0.0 : arma::abs(g).max();
    std::string guard = "none";

    if (!g.is_finite()) {
      note = "gradient not finite; stopped";
      stopped_by = "failed";
      guard = "gradient non-finite";
      if (keep_trace) {
        tr_iter.push_back(it); tr_value.push_back(f);
        tr_gnorm.push_back(gnorm); tr_step.push_back(0.0);
        tr_guard.push_back(guard);
      }
      break;
    }

    const arma::vec d = dir->compute(*obj, x, g, f, guard);

    LineSearchResult res = run_line_search(ls, *obj, x, f, g, d, step);
    if (!res.ok) {
      note = "the line search found no acceptable step";
      stopped_by = "failed";
      if (res.guard != "none") guard = res.guard;
      if (keep_trace) {
        tr_iter.push_back(it); tr_value.push_back(f);
        tr_gnorm.push_back(gnorm); tr_step.push_back(0.0);
        tr_guard.push_back(guard);
      }
      break;
    }
    if (guard == "none" && res.guard != "none") guard = res.guard;

    const arma::vec x_old = x;
    const arma::vec g_old = g;
    const double f_old = f;

    x = res.x_new;
    f = res.f_new;
    // The Wolfe search already evaluated the gradient at the accepted point.
    g = res.g_new_valid ? res.g_new : obj->grad(x);

    // The secant pair, for the methods that remember. Passing guard through
    // lets a skipped update be reported at the iteration where it happened.
    dir->update(x - x_old, g - g_old, guard);

    if (keep_trace) {
      tr_iter.push_back(it); tr_value.push_back(f);
      tr_gnorm.push_back(gnorm); tr_step.push_back(res.step);
      tr_guard.push_back(guard);
    }

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << f
                  << std::setw(14) << arma::abs(g).max()
                  << std::setw(12) << res.step
                  << "  " << guard << "\n";
    }

    if (criterion_met(crit_fn, criterion, it, f, f_old, x, x_old, g, true, true)) {
      converged = true;
      stopped_by = "criterion";
      break;
    }

    if (obj->n_value >= max_eval) {
      note = "evaluation budget exhausted";
      stopped_by = "max_eval";
      break;
    }
  }
  if (it > maxit) it = maxit;

  if (verbose) {
    Rcpp::Rcout << "  stopped after " << it << " iterations ("
                << stopped_by << "), objective " << f << "\n";
  }

  // SEXP, not Rcpp::List: assigning R_NilValue to an Rcpp::List builds an empty
  // list rather than NULL.
  SEXP trace = R_NilValue;
  if (keep_trace) {
    trace = Rcpp::DataFrame::create(
      Rcpp::Named("iteration") = tr_iter,
      Rcpp::Named("value")     = tr_value,
      Rcpp::Named("gnorm")     = tr_gnorm,
      Rcpp::Named("step")      = tr_step,
      Rcpp::Named("safeguard") = tr_guard,
      Rcpp::Named("stringsAsFactors") = false
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(x),
    Rcpp::Named("value")      = f,
    Rcpp::Named("gradient")   = as_r_vector(g),
    Rcpp::Named("iterations") = it,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("stopped_by") = stopped_by,
    Rcpp::Named("message")    = note,
    Rcpp::Named("n_value")    = obj->n_value,
    Rcpp::Named("n_grad")     = obj->n_grad,
    Rcpp::Named("n_hess")     = obj->n_hess,
    Rcpp::Named("trace")      = trace
  );
}
