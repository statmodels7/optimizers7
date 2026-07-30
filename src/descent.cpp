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
#include "bounded.h"
#include "loop_support.h"

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
                       Rcpp::List bounds) {

  std::unique_ptr<Objective> inner(make_objective(spec));

  // A box is removed rather than enforced: the objective is wrapped so that the
  // algorithm sees an unconstrained problem in eta and never learns a box was
  // involved.
  const bool bounded = bounds.size() > 0;
  std::vector<Coord> coords;
  std::unique_ptr<BoundedObjective> wrap;
  Objective* obj_ptr = inner.get();
  if (bounded) {
    coords = make_coords(bounds);
    wrap.reset(new BoundedObjective(inner.get(), coords));
    obj_ptr = wrap.get();
    par = to_eta(coords, par);
  }
  Objective* obj = obj_ptr;
  std::unique_ptr<Direction> dir(make_direction(method, par.n_elem));
  const LineSearchSpec ls = LineSearchSpec::from_list(line_search);

  arma::vec x = par;
  double f = obj->value(x);
  if (!std::isfinite(f)) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  Trace tr;

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
      if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
      break;
    }

    const arma::vec d = dir->compute(*obj, x, g, f, guard);

    LineSearchResult res = run_line_search(ls, *obj, x, f, g, d, step);
    if (!res.ok) {
      note = "the line search found no acceptable step";
      stopped_by = "failed";
      if (res.guard != "none") guard = res.guard;
      if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
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

    if (keep_trace) tr.push(it, f, gnorm, res.step, guard);

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << f
                  << std::setw(14) << arma::abs(g).max()
                  << std::setw(12) << res.step
                  << "  " << guard << "\n";
    }

    if (ask_criterion(crit_fn, criterion, it, f, f_old, x, x_old, g, true, true)) {
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

  // Reported on the scale the user thinks in. `par` and `gradient` must refer
  // to the same thing, so the gradient is mapped back too rather than left as
  // the eta-scale one the criterion happened to use.
  arma::vec par_out = bounded ? wrap->report_par(x) : x;
  arma::vec g_out   = bounded ? wrap->report_grad(x, g) : g;

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(par_out),
    Rcpp::Named("value")      = f,
    Rcpp::Named("gradient")   = as_r_vector(g_out),
    Rcpp::Named("iterations") = it,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("stopped_by") = stopped_by,
    Rcpp::Named("message")    = note,
    Rcpp::Named("n_value")    = obj->n_value,
    Rcpp::Named("n_grad")     = obj->n_grad,
    Rcpp::Named("n_hess")     = obj->n_hess,
    Rcpp::Named("trace")      = tr.as_r(keep_trace)
  );
}
