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
#include <deque>

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

  // The window a nonmonotone search compares against. Empty when the search is
  // monotone, in which case f_ref below is just f and nothing changes.
  const int mem = ls.memory;
  std::deque<double> recent;
  if (mem > 0) recent.push_back(f);

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

    // Grippo, Lampariello and Lucidi: the reference is the LARGEST of the last
    // mem + 1 values rather than the current one. A step that makes the
    // objective worse than where it stands is then admissible, provided it
    // improves on the worst of recent memory -- which is exactly what a method
    // whose efficiency comes from such steps needs, and exactly what an Armijo
    // condition forbids.
    double f_ref = f;
    if (mem > 0) for (double v : recent) f_ref = std::max(f_ref, v);

    LineSearchResult res = run_line_search(ls, *obj, x, f, g, d, step, f_ref);
    if (!res.ok) {
      // A search that cannot improve on the current point is a failure only if
      // the point is not already an answer, and it routinely is: a caller that
      // starts from a closed-form estimate hands the optimizer the solution.
      // The rule is asked with no previous iterate, so a criterion reading a
      // change in the objective returns FALSE and only the state at x -- its
      // gradient, its stationarity -- can declare the run finished.
      if (ask_criterion(crit_fn, criterion, it, f, f, x, x, g, true, false)) {
        converged = true;
        stopped_by = "criterion";
        if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
        break;
      }
      // A search that stopped because the decrease it must verify fell below
      // what the objective can RESOLVE has not failed: it has established that
      // no step along a descent direction improves on this point by anything
      // this objective can show. That is convergence to the accuracy the
      // objective has, and it is reported in those words rather than as a
      // criterion being met, so a reader is told the weaker true thing and not
      // the stronger false one. It is reachable only where the caller declared
      // a resolution; the default of zero leaves this branch dead.
      if (res.unresolvable) {
        converged = true;
        stopped_by = "resolution";
        if (res.guard != "none") guard = res.guard;
        if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
        break;
      }
      note = "the line search found no acceptable step";
      stopped_by = "failed";
      if (res.guard != "none") guard = res.guard;
      if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
      break;
    }
    if (guard == "none" && res.guard != "none") guard = res.guard;

    // A step that moves the iterate by no more than rounding has not moved it,
    // and is handled as the failed search above. Where the objective is
    // rounding noise (a filter that has exploded: a gradient of 1e27 and a
    // Hessian of 1e55, so that the repaired direction is 1e-20 long) such a
    // step changes the objective by the noise and passes Armijo whenever the
    // noise is downwards, and a rule reading the change in the parameters or
    // in the objective then fired: converged at |g| = 4e27 on one platform,
    // not on another, at the same point.
    const double dx_max = arma::abs(res.x_new - x).max();
    const double x_max = x.is_empty() ? 0.0 : arma::abs(x).max();
    if (dx_max <= 4.0 * arma::datum::eps * std::max(1.0, x_max)) {
      if (ask_criterion(crit_fn, criterion, it, f, f, x, x, g, true, false)) {
        converged = true;
        stopped_by = "criterion";
        if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
        break;
      }
      note = "the step did not move the iterate beyond rounding";
      stopped_by = "failed";
      guard = "step below rounding";
      if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
      break;
    }

    const arma::vec x_old = x;
    const arma::vec g_old = g;
    const double f_old = f;

    x = res.x_new;
    f = res.f_new;

    if (mem > 0) {
      recent.push_back(f);
      while (static_cast<int>(recent.size()) > mem + 1) recent.pop_front();
    }
    // The Wolfe search already evaluated the gradient at the accepted point.
    g = res.g_new_valid ? res.g_new : obj->grad(x);

    // The secant pair, for the methods that remember. Passing guard through
    // lets a skipped update be reported at the iteration where it happened.
    dir->update(x - x_old, g - g_old, guard);

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << f
                  << std::setw(14) << arma::abs(g).max()
                  << std::setw(12) << res.step
                  << "  " << guard << "\n";
    }

    if (ask_criterion(crit_fn, criterion, it, f, f_old, x, x_old, g, true, true)) {
      // A rule on the change in the parameters or in the objective does not
      // declare convergence on a step along which the first-order change at
      // the new point, |g' dx|, exceeds max(1, |f|): the objective does not
      // resolve such a step, which passed Armijo on noise. The step is undone
      // and handled as the failed search above. Measured, the ratio
      // |g' dx| / max(1, |f|) at a firing that the gradient rule alone would
      // not make is at most 1.8e-5 over the suites of optimizers7,
      // distributions7 and statmodels7, and at least 2e4 on objectives that
      // are rounding noise (steps of up to 3e-11, past the test above).
      const double moved = std::abs(arma::dot(g, x - x_old));
      if (moved > std::max(1.0, std::abs(f))) {
        x = x_old;
        f = f_old;
        g = g_old;
        if (ask_criterion(crit_fn, criterion, it, f, f, x, x, g, true, false)) {
          converged = true;
          stopped_by = "criterion";
          if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
          break;
        }
        note = "the step that met the stopping rule is not resolved by the objective";
        stopped_by = "failed";
        guard = "step not resolved";
        if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
        break;
      }
      // A rule on the change does not close the run at a point where the
      // direction can tell from the curvature that the point is not a minimum
      // (newton(): a Hessian without a Cholesky factor, or a Newton decrement
      // above stop_decrement), nor, for a direction without a Hessian, where
      // the step just taken still predicts a decrease |g'dx| / 2 above
      // stop_decrement at the new point. Only a gradient rule can then close
      // the run, and the rule is asked for that at the new point with no
      // previous iterate.
      std::string why;
      const double step_dec = dir->step_stop_decrement();
      if ((dir->vets_change_stop() || std::isfinite(step_dec)) &&
          !ask_criterion(crit_fn, criterion, it, f, f, x, x, g, true, false)) {
        bool refused = false;
        if (dir->vets_change_stop() &&
            !dir->accepts_change_stop(*obj, x, g, why)) {
          refused = true;
        } else if (std::isfinite(step_dec) && 0.5 * moved > step_dec) {
          std::ostringstream os;
          os << "the last step predicts a decrease of " << 0.5 * moved
             << ", above stop_decrement = " << step_dec;
          why = os.str();
          refused = true;
        }
        if (refused) {
          note = "the stopping rule was met on a change at a point that is not a minimum: " + why;
          stopped_by = "failed";
          guard = "change stop refused";
          if (keep_trace) tr.push(it, f, gnorm, res.step, guard);
          break;
        }
      }
      if (keep_trace) tr.push(it, f, gnorm, res.step, guard);
      converged = true;
      stopped_by = "criterion";
      break;
    }

    if (keep_trace) tr.push(it, f, gnorm, res.step, guard);

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
