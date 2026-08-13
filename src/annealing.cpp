// Simulated annealing with an adaptive step.
//
// The loop is Corana, Marchesi, Martini and Ridella (1987): the parameters are
// moved ONE COORDINATE AT A TIME, and the step of each coordinate is adjusted
// every few sweeps to hold that coordinate's acceptance rate near a target.
// That adaptation is the whole reason the method is usable on a statistical
// objective, where the coordinates are on scales orders of magnitude apart -- a
// log-scale, a coefficient and a partial autocorrelation on a rhobit chart all
// sit in the same vector, and a single step length is wrong for all three.
//
// The proposal itself is either uniform on the coordinate's step (Corana) or
// Cauchy (Szu and Hartley 1987, the q = 2 member of the Tsallis family), whose
// heavy tail lets a run leave a basin the uniform proposal would have to walk
// out of.
//
// What is NOT here, deliberately: the general Tsallis visiting distribution at
// arbitrary q, whose generator would have to be transcribed and could not be
// validated against anything the package already has. The Cauchy case can be,
// and is.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "objective.h"
#include "bounded.h"
#include "loop_support.h"

using namespace optimizers7;

namespace {

// +Inf rather than NaN, so a proposal into a region where the objective does
// not exist loses every comparison instead of poisoning them.
inline double safe_value(Objective& obj, const arma::vec& x) {
  const double v = obj.value(x);
  return std::isfinite(v) ? v : R_PosInf;
}

// The Metropolis rule. Written so that a non-finite proposal and a frozen
// temperature both reject rather than divide by zero.
inline bool accept(double f_new, double f_old, double temp) {
  if (f_new <= f_old) return true;
  if (!std::isfinite(f_new) || temp <= 0.0) return false;
  return R::unif_rand() < std::exp(-(f_new - f_old) / temp);
}

} // namespace


// [[Rcpp::export]]
Rcpp::List sa_run(Rcpp::List spec,
                  arma::vec par,
                  SEXP criterion,
                  SEXP crit_fn,
                  bool cauchy,
                  double t0,
                  double cooling,
                  int cycles,
                  int steps,
                  double step,
                  double target_accept,
                  double adjust,
                  int n_eps,
                  int maxit,
                  int max_eval,
                  bool verbose,
                  int refresh,
                  bool keep_trace,
                  Rcpp::List bounds) {

  std::unique_ptr<Objective> inner(make_objective(spec));

  const bool bounded = bounds.size() > 0;
  std::vector<Coord> coords;
  std::unique_ptr<BoundedObjective> wrap;
  Objective* obj = inner.get();
  if (bounded) {
    coords = make_coords(bounds);
    wrap.reset(new BoundedObjective(inner.get(), coords));
    obj = wrap.get();
    par = to_eta(coords, par);
  }

  const arma::uword p = par.n_elem;
  arma::vec x = par;
  double f = safe_value(*obj, x);
  if (!std::isfinite(f)) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  // The best point SEEN, which is what the run returns. An annealing run wanders
  // by construction and its last iterate is a sample from a distribution, not an
  // answer; reporting it would report noise.
  arma::vec best_x = x;
  double best_f = f;

  // One step length per coordinate, started relative to where that coordinate
  // is so that a parameter measured in thousands is not proposed in tenths.
  arma::vec v(p);
  for (arma::uword j = 0; j < p; ++j) {
    v[j] = step * std::max(1.0, std::abs(x[j]));
  }

  // The initial temperature is CALIBRATED unless the caller set it. A fixed
  // number cannot serve an objective whose scale is unknown: at 1e6 every
  // proposal is accepted and the run is a random walk, at 1e-6 none is and it is
  // a bad local search. Sampling the objective's own variation and asking for a
  // stated initial acceptance rate is the standard repair, and it costs a few
  // evaluations that are reported like any other.
  std::string note = "";
  if (!(t0 > 0.0)) {
    double total = 0.0;
    int m = 0;
    const int n_calib = std::max(10, static_cast<int>(10 * p));
    for (int k = 0; k < n_calib; ++k) {
      const arma::uword j = static_cast<arma::uword>(
        std::min<double>(p - 1, std::floor(R::unif_rand() * p)));
      arma::vec trial = x;
      trial[j] += v[j] * (2.0 * R::unif_rand() - 1.0);
      const double ft = safe_value(*obj, trial);
      if (std::isfinite(ft)) { total += std::abs(ft - f); ++m; }
    }
    const double mean_jump = (m > 0) ? total / m : 0.0;
    // exp(-d/T) = 0.8 at the average uphill move
    t0 = (mean_jump > 0.0) ? mean_jump / (-std::log(0.8)) : 1.0;
    note = "initial temperature calibrated to " + std::to_string(t0);
  }
  double temp = t0;

  Trace tr;
  bool converged = false;
  std::string stopped_by = "maxit";
  int it = 0;

  // Corana's own termination rule, and BOTH of its halves. The first is that
  // the value the WALK ends each temperature level at has not moved over the
  // last few levels; the second is that it sits at the best point found. The
  // second half is not decoration: measured without it, the rule reads the
  // best-so-far, which stops improving long before the walk settles, so a run
  // stopped at the seventh level whatever its budget and returned a point two
  // tenths from the answer. It is reported as the stationarity measure, so
  // crit_stationary() IS that rule and no second convention is invented.
  // Before enough levels have run there is nothing to compare against and the
  // measure is infinite, which no tolerance can pass.
  std::vector<double> recent;
  double stat = R_PosInf;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective  stationarity        step  safeguard\n";
  }

  double f_prev = f;
  arma::vec x_prev = x;
  bool have_prev = false;
  arma::vec n_acc(p);

  for (it = 1; it <= maxit; ++it) {
    Rcpp::checkUserInterrupt();

    if (have_prev &&
        ask_criterion(crit_fn, criterion, it - 1, best_f, f_prev, best_x, x_prev,
                      arma::vec(), false, true, stat, true)) {
      converged = true;
      stopped_by = "criterion";
      it = it - 1;
      break;
    }
    f_prev = best_f;
    x_prev = best_x;
    have_prev = true;

    bool out_of_budget = false;
    for (int cyc = 0; cyc < cycles && !out_of_budget; ++cyc) {
      n_acc.zeros();
      for (int m = 0; m < steps && !out_of_budget; ++m) {
        for (arma::uword j = 0; j < p; ++j) {
          const double dev = cauchy ? R::rcauchy(0.0, 1.0)
                                    : (2.0 * R::unif_rand() - 1.0);
          arma::vec trial = x;
          trial[j] += v[j] * dev;
          const double ft = safe_value(*obj, trial);
          if (accept(ft, f, temp)) {
            x = trial;
            f = ft;
            n_acc[j] += 1.0;
            if (f < best_f) { best_f = f; best_x = x; }
          }
          if (obj->n_value >= max_eval) {
            out_of_budget = true;
            stopped_by = "max_eval";
            break;
          }
        }
      }

      // Corana's step adjustment: hold each coordinate's acceptance rate inside
      // a band around the target. Too many acceptances means the steps are too
      // short to explore, too few that they are too long to land anywhere.
      const double lo = target_accept - 0.1;
      const double hi = target_accept + 0.1;
      for (arma::uword j = 0; j < p; ++j) {
        const double ratio = n_acc[j] / static_cast<double>(steps);
        if (ratio > hi) {
          v[j] *= 1.0 + adjust * (ratio - hi) / (1.0 - hi);
        } else if (ratio < lo) {
          v[j] /= 1.0 + adjust * (lo - ratio) / lo;
        }
        if (!std::isfinite(v[j]) || v[j] <= 0.0) v[j] = step;
      }
    }

    recent.push_back(f);
    if (static_cast<int>(recent.size()) > n_eps) {
      recent.erase(recent.begin());
      stat = std::abs(f - best_f);
      for (size_t k = 0; k < recent.size(); ++k) {
        stat = std::max(stat, std::abs(f - recent[k]));
      }
    }

    if (keep_trace) {
      tr.push_stat(it, best_f, stat, arma::mean(v),
                   out_of_budget ? "budget" : "none");
    }
    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << "  " << it << "  " << best_f << "  " << stat << "  "
                  << arma::mean(v) << "\n";
    }

    if (out_of_budget) break;
    temp *= cooling;
  }
  if (it > maxit) it = maxit;

  // The run ends where the best point is, not where the walk happened to be.
  // Whether it CONVERGED is a separate question and is answered by the rule
  // above, never by the schedule having finished: a temperature reaching zero
  // says the budget is spent, not that the point is any good.
  if (verbose) {
    Rcpp::Rcout << "  stopped after " << it << " temperature levels ("
                << stopped_by << "), best objective " << best_f << "\n";
  }

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(bounded ? wrap->report_par(best_x)
                                                    : best_x),
    Rcpp::Named("value")      = best_f,
    Rcpp::Named("gradient")   = R_NilValue,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("stopped_by") = stopped_by,
    Rcpp::Named("message")    = note,
    Rcpp::Named("n_value")    = inner->n_value,
    Rcpp::Named("n_grad")     = inner->n_grad,
    Rcpp::Named("n_hess")     = inner->n_hess,
    Rcpp::Named("trace")      = tr.as_r(keep_trace)
  );
}
