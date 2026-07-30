#ifndef OPTIMIZERS7_LINE_SEARCH_H
#define OPTIMIZERS7_LINE_SEARCH_H

// The line search, written once for every descent method in the package.
//
// A descent method produces a DIRECTION; how far to go along it is a separate
// question with its own theory, and answering it badly is how a correct
// direction turns into a run that oscillates, stalls, or steps into a region
// where the objective is not defined. Newton, BFGS and L-BFGS will all call
// this; only the direction they hand it differs.

#include <RcppArmadillo.h>
#include "objective.h"

namespace optimizers7 {

struct LineSearchResult {
  double step  = 0.0;      // the accepted step length
  double f_new = 0.0;
  arma::vec x_new;
  arma::vec g_new;         // gradient at x_new; Wolfe needs it, and the caller
  bool g_new_valid = false;// can reuse it instead of recomputing
  bool ok = false;
  int evaluations = 0;
  std::string guard = "none";
};


// Common preamble: the directional derivative, and the check that the direction
// is one of descent at all.
//
// A method whose direction points uphill has already gone wrong -- a Newton step
// through an indefinite Hessian is the usual way -- and no line search can
// rescue it, because every step along that direction increases the objective.
// Detecting it here means the caller can fall back to steepest descent instead
// of grinding through thirty shrinkages before giving up.
inline double directional_derivative(const arma::vec& g, const arma::vec& d) {
  return arma::dot(g, d);
}


// --- Armijo backtracking ----------------------------------------------------
//
// Halve until sufficient decrease. Cheap, needs no gradient at the trial points,
// and enough for a method that only needs to make progress.
//
// The condition is f(x + s d) <= f(x) + c1 s g'd, NOT merely f_new <= f. The
// weak version accepts a step that leaves the objective exactly unchanged, and
// on a quadratic with step 1 the update reflects through the minimum and does
// exactly that: the iterate then oscillates forever while a criterion watching
// the objective sees no change and reports convergence at a point that is not a
// minimum. Sufficient decrease is what forbids it.
inline LineSearchResult armijo_search(Objective& obj,
                                      const arma::vec& x, double f,
                                      const arma::vec& g, const arma::vec& d,
                                      double step0, double c1, double shrink,
                                      int max_step) {
  LineSearchResult out;
  const double dg = directional_derivative(g, d);

  if (!(dg < 0.0)) {
    out.guard = "not a descent direction";
    return out;
  }

  double s = step0;
  for (int k = 0; k < max_step; ++k) {
    arma::vec xt = x + s * d;
    const double ft = obj.value(xt);
    ++out.evaluations;

    if (std::isfinite(ft) && ft <= f + c1 * s * dg) {
      out.step = s;
      out.f_new = ft;
      out.x_new = xt;
      out.ok = true;
      if (k > 0) out.guard = "step shortened";
      return out;
    }
    // A non-finite trial is a refusal, not a failure: the optimiser is entitled
    // to probe a point where the objective does not exist, and the line search
    // simply does not go there.
    if (!std::isfinite(ft)) out.guard = "objective non-finite";
    s *= shrink;
  }

  if (out.guard == "none") out.guard = "no acceptable step";
  return out;
}


// --- strong Wolfe -----------------------------------------------------------
//
// Sufficient decrease AND a curvature condition:
//
//     f(x + s d) <= f(x) + c1 s g'd                     (Armijo)
//     |g(x + s d)'d| <= c2 |g'd|                        (strong curvature)
//
// The second is what Armijo alone cannot give and what BFGS needs: it rules out
// a step so short that the gradient has barely changed, which is precisely the
// case where the secant pair (s, y) carries no curvature information and the
// BFGS update either does nothing useful or has to be skipped. A quasi-Newton
// method run on Armijo alone loses its convergence guarantee for exactly this
// reason, so this is not a refinement but a requirement.
//
// Nocedal and Wright's bracketing-and-zoom scheme. Bisection inside the zoom
// rather than interpolation: a few more evaluations, and it cannot be defeated
// by a badly shaped interval.
inline LineSearchResult wolfe_search(Objective& obj,
                                     const arma::vec& x, double f,
                                     const arma::vec& g, const arma::vec& d,
                                     double step0, double c1, double c2,
                                     int max_step) {
  LineSearchResult out;
  const double dphi0 = directional_derivative(g, d);

  if (!(dphi0 < 0.0)) {
    out.guard = "not a descent direction";
    return out;
  }

  auto phi = [&](double a, double& fa, arma::vec& xa) {
    xa = x + a * d;
    fa = obj.value(xa);
    ++out.evaluations;
  };

  double a_prev = 0.0, f_prev = f;
  double a = step0;
  arma::vec xa, ga;
  double fa = f;

  bool bracketed = false;
  double lo = 0.0, hi = 0.0, f_lo = f;

  for (int k = 0; k < max_step; ++k) {
    phi(a, fa, xa);

    if (!std::isfinite(fa)) {
      // Shrink towards the last point known to be finite rather than treating
      // it as a bracket: there is nothing to interpolate between when one end
      // is not a number.
      out.guard = "objective non-finite";
      a = 0.5 * (a_prev + a);
      if (a <= 0.0) break;
      continue;
    }

    if (fa > f + c1 * a * dphi0 || (k > 0 && fa >= f_prev)) {
      lo = a_prev; hi = a; f_lo = f_prev; bracketed = true;
      break;
    }

    ga = obj.grad(xa);
    const double dphi = directional_derivative(ga, d);

    if (std::abs(dphi) <= -c2 * dphi0) {
      out.step = a; out.f_new = fa; out.x_new = xa;
      out.g_new = ga; out.g_new_valid = true; out.ok = true;
      if (k > 0) out.guard = "step adjusted";
      return out;
    }

    if (dphi >= 0.0) {
      lo = a; hi = a_prev; f_lo = fa; bracketed = true;
      break;
    }

    a_prev = a; f_prev = fa;
    a *= 2.0;
  }

  if (!bracketed) {
    if (out.guard == "none") out.guard = "no acceptable step";
    return out;
  }

  // zoom: shrink the bracket until a point satisfies both conditions
  for (int k = 0; k < max_step; ++k) {
    const double aj = 0.5 * (lo + hi);
    if (!(std::abs(hi - lo) > 0.0)) break;

    double fj;
    arma::vec xj;
    phi(aj, fj, xj);

    if (!std::isfinite(fj) || fj > f + c1 * aj * dphi0 || fj >= f_lo) {
      hi = aj;
      if (!std::isfinite(fj)) out.guard = "objective non-finite";
      continue;
    }

    arma::vec gj = obj.grad(xj);
    const double dphij = directional_derivative(gj, d);

    if (std::abs(dphij) <= -c2 * dphi0) {
      out.step = aj; out.f_new = fj; out.x_new = xj;
      out.g_new = gj; out.g_new_valid = true; out.ok = true;
      out.guard = (out.guard == "none") ? "step adjusted" : out.guard;
      return out;
    }

    if (dphij * (hi - lo) >= 0.0) hi = lo;
    lo = aj; f_lo = fj;
  }

  // The bracket collapsed without meeting the curvature condition. The Armijo
  // end of it is still a decrease, so it is returned rather than failing the
  // iteration -- but it is reported, because a quasi-Newton method may want to
  // skip its update when the curvature condition was not met.
  if (lo > 0.0 && f_lo < f) {
    arma::vec xl = x + lo * d;
    out.step = lo; out.f_new = f_lo; out.x_new = xl;
    out.ok = true;
    out.guard = "curvature condition not met";
    return out;
  }

  if (out.guard == "none") out.guard = "no acceptable step";
  return out;
}


// The settings, as they arrive from the R side.
struct LineSearchSpec {
  std::string type = "armijo";
  double c1 = 1e-4;
  double c2 = 0.9;
  double shrink = 0.5;
  int max_step = 30;

  static LineSearchSpec from_list(Rcpp::List ls) {
    LineSearchSpec s;
    s.type     = Rcpp::as<std::string>(ls["type"]);
    s.c1       = Rcpp::as<double>(ls["c1"]);
    s.c2       = Rcpp::as<double>(ls["c2"]);
    s.shrink   = Rcpp::as<double>(ls["shrink"]);
    s.max_step = Rcpp::as<int>(ls["max_step"]);
    return s;
  }
};

inline LineSearchResult run_line_search(const LineSearchSpec& spec,
                                        Objective& obj,
                                        const arma::vec& x, double f,
                                        const arma::vec& g, const arma::vec& d,
                                        double step0) {
  if (spec.type == "wolfe") {
    return wolfe_search(obj, x, f, g, d, step0, spec.c1, spec.c2, spec.max_step);
  }
  return armijo_search(obj, x, f, g, d, step0, spec.c1, spec.shrink,
                       spec.max_step);
}

} // namespace optimizers7

#endif
