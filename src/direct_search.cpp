// The two methods that never look at a derivative.
//
// They share a loop shape with nothing else in the package: there is no
// direction, no line search, and no gradient, so there is nothing for
// descent.cpp to give them. What they do share is the reporting, through
// loop_support.h, and the stationarity measure they report in place of a
// gradient norm -- see crit_stationary() on the R side for why that had to be a
// new quantity rather than a reuse of an existing one.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "objective.h"
#include "bounded.h"
#include "loop_support.h"

using namespace optimizers7;

namespace {

// A value that is not finite must never win a comparison. Returning +Inf rather
// than propagating the NaN is what lets a method walk up to the edge of its
// objective's domain, be told there is nothing there, and carry on.
inline double safe_value(Objective& obj, const arma::vec& x) {
  const double v = obj.value(x);
  return std::isfinite(v) ? v : R_PosInf;
}

// How far the simplex is from being flat: |det E| divided by the product of the
// edge lengths, where E holds the edges from the best vertex. It is 1 for a
// right-angled simplex and 0 for one that has collapsed into a lower dimension,
// and it is scale-free, so the same threshold means the same thing however
// large the simplex is.
double simplex_conditioning(const arma::mat& V) {
  const arma::uword p = V.n_rows;
  if (p == 0) return 1.0;
  arma::mat E(p, p);
  double prod = 1.0;
  for (arma::uword j = 0; j < p; ++j) {
    E.col(j) = V.col(j + 1) - V.col(0);
    const double nj = arma::norm(E.col(j), 2);
    if (nj <= 0.0) return 0.0;
    prod *= nj;
    E.col(j) /= nj;
  }
  return std::abs(arma::det(E));
}

// An axis-aligned simplex of a given edge length, used both to start and to
// restart. Restarting with the edge length already reached keeps the scale the
// run has earned instead of throwing it away.
arma::mat build_simplex(const arma::vec& x0, double edge) {
  const arma::uword p = x0.n_elem;
  arma::mat V(p, p + 1);
  V.col(0) = x0;
  for (arma::uword j = 0; j < p; ++j) {
    V.col(j + 1) = x0;
    V(j, j + 1) += edge * std::max(1.0, std::abs(x0[j]));
  }
  return V;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List nelder_mead_run(Rcpp::List spec,
                           arma::vec par,
                           SEXP criterion,
                           SEXP crit_fn,
                           double step,
                           bool adaptive,
                           int max_restarts,
                           double degenerate_tol,
                           Rcpp::NumericMatrix start_simplex,
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
  const double n = static_cast<double>(p);

  // Gao and Han's dimension-dependent coefficients. At p = 2 they reduce
  // EXACTLY to the classical 1, 2, 1/2, 1/2, so switching them on costs nothing
  // on a small problem and rescues a large one, where a fixed expansion of 2
  // makes the simplex chase the wrong direction.
  const double c_ref = 1.0;
  const double c_exp = adaptive ? 1.0 + 2.0 / n : 2.0;
  const double c_con = adaptive ? 0.75 - 1.0 / (2.0 * n) : 0.5;
  const double c_shr = adaptive ? 1.0 - 1.0 / n : 0.5;

  arma::mat V;
  if (start_simplex.nrow() > 0) {
    // Given by the user, one vertex per ROW, which is how anyone writes one
    // down; stored here one per column, which is how the arithmetic wants it.
    arma::mat S = Rcpp::as<arma::mat>(start_simplex);
    V = S.t();
  } else {
    V = build_simplex(par, step);
  }

  arma::vec fv(p + 1);
  for (arma::uword j = 0; j <= p; ++j) fv[j] = safe_value(*obj, V.col(j));
  if (!std::isfinite(fv[0])) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  Trace tr;
  bool converged = false;
  std::string stopped_by = "maxit";
  std::string note = "";
  int it = 0;
  int restarts = 0;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective  stationarity        step  safeguard\n";
  }

  double f_prev = fv.min();
  arma::vec x_prev = V.col(0);
  bool have_prev = false;

  for (it = 1; it <= maxit; ++it) {
    Rcpp::checkUserInterrupt();

    // Sort the vertices, best first. Everything below names them by rank.
    arma::uvec ord = arma::sort_index(fv);
    V = V.cols(ord);
    fv = fv(ord);

    double diameter = 0.0;
    for (arma::uword j = 1; j <= p; ++j) {
      diameter = std::max(diameter, arma::norm(V.col(j) - V.col(0), 2));
    }

    std::string guard = "none";

    // THE McKINNON SAFEGUARD.
    //
    // Nelder-Mead can converge to a point that is not a minimizer at all:
    // McKinnon (1998) exhibited a strictly convex function on which it performs
    // inside contractions for ever, the simplex flattening onto a line through
    // a non-stationary point while every ordinary stopping rule reports
    // success. What goes wrong is visible in the SHAPE of the simplex rather
    // than in any value, so that is what is watched: when the simplex has
    // collapsed towards a lower dimension it is rebuilt, right-angled, at the
    // current best vertex and at the size it had reached.
    if (restarts < max_restarts && diameter > 0.0 &&
        simplex_conditioning(V) < degenerate_tol) {
      const arma::vec best = V.col(0);
      V = build_simplex(best, diameter);
      for (arma::uword j = 0; j <= p; ++j) fv[j] = safe_value(*obj, V.col(j));
      ord = arma::sort_index(fv);
      V = V.cols(ord);
      fv = fv(ord);
      ++restarts;
      guard = "restart";
    }

    if (have_prev &&
        ask_criterion(crit_fn, criterion, it - 1, fv[0], f_prev, V.col(0),
                      x_prev, arma::vec(), false, true, diameter, true)) {
      converged = true;
      stopped_by = "criterion";
      it = it - 1;
      break;
    }
    f_prev = fv[0];
    x_prev = V.col(0);
    have_prev = true;

    // The centroid of every vertex but the worst.
    arma::vec c(p, arma::fill::zeros);
    for (arma::uword j = 0; j < p; ++j) c += V.col(j);
    c /= n;

    const arma::vec worst = V.col(p);
    const arma::vec xr = c + c_ref * (c - worst);
    const double fr = safe_value(*obj, xr);

    if (fr < fv[0]) {
      const arma::vec xe = c + c_exp * (xr - c);
      const double fe = safe_value(*obj, xe);
      if (fe < fr) { V.col(p) = xe; fv[p] = fe; if (guard == "none") guard = "expand"; }
      else         { V.col(p) = xr; fv[p] = fr; if (guard == "none") guard = "reflect"; }
    } else if (fr < fv[p - 1]) {
      V.col(p) = xr; fv[p] = fr;
      if (guard == "none") guard = "reflect";
    } else {
      bool shrink = true;
      if (fr < fv[p]) {
        const arma::vec xc = c + c_con * (xr - c);      // outside contraction
        const double fc = safe_value(*obj, xc);
        if (fc <= fr) { V.col(p) = xc; fv[p] = fc; shrink = false;
                        if (guard == "none") guard = "contract out"; }
      } else {
        const arma::vec xc = c - c_con * (c - worst);   // inside contraction
        const double fc = safe_value(*obj, xc);
        if (fc < fv[p]) { V.col(p) = xc; fv[p] = fc; shrink = false;
                          if (guard == "none") guard = "contract in"; }
      }
      if (shrink) {
        for (arma::uword j = 1; j <= p; ++j) {
          V.col(j) = V.col(0) + c_shr * (V.col(j) - V.col(0));
          fv[j] = safe_value(*obj, V.col(j));
        }
        if (guard == "none") guard = "shrink";
      }
    }

    if (keep_trace) tr.push_stat(it, fv.min(), diameter, diameter, guard);

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << fv.min()
                  << std::setw(14) << diameter
                  << std::setw(12) << diameter
                  << "  " << guard << "\n";
    }

    if (inner->n_value >= max_eval) {
      note = "evaluation budget exhausted";
      stopped_by = "max_eval";
      break;
    }
  }
  if (it > maxit) it = maxit;

  arma::uvec ord = arma::sort_index(fv);
  V = V.cols(ord);
  fv = fv(ord);
  arma::vec x = V.col(0);

  if (restarts > 0) {
    const std::string s = "simplex restarted " + std::to_string(restarts) +
      (restarts == 1 ? " time" : " times");
    note = note.empty() ? s : note + "; " + s;
  }

  if (verbose) {
    Rcpp::Rcout << "  stopped after " << it << " iterations ("
                << stopped_by << "), objective " << fv[0] << "\n";
  }

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(bounded ? wrap->report_par(x) : x),
    Rcpp::Named("value")      = fv[0],
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


// [[Rcpp::export]]
Rcpp::List compass_run(Rcpp::List spec,
                       arma::vec par,
                       SEXP criterion,
                       SEXP crit_fn,
                       double step,
                       bool random_directions,
                       bool opportunistic,
                       double expand,
                       double shrink,
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

  // The poll size, on the parameter scale and started relative to where we are,
  // so that a problem measured in thousands is not polled in steps of a tenth.
  double delta = step * std::max(1.0, arma::abs(x).max());
  const double delta0 = delta;

  Trace tr;
  bool converged = false;
  std::string stopped_by = "maxit";
  std::string note = "";
  int it = 0;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective  stationarity        step  safeguard\n";
  }

  double f_prev = f;
  arma::vec x_prev = x;
  bool have_prev = false;

  for (it = 1; it <= maxit; ++it) {
    Rcpp::checkUserInterrupt();

    if (have_prev &&
        ask_criterion(crit_fn, criterion, it - 1, f, f_prev, x, x_prev,
                      arma::vec(), false, true, delta, true)) {
      converged = true;
      stopped_by = "criterion";
      it = it - 1;
      break;
    }
    f_prev = f;
    x_prev = x;
    have_prev = true;

    // The poll directions: a positive spanning set, so that if x is not
    // stationary at least one of them goes downhill.
    //
    // Coordinate: the 2p signed axes. Cheap and deterministic, and enough for a
    // theorem when f is continuously differentiable.
    //
    // Random: a fresh random orthonormal basis, taken plus and minus. The
    // theorem for a merely LIPSCHITZ f -- which is the case these methods exist
    // for -- needs the directions used over the whole run to be dense in the
    // sphere, and drawing a new basis at every poll gives that with probability
    // one. This is the idea of MADS rather than LTMADS as published; what it
    // shares is the property the proof actually rests on.
    arma::mat D;
    if (random_directions) {
      arma::mat A(p, p);
      for (arma::uword i = 0; i < p; ++i)
        for (arma::uword j = 0; j < p; ++j) A(i, j) = R::norm_rand();
      arma::mat Q, R;
      if (!arma::qr(Q, R, A)) Q = arma::eye(p, p);
      D = Q;
    } else {
      D = arma::eye(p, p);
    }

    bool improved = false;
    arma::vec best_x = x;
    double best_f = f;
    std::string guard = "none";

    for (arma::uword j = 0; j < p && !(improved && opportunistic); ++j) {
      for (int s = 0; s < 2 && !(improved && opportunistic); ++s) {
        const arma::vec trial = x + (s == 0 ? delta : -delta) * D.col(j);
        const double ft = safe_value(*obj, trial);
        if (ft < best_f) { best_f = ft; best_x = trial; improved = true; }
      }
    }

    if (improved) {
      x = best_x;
      f = best_f;
      delta *= expand;
      // The poll size must not run away: with an expansion at every success a
      // long descent would end up polling at a radius unrelated to the problem.
      if (delta > delta0 * 1e4) { delta = delta0 * 1e4; guard = "poll capped"; }
      else guard = "poll expanded";
    } else {
      delta *= shrink;
      guard = "poll shrunk";
      // Below this the poll points are not distinct from x in double precision,
      // so every further iteration would evaluate the same point.
      if (delta < 1e-300) {
        note = "the poll size underflowed";
        stopped_by = "failed";
        if (keep_trace) tr.push_stat(it, f, delta, delta, "poll underflow");
        break;
      }
    }

    if (keep_trace) tr.push_stat(it, f, delta, delta, guard);

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << f
                  << std::setw(14) << delta
                  << std::setw(12) << delta
                  << "  " << guard << "\n";
    }

    if (inner->n_value >= max_eval) {
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

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(bounded ? wrap->report_par(x) : x),
    Rcpp::Named("value")      = f,
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
