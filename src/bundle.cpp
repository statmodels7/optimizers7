// The proximal bundle method: the one algorithm here built for a function with
// a kink in it.
//
// WHY NEITHER OF THE OTHER TWO WILL DO.
//
// Nelder-Mead and compass search never look at a derivative, so a kink does not
// mislead them -- but neither of them uses any information beyond an ordering of
// values, and both are correspondingly slow. A descent method fails for the
// opposite reason: at the minimum of |x| the subgradient it evaluates is +-1,
// so the gradient never becomes small, the line search finds no acceptable
// step, and the run stops on a failure while standing exactly on the answer.
//
// The bundle method uses the derivatives and is not fooled by them. It keeps a
// COLLECTION of subgradients from the points it has visited, and builds from
// them a piecewise-linear model of f -- the maximum of the linearisations -- which
// represents a kink as a kink instead of trying to smooth one away. The step is
// then the minimiser of that model plus a proximal term, and the term is what
// keeps the step inside the region where the model is believed.
//
// Written for a CONVEX f, which is the case with a theorem. It runs on a
// non-convex one and the safeguards below say what is clipped to make that
// possible, but the guarantee does not survive.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "objective.h"
#include "bounded.h"
#include "loop_support.h"

using namespace optimizers7;

namespace {

// Euclidean projection onto the unit simplex, by the sort-based algorithm. The
// bundle is small, so an exact projection costs nothing and spares us wondering
// whether an approximate one is why a run stalled.
arma::vec project_simplex(arma::vec y) {
  const arma::uword k = y.n_elem;
  if (k == 1) { y[0] = 1.0; return y; }
  arma::vec u = arma::sort(y, "descend");
  double css = 0.0, theta = 0.0;
  arma::uword rho = 0;
  for (arma::uword j = 0; j < k; ++j) {
    css += u[j];
    const double t = (css - 1.0) / static_cast<double>(j + 1);
    if (u[j] - t > 0.0) { rho = j; theta = t; }
  }
  (void) rho;
  arma::vec out = y - theta;
  out.transform([](double v) { return v > 0.0 ? v : 0.0; });
  const double s = arma::accu(out);
  if (s > 0.0) out /= s; else out.fill(1.0 / static_cast<double>(k));
  return out;
}

// The dual of the bundle subproblem:
//
//   min over lambda in the simplex of  (t/2) || G lambda ||^2 + alpha' lambda
//
// which is what the primal
//
//   min over d of  max_j ( -alpha_j + g_j' d ) + (1/2t) || d ||^2
//
// becomes once the max is written as a convex combination. Solved by FISTA,
// because the feasible set is a simplex whose projection is exact and cheap, so
// an accelerated projected gradient reaches a far tighter solution than the
// stopping tolerance needs in a few hundred iterations on a problem of at most
// twenty variables.
arma::vec solve_dual(const arma::mat& G, const arma::vec& alpha, double t,
                     int iters, double tol) {
  const arma::uword k = G.n_cols;
  if (k == 1) return arma::vec(1, arma::fill::ones);

  const arma::mat Q = t * (G.t() * G);
  arma::vec ev;
  double L = 1.0;
  if (arma::eig_sym(ev, Q) && ev.n_elem) L = std::max(ev.max(), 1e-12);
  else L = std::max(arma::norm(Q, "fro"), 1e-12);

  arma::vec lam(k, arma::fill::ones);
  lam /= static_cast<double>(k);
  arma::vec y = lam, lam_old = lam;
  double theta = 1.0;

  for (int i = 0; i < iters; ++i) {
    const arma::vec grad = Q * y + alpha;
    lam_old = lam;
    lam = project_simplex(y - grad / L);
    const double theta_new = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * theta * theta));
    y = lam + ((theta - 1.0) / theta_new) * (lam - lam_old);
    theta = theta_new;
    if (arma::norm(lam - lam_old, "inf") < tol) break;
  }
  return lam;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List bundle_run(Rcpp::List spec,
                      arma::vec par,
                      SEXP criterion,
                      SEXP crit_fn,
                      double t0,
                      double t_min,
                      double t_max,
                      double m_serious,
                      int bundle_size,
                      int qp_iters,
                      double qp_tol,
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
  arma::vec xhat = par;
  double fhat = obj->value(xhat);
  if (!std::isfinite(fhat)) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  // Each bundle element is stored as the AFFINE FUNCTION it defines,
  // l_j(x) = a_j + g_j' x, rather than as the point it came from. The
  // linearisation error at the current centre is then f(xhat) - l_j(xhat), one
  // subtraction, and nothing has to be recomputed when the centre moves. It
  // also lets the aggregate be stored exactly like an ordinary element, which
  // is what makes bundle compression a one-line operation below.
  std::vector<arma::vec> Gs;
  std::vector<double>    as;

  // t0 IS A STEP LENGTH, NOT A BARE MULTIPLIER.
  //
  // The step is d = -t p, so a raw t makes the first step as long as the
  // gradient happens to be big -- 233 units on Rosenbrock from its usual start.
  // That lands where the objective is 1e11 and its gradient 1e15; the next null
  // step halves t while the gradient has squared, so t can never catch up, and
  // the run spends its whole budget taking rejected steps and returns the point
  // it started from. Measured before the fix: 0 serious steps in 3000
  // iterations on rosenbrock, beale and powell alike, with the aggregate
  // subgradient at 1e196.
  //
  // Dividing by the gradient's size makes the first step of length t0 in
  // PARAMETER space, which is the same normalisation a line search performs
  // when it starts at 1 along a unit direction. The bounds move with it, since
  // they bound the same quantity. After it: all eight battery problems solved.
  double tscale = 1.0;
  {
    const arma::vec g0 = obj->grad(xhat);
    Gs.push_back(g0);
    as.push_back(fhat - arma::dot(g0, xhat));
    const double gn = g0.is_empty() ? 0.0 : arma::abs(g0).max();
    if (std::isfinite(gn) && gn > 1.0) tscale = 1.0 / gn;
  }
  t0 *= tscale;
  t_min *= tscale;
  t_max *= tscale;

  double t = t0;
  Trace tr;
  bool converged = false;
  std::string stopped_by = "maxit";
  std::string note = "";
  int it = 0;
  int serious = 0, null_steps = 0;
  arma::vec p_agg(p, arma::fill::zeros);
  double stat = R_PosInf;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective  stationarity        step  safeguard\n";
  }

  double f_prev = fhat;
  arma::vec x_prev = xhat;
  bool have_prev = false;

  for (it = 1; it <= maxit; ++it) {
    Rcpp::checkUserInterrupt();

    const arma::uword k = Gs.size();
    arma::mat G(p, k);
    arma::vec alpha(k);
    for (arma::uword j = 0; j < k; ++j) {
      G.col(j) = Gs[j];
      // Non-negative by convexity. On a non-convex objective it can come out
      // negative, and clipping is what lets the method run there at all -- it
      // is the standard convexification, and it is also exactly where the
      // convergence theory stops applying.
      alpha[j] = std::max(fhat - as[j] - arma::dot(Gs[j], xhat), 0.0);
    }

    const arma::vec lam = solve_dual(G, alpha, t, qp_iters, qp_tol);
    p_agg = G * lam;
    const double a_agg = arma::dot(alpha, lam);

    // The predicted decrease, used for the descent test: how much better the
    // model believes the step it is about to propose will be.
    const double v = -(t * arma::dot(p_agg, p_agg) + a_agg);

    // The OPTIMALITY ESTIMATE, which is what the stopping rule watches, and it
    // is deliberately not -v. Both are non-negative and both vanish when 0 lies
    // in the convex hull of the collected subgradients with no linearisation
    // error -- the computable stand-in for 0 being in the subdifferential -- but
    // -v carries a factor of t, and t is halved at every null step. So -v can
    // be driven below any tolerance by the TRUST PARAMETER SHRINKING rather
    // than by the point becoming stationary, and the run then reports success
    // while standing somewhere the model still says is steep.
    //
    // Not hypothetical: on an objective finite nowhere but its starting point,
    // every trial is refused, t halves down past 1e-8, and -v follows it under
    // the tolerance. The run reported converged = TRUE at the point it began,
    // having taken zero serious and zero null steps. Dropping t leaves the
    // model's own claim, which no amount of shrinking can flatter.
    stat = arma::dot(p_agg, p_agg) + a_agg;

    // The second half of the same guard. Scaling t0 stops the runaway from
    // starting, but nothing stops an objective that grows fast enough from
    // producing a subgradient whose square overflows -- and once G'G holds an
    // infinity the subproblem returns nonsense and every iteration after it is
    // wasted. Ending here with a diagnosis is worth far more than spending the
    // whole budget and reporting the starting point as though it were an answer.
    if (!p_agg.is_finite() || !std::isfinite(v)) {
      note = "the model diverged: a subgradient overflowed. Try a smaller t0.";
      stopped_by = "failed";
      if (keep_trace) tr.push_stat(it, fhat, R_PosInf, 0.0, "diverged");
      break;
    }

    if (have_prev &&
        ask_criterion(crit_fn, criterion, it - 1, fhat, f_prev, xhat, x_prev,
                      arma::vec(), false, true, stat, true)) {
      converged = true;
      stopped_by = "criterion";
      it = it - 1;
      break;
    }
    f_prev = fhat;
    x_prev = xhat;
    have_prev = true;

    const arma::vec d = -t * p_agg;
    const arma::vec y = xhat + d;
    const double fy = obj->value(y);
    std::string guard;

    if (std::isfinite(fy)) {
      const arma::vec gy = obj->grad(y);
      Gs.push_back(gy);
      as.push_back(fy - arma::dot(gy, y));

      // A serious step needs a decrease that is a fixed FRACTION of the one the
      // model promised -- the same bargain as an Armijo condition, and for the
      // same reason: accepting any decrease at all lets a sequence of ever
      // tinier improvements masquerade as progress.
      if (fy <= fhat + m_serious * v) {
        xhat = y;
        fhat = fy;
        ++serious;
        t = std::min(t * 2.0, t_max);
        guard = "serious";
      } else {
        ++null_steps;
        // The model was too optimistic over that distance, so trust it less
        // far. Kiwiel's rule chooses the factor by a curvature estimate; this
        // is the safeguarded halving, which is cruder and cannot run away.
        t = std::max(t * 0.5, t_min);
        guard = "null";
      }
    } else {
      // Stepped outside the domain. Nothing is learnt, so nothing is added --
      // only the trust region shrinks.
      const double t_new = std::max(t * 0.5, t_min);
      if (t_new == t) {
        // Already at the floor, so halving can do nothing more and every
        // further iteration would propose the same rejected step. The trigger
        // needs no chosen constant: t_min IS the limit of the only control the
        // method has left, and reaching it while still being refused is the
        // definition of stuck. Without this the run spends its whole budget
        // re-proposing one step and reports the point it started from.
        note = "every trial step left the domain, even at the smallest "
               "trust region";
        stopped_by = "failed";
        if (keep_trace) tr.push_stat(it, fhat, stat, 0.0, "trust floor");
        break;
      }
      t = t_new;
      guard = "step rejected";
    }

    // Compression. When the bundle is full the OLDEST elements are replaced by
    // the aggregate, which is the affine function the dual solution just
    // produced. Discarding them outright would lose the information they carry
    // and can stall the method; the aggregate keeps a summary of all of it in
    // one element, which is what makes a bounded bundle safe.
    if (static_cast<int>(Gs.size()) > bundle_size) {
      const double a_of_agg = fhat - a_agg - arma::dot(p_agg, xhat);
      std::vector<arma::vec> Gn;
      std::vector<double> an;
      Gn.push_back(p_agg);
      an.push_back(a_of_agg);
      const int keep = bundle_size - 1;
      for (int j = static_cast<int>(Gs.size()) - keep;
           j < static_cast<int>(Gs.size()); ++j) {
        Gn.push_back(Gs[j]);
        an.push_back(as[j]);
      }
      Gs.swap(Gn);
      as.swap(an);
      if (guard == "null") guard = "null, compressed";
    }

    if (keep_trace) tr.push_stat(it, fhat, stat, arma::norm(d, 2), guard);

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << fhat
                  << std::setw(14) << stat
                  << std::setw(12) << arma::norm(d, 2)
                  << "  " << guard << "\n";
    }

    if (inner->n_value >= max_eval) {
      note = "evaluation budget exhausted";
      stopped_by = "max_eval";
      break;
    }
  }
  if (it > maxit) it = maxit;

  {
    const std::string s = std::to_string(serious) + " serious, " +
      std::to_string(null_steps) + " null";
    note = note.empty() ? s : note + "; " + s;
  }

  if (verbose) {
    Rcpp::Rcout << "  stopped after " << it << " iterations ("
                << stopped_by << "), objective " << fhat << "\n";
  }

  // The AGGREGATE subgradient, not a subgradient at the point. At a solution
  // this goes to zero while any individual subgradient does not, so it is the
  // one quantity of that shape worth reporting.
  arma::vec g_out = bounded ? wrap->report_grad(xhat, p_agg) : p_agg;

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(bounded ? wrap->report_par(xhat)
                                                    : xhat),
    Rcpp::Named("value")      = fhat,
    Rcpp::Named("gradient")   = as_r_vector(g_out),
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
