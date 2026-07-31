// Adam, and the reason it needs a loop of its own.
//
// Every other method in the package is a direction plus a line search, so they
// share descent.cpp and differ only in the Direction they carry. Adam is not
// one of those, on two counts, and forcing it into that shape would have
// changed the algorithm rather than reused the code:
//
//   * it takes no line search. Its step length is set by the accumulated
//     second-moment estimate, and imposing a sufficient-decrease test on top
//     would be a different algorithm wearing Adam's name;
//   * it is not a descent method. The objective is allowed to rise, and that
//     freedom is most of why it copes with a gradient that points downhill
//     only on average.
//
// WHAT IT DELIBERATELY DOES NOT DO IS DRAW ITS OWN SUBSAMPLES. An optimiser
// does not know what an observation is, and a version that did would have to be
// told -- which means a second kind of objective, a second kind of criterion,
// and a rule for which combinations are allowed. All of that belonged to
// whichever caller already knows the answer. Adam here optimises the objective
// it is given; if that objective is stochastic, Adam is stochastic, and neither
// side needs to negotiate about it.
//
// What is shared is shared: the stopping rule, the trace and the reporting all
// come from loop_support.h, so the two loops cannot drift in what they report.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "objective.h"
#include "bounded.h"
#include "loop_support.h"

using namespace optimizers7;

// [[Rcpp::export]]
Rcpp::List adam_run(Rcpp::List spec,
                    arma::vec par,
                    SEXP criterion,
                    SEXP crit_fn,
                    double alpha,
                    double beta1,
                    double beta2,
                    double eps,
                    double decay,
                    bool amsgrad,
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

  double f = obj->value(x);
  if (!std::isfinite(f)) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  arma::vec mv(p, arma::fill::zeros);   // first moment
  arma::vec vv(p, arma::fill::zeros);   // second moment
  arma::vec vmax(p, arma::fill::zeros); // the AMSGrad running maximum

  Trace tr;
  bool converged = false;
  std::string stopped_by = "maxit";
  std::string note = "";
  int it = 0;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective     |grad|max        step  safeguard\n";
  }

  // Bookkeeping for the criterion, which is asked at the TOP of the next
  // iteration rather than the bottom of this one. The reason is that Adam
  // computes its gradient before stepping, so a rule asked at the bottom would
  // be shown the gradient at the point the run had just LEFT. Deferring by one
  // pass costs nothing -- the gradient is needed anyway -- and lets crit_grad()
  // mean what it says.
  bool have_prev = false;
  double f_prev = 0.0;
  arma::vec x_prev(p, arma::fill::zeros);

  for (it = 1; it <= maxit; ++it) {
    Rcpp::checkUserInterrupt();

    const arma::vec g = obj->grad(x);
    std::string guard = "none";
    const double gnorm = g.is_empty() ? 0.0 : arma::abs(g).max();

    if (have_prev &&
        ask_criterion(crit_fn, criterion, it - 1, f, f_prev, x, x_prev,
                      g, true, true)) {
      converged = true;
      stopped_by = "criterion";
      it = it - 1;
      break;
    }

    if (!g.is_finite()) {
      note = "gradient not finite; stopped";
      stopped_by = "failed";
      guard = "gradient non-finite";
      if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
      break;
    }

    // The learning rate schedule. Constant by default, which is Adam as
    // published; with decay > 0 it is O(1/t), the Robbins-Monro condition a
    // run on a noisy objective needs to settle rather than rattle around the
    // optimum.
    const double at = alpha / (1.0 + decay * (it - 1));

    mv = beta1 * mv + (1.0 - beta1) * g;
    vv = beta2 * vv + (1.0 - beta2) * (g % g);

    // Bias correction. Both moments start at zero, so early on they are pulled
    // towards it; dividing by 1 - beta^t undoes exactly that and matters most
    // in the first few iterations, where an uncorrected Adam barely moves.
    const double c1 = 1.0 - std::pow(beta1, static_cast<double>(it));
    const double c2 = 1.0 - std::pow(beta2, static_cast<double>(it));

    arma::vec second = vv;
    if (amsgrad) {
      vmax = arma::max(vmax, vv);
      second = vmax;
      if (arma::any(vmax > vv)) guard = "amsgrad floor";
    }

    const arma::vec mhat = mv / c1;
    const arma::vec vhat = second / c2;
    const arma::vec dx = at * (mhat / (arma::sqrt(vhat) + eps));

    if (!dx.is_finite()) {
      note = "the update is not finite; stopped";
      stopped_by = "failed";
      guard = "update non-finite";
      if (keep_trace) tr.push(it, f, gnorm, 0.0, guard);
      break;
    }

    x_prev = x;
    f_prev = f;
    have_prev = true;
    x -= dx;

    // Evaluated every iteration, though the algorithm itself never reads it.
    // Making it conditional saved one call and cost a silent defect: with the
    // value left stale, a rule comparing f_new against f_old would compare a
    // number against itself and fire at once. One evaluation is what every
    // other method here already pays.
    f = obj->value(x);
    const double step = arma::abs(dx).max();

    if (keep_trace) tr.push(it, f, gnorm, step, guard);

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << f
                  << std::setw(14) << gnorm
                  << std::setw(12) << step
                  << "  " << guard << "\n";
    }

    if (inner->n_value >= max_eval) {
      note = "evaluation budget exhausted";
      stopped_by = "max_eval";
      break;
    }
  }
  if (it > maxit) it = maxit;

  // The value and gradient at the point actually reached, rather than whatever
  // the last iteration happened to have computed.
  const double f_end = obj->value(x);
  const arma::vec g_end = obj->grad(x);

  if (verbose) {
    Rcpp::Rcout << "  stopped after " << it << " iterations ("
                << stopped_by << "), objective " << f_end << "\n";
  }

  arma::vec par_out = bounded ? wrap->report_par(x) : x;
  arma::vec g_out   = bounded ? wrap->report_grad(x, g_end) : g_end;

  return Rcpp::List::create(
    Rcpp::Named("par")        = as_r_vector(par_out),
    Rcpp::Named("value")      = f_end,
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
