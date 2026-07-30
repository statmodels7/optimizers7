// The first algorithm, and the one whose job is to prove the frame.
//
// Gradient descent with backtracking is not interesting as an optimiser. It is
// here because everything around it is the part that has to be right: the
// objective abstraction, the callback into R for the stopping rule, the trace,
// the reporting, the interrupt handling and the evaluation budget. Every
// algorithm that follows slots into exactly this shape.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "objective.h"
#include "line_search.h"

using namespace optimizers7;

// Ask the R-side criterion whether to stop. One callback per iteration, not per
// evaluation, so the cost is negligible beside the objective -- and in exchange
// a user-written criterion works exactly like a built-in one.
static bool criterion_met(SEXP crit_fn, SEXP criterion, int iter,
                          double f_new, double f_old,
                          const arma::vec& x_new, const arma::vec& x_old,
                          const arma::vec& g, bool have_g, bool have_old) {
  // Built with the optional slots empty and filled in afterwards: a ternary
  // whose branches are a NumericVector and R_NilValue has no common type, and
  // forcing one loses the protection the Rcpp wrapper provides.
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
  if (have_g) {
    state["gradient"] = as_r_vector(g);
  }
  Rcpp::LogicalVector out = Rcpp::Function(crit_fn)(criterion, state);
  return out.size() == 1 && out[0] == TRUE;
}

// [[Rcpp::export]]
Rcpp::List gd_run(Rcpp::List spec,
                  arma::vec par,
                  SEXP criterion,
                  SEXP crit_fn,
                  int maxit,
                  int max_eval,
                  bool verbose,
                  int refresh,
                  bool keep_trace,
                  double step,
                  Rcpp::List line_search) {

  std::unique_ptr<Objective> obj(make_objective(spec));
  const LineSearchSpec ls = LineSearchSpec::from_list(line_search);

  arma::vec x = par;
  double f = obj->value(x);
  if (!std::isfinite(f)) {
    Rcpp::stop("The objective is not finite at the starting value.");
  }

  std::vector<int>    tr_iter;
  std::vector<double> tr_value;
  std::vector<double> tr_gnorm;
  std::vector<double> tr_step;
  std::vector<std::string> tr_guard;

  bool converged = false;
  std::string stopped_by = "maxit";
  std::string note = "";
  int it = 0;

  if (verbose) {
    Rcpp::Rcout << "  iter        objective     |grad|max        step  safeguard\n";
  }

  for (it = 1; it <= maxit; ++it) {
    // A C++ loop that ignores interrupts is a session that cannot be stopped.
    Rcpp::checkUserInterrupt();

    arma::vec g = obj->grad(x);
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

    // The direction is the steepest-descent one; how far to go along it is
    // the line search's business, written once in line_search.h and shared
    // with every method that will follow.
    const arma::vec d = -g;
    LineSearchResult step_out = run_line_search(ls, *obj, x, f, g, d, step);

    if (!step_out.ok) {
      note = "the line search found no acceptable step";
      stopped_by = "failed";
      guard = step_out.guard;
      if (keep_trace) {
        tr_iter.push_back(it); tr_value.push_back(f);
        tr_gnorm.push_back(gnorm); tr_step.push_back(0.0);
        tr_guard.push_back(guard);
      }
      break;
    }
    if (step_out.guard != "none") guard = step_out.guard;

    const double s = step_out.step;
    const arma::vec x_new = step_out.x_new;
    const double f_new = step_out.f_new;

    const arma::vec x_old = x;
    const double f_old = f;
    x = x_new;
    f = f_new;

    if (keep_trace) {
      tr_iter.push_back(it); tr_value.push_back(f);
      tr_gnorm.push_back(gnorm); tr_step.push_back(s);
      tr_guard.push_back(guard);
    }

    if (verbose && refresh > 0 && (it % refresh == 0)) {
      Rcpp::Rcout << std::setw(6) << it
                  << std::setw(17) << f
                  << std::setw(14) << gnorm
                  << std::setw(12) << s
                  << "  " << guard << "\n";
    }

    // The Wolfe search has already evaluated the gradient at the accepted
    // point; recomputing it would double the gradient count for nothing.
    arma::vec g_new = step_out.g_new_valid ? step_out.g_new : obj->grad(x);
    if (criterion_met(crit_fn, criterion, it, f, f_old, x, x_old,
                      g_new, true, true)) {
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

  arma::vec g_final = obj->grad(x);

  if (verbose) {
    Rcpp::Rcout << "  stopped after " << it << " iterations ("
                << stopped_by << "), objective " << f << "\n";
  }

  // SEXP, not Rcpp::List: assigning R_NilValue to an Rcpp::List builds an
  // EMPTY LIST rather than NULL, so the R side would see list() where it tests
  // for NULL, and trace$value would be NULL inside code that expected a column.
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
    Rcpp::Named("gradient")   = as_r_vector(g_final),
    Rcpp::Named("iterations") = it,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("stopped_by") = stopped_by,
    Rcpp::Named("message")    = note,
    Rcpp::Named("n_value")    = obj->n_value,
    Rcpp::Named("n_grad")     = obj->n_grad,
    Rcpp::Named("trace")      = trace
  );
}
