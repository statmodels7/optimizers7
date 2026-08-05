#ifndef OPTIMIZERS7_OBJECTIVE_H
#define OPTIMIZERS7_OBJECTIVE_H

// The one interface every algorithm in this package is written against.
//
// An algorithm asks an Objective for values and gradients and knows nothing
// else about it: not whether it is an R closure, a sum over observations, or a
// compiled function reached through a pointer. That is what lets each algorithm
// be written once while the objective arrives in three shapes.
//
// Evaluation counts live here rather than in the algorithms, so that every
// algorithm reports them the same way and none can forget to.

#include <RcppArmadillo.h>

namespace optimizers7 {

// Rcpp::wrap() on an arma::vec produces a one-COLUMN MATRIX, not a plain
// numeric vector. Handing that to a user's R function is a real bug and a
// quiet one: `y[idx] - par` then does vector-array arithmetic, R warns about
// recycling, and code that indexes or names the parameter behaves oddly for
// reasons that point nowhere near here. Everything crossing back into R goes
// through this instead.
inline Rcpp::NumericVector as_r_vector(const arma::vec& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

class Objective {
public:
  virtual ~Objective() {}

  // The value at x. Implementations return a non-finite number rather than
  // throwing when the objective is undefined there: an optimizer must be able
  // to propose a bad point and be told so, which is a normal event in a line
  // search and not an error.
  virtual double value(const arma::vec& x) = 0;

  // The gradient at x. Only called when has_gradient() is true; otherwise the
  // caller uses fd_gradient() below.
  virtual arma::vec gradient(const arma::vec& x) = 0;

  // The Hessian at x. Only called when has_hessian() is true.
  virtual arma::mat hessian(const arma::vec& x) {
    Rcpp::stop("This objective supplies no Hessian.");
  }

  virtual bool has_gradient() const = 0;
  virtual bool has_hessian() const { return false; }

  int n_value = 0;
  int n_grad  = 0;
  int n_hess  = 0;

  // Central-difference gradient, for an objective that supplies none. One
  // differentiation, never nested: the same rule the rest of the toolkit
  // follows, and the reason it matters is that composing finite differences
  // multiplies their error.
  arma::vec fd_gradient(const arma::vec& x) {
    const arma::uword p = x.n_elem;
    arma::vec g(p, arma::fill::zeros);
    const double eps = std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0);
    arma::vec xp = x, xm = x;
    for (arma::uword j = 0; j < p; ++j) {
      const double h = eps * std::max(1.0, std::abs(x[j]));
      xp[j] = x[j] + h;
      xm[j] = x[j] - h;
      const double fp = value(xp);
      const double fm = value(xm);
      xp[j] = x[j];
      xm[j] = x[j];
      g[j] = (fp - fm) / (2.0 * h);
    }
    return g;
  }

  // What an algorithm should call: the analytic gradient when there is one and
  // the difference otherwise, so no algorithm has to branch.
  arma::vec grad(const arma::vec& x) {
    return has_gradient() ? gradient(x) : fd_gradient(x);
  }

  // Central-difference Hessian, obtained by differentiating grad() ONCE. When
  // the objective supplies an analytic gradient that is a single numerical
  // differentiation and the result is good; when it does not, grad() is itself
  // a difference and this is the one place in the package where two are
  // composed. It is unavoidable -- a Hessian from values alone has to be -- but
  // it is why newton() warns that a density-only objective is better served by
  // bfgs(), which never needs one.
  arma::mat fd_hessian(const arma::vec& x) {
    const arma::uword p = x.n_elem;
    arma::mat H(p, p, arma::fill::zeros);
    const double eps = std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0);
    arma::vec xp = x, xm = x;
    for (arma::uword j = 0; j < p; ++j) {
      const double h = eps * std::max(1.0, std::abs(x[j]));
      xp[j] = x[j] + h;
      xm[j] = x[j] - h;
      const arma::vec gp = grad(xp);
      const arma::vec gm = grad(xm);
      xp[j] = x[j];
      xm[j] = x[j];
      H.col(j) = (gp - gm) / (2.0 * h);
    }
    // Symmetrize: the two triangles differ by the differencing error, and a
    // Hessian that is not exactly symmetric breaks eig_sym further on.
    return 0.5 * (H + H.t());
  }

  arma::mat hess(const arma::vec& x) {
    return has_hessian() ? hessian(x) : fd_hessian(x);
  }
};


// --- an ordinary R function -------------------------------------------------
//
// Every evaluation is a callback into R. That is the cost of accepting an
// arbitrary R closure, and it is why a compiled objective is worth offering.

class RObjective : public Objective {
public:
  RObjective(SEXP fn, SEXP gr, SEXP he)
    : fn_(fn), gr_(gr), he_(he),
      has_gr_(gr != R_NilValue), has_he_(he != R_NilValue) {}

  double value(const arma::vec& x) {
    ++n_value;
    Rcpp::NumericVector out = Rcpp::Function(fn_)(as_r_vector(x));
    if (out.size() != 1) {
      Rcpp::stop("The objective must return a single number.");
    }
    return out[0];
  }

  arma::vec gradient(const arma::vec& x) {
    ++n_grad;
    Rcpp::NumericVector out = Rcpp::Function(gr_)(as_r_vector(x));
    return Rcpp::as<arma::vec>(out);
  }

  arma::mat hessian(const arma::vec& x) {
    ++n_hess;
    Rcpp::NumericMatrix out = Rcpp::Function(he_)(as_r_vector(x));
    return Rcpp::as<arma::mat>(out);
  }

  bool has_gradient() const { return has_gr_; }
  bool has_hessian() const { return has_he_; }

private:
  SEXP fn_;
  SEXP gr_;
  SEXP he_;
  bool has_gr_;
  bool has_he_;
};


// Build the Objective from the handle as_objective() produced.
//
// There is one shape now, and there were three. A finite sum evaluable on a
// subset of its terms went when Adam stopped drawing minibatches; an objective
// behind a pair of bare C++ function pointers went because it could not carry
// data -- double(*)(const arma::vec&) has no closure, so any real statistical
// objective had to keep y and X in globals, which is not reentrant and means one
// model at a time. Measured, it was also SLOWER than vectorized R below about
// twenty thousand observations and never better than about twice as fast, since
// R's matrix arithmetic and Armadillo's call the same BLAS.
//
// The dispatch stays a generic on the R side even with one method registered,
// because that is the extension point: a caller with its own kind of objective
// registers a method and every algorithm accepts it. What went was a second
// shipped shape that nothing could use.
inline Objective* make_objective(Rcpp::List spec) {
  std::string kind = Rcpp::as<std::string>(spec["kind"]);
  if (kind == "r") {
    return new RObjective(spec["fn"], spec["gr"], spec["he"]);
  }
  Rcpp::stop("Unknown objective kind '%s'.", kind);
  return nullptr;
}

} // namespace optimizers7

#endif
