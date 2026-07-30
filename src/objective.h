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
  // throwing when the objective is undefined there: an optimiser must be able
  // to propose a bad point and be told so, which is a normal event in a line
  // search and not an error.
  virtual double value(const arma::vec& x) = 0;

  // The gradient at x. Only called when has_gradient() is true; otherwise the
  // caller uses fd_gradient() below.
  virtual arma::vec gradient(const arma::vec& x) = 0;

  virtual bool has_gradient() const = 0;

  int n_value = 0;
  int n_grad  = 0;

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
};


// --- an ordinary R function -------------------------------------------------
//
// Every evaluation is a callback into R. That is the cost of accepting an
// arbitrary R closure, and it is why a compiled objective is worth offering.

class RObjective : public Objective {
public:
  RObjective(SEXP fn, SEXP gr) : fn_(fn), gr_(gr), has_gr_(gr != R_NilValue) {}

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

  bool has_gradient() const { return has_gr_; }

private:
  SEXP fn_;
  SEXP gr_;
  bool has_gr_;
};


// --- a compiled objective ---------------------------------------------------
//
// Function pointers, called directly. The loop never returns to R, which is the
// whole reason the algorithms are in C++ at all.

typedef double   (*fn_ptr_t)(const arma::vec&);
typedef arma::vec (*gr_ptr_t)(const arma::vec&);

class CppObjective : public Objective {
public:
  CppObjective(SEXP fn, SEXP gr) {
    Rcpp::XPtr<fn_ptr_t> f(fn);
    fn_ = *f;
    if (gr != R_NilValue) {
      Rcpp::XPtr<gr_ptr_t> g(gr);
      gr_ = *g;
      has_gr_ = true;
    }
  }

  double value(const arma::vec& x) {
    ++n_value;
    return fn_(x);
  }

  arma::vec gradient(const arma::vec& x) {
    ++n_grad;
    return gr_(x);
  }

  bool has_gradient() const { return has_gr_; }

private:
  fn_ptr_t fn_ = nullptr;
  gr_ptr_t gr_ = nullptr;
  bool has_gr_ = false;
};


// --- a finite sum over observations -----------------------------------------
//
// Evaluated on the whole sample here. The point of the class is that it KNOWS
// it is a sum, so a stochastic method can ask it for a subsample; the frame
// carries that knowledge from the start even though the only algorithm in the
// skeleton never uses it.

class FiniteSumObjective : public Objective {
public:
  FiniteSumObjective(SEXP fn, SEXP gr, int n)
    : fn_(fn), gr_(gr), n_(n), has_gr_(gr != R_NilValue) {}

  double value(const arma::vec& x) {
    ++n_value;
    Rcpp::NumericVector out = Rcpp::Function(fn_)(as_r_vector(x), all_idx());
    return out[0];
  }

  arma::vec gradient(const arma::vec& x) {
    ++n_grad;
    Rcpp::NumericVector out = Rcpp::Function(gr_)(as_r_vector(x), all_idx());
    return Rcpp::as<arma::vec>(out);
  }

  // A subsample of the terms, for a stochastic method. Draws through R's RNG so
  // that set.seed() governs the run and it can be reproduced.
  Rcpp::IntegerVector sample_idx(int m) const {
    Rcpp::Function sample_int("sample.int", R_BaseNamespace);
    return sample_int(n_, m);
  }

  double value_on(const arma::vec& x, Rcpp::IntegerVector idx) {
    ++n_value;
    Rcpp::NumericVector out = Rcpp::Function(fn_)(as_r_vector(x), idx);
    return out[0];
  }

  arma::vec gradient_on(const arma::vec& x, Rcpp::IntegerVector idx) {
    ++n_grad;
    Rcpp::NumericVector out = Rcpp::Function(gr_)(as_r_vector(x), idx);
    return Rcpp::as<arma::vec>(out);
  }

  bool has_gradient() const { return has_gr_; }
  int n_terms() const { return n_; }

private:
  Rcpp::IntegerVector all_idx() const {
    Rcpp::IntegerVector idx(n_);
    for (int i = 0; i < n_; ++i) idx[i] = i + 1;
    return idx;
  }

  SEXP fn_;
  SEXP gr_;
  int n_;
  bool has_gr_;
};


// Build the right Objective from the handle as_objective() produced. The single
// place in the package that knows the three shapes apart.
inline Objective* make_objective(Rcpp::List spec) {
  std::string kind = Rcpp::as<std::string>(spec["kind"]);
  if (kind == "r") {
    return new RObjective(spec["fn"], spec["gr"]);
  } else if (kind == "cpp") {
    return new CppObjective(spec["fn_ptr"], spec["gr_ptr"]);
  } else if (kind == "finite_sum") {
    return new FiniteSumObjective(spec["fn_idx"], spec["gr_idx"],
                                  Rcpp::as<int>(spec["n"]));
  }
  Rcpp::stop("Unknown objective kind '%s'.", kind);
  return nullptr;
}

} // namespace optimizers7

#endif
