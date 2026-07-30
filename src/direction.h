#ifndef OPTIMIZERS7_DIRECTION_H
#define OPTIMIZERS7_DIRECTION_H

// What separates one descent method from another.
//
// Steepest descent, Newton, BFGS and L-BFGS share their entire iteration: take
// a direction, line search along it, accept, test the stopping rule. The only
// thing that differs is how the direction is produced and what the method
// remembers between iterations. That is this interface, and it is why the loop
// in descent.cpp is written once.

#include <RcppArmadillo.h>
#include "objective.h"

namespace optimizers7 {

class Direction {
public:
  virtual ~Direction() {}

  // The search direction at x, given the gradient there. `guard` is set when a
  // safeguard fires, and is carried into the trace and the verbose line.
  virtual arma::vec compute(Objective& obj, const arma::vec& x,
                            const arma::vec& g, double f,
                            std::string& guard) = 0;

  // Told about the step just taken, for the methods that remember. s = x_new -
  // x_old and y = g_new - g_old.
  virtual void update(const arma::vec& s, const arma::vec& y,
                      std::string& guard) {}

  virtual bool needs_hessian() const { return false; }
};


// --- steepest descent -------------------------------------------------------

class SteepestDescent : public Direction {
public:
  arma::vec compute(Objective&, const arma::vec&, const arma::vec& g, double,
                    std::string&) {
    return -g;
  }
};


// --- Newton -----------------------------------------------------------------
//
// The direction solves H d = -g. Away from a minimum H need not be positive
// definite, and then the solution is not a descent direction at all: it points
// at a saddle or a maximum. Modifying H until it is positive definite is what
// makes Newton usable rather than merely fast near the answer.
//
// Two strategies, both starting with a Cholesky attempt because when H is
// already positive definite that is both the test and the solve:
//
//   "eigen"  decompose and raise every eigenvalue below a floor to it. Costs a
//            symmetric eigendecomposition, and produces the best-conditioned
//            modification: the direction is the Newton one in the subspace
//            where the curvature is trustworthy and gradient-like elsewhere.
//
//   "ridge"  add tau*I with tau doubling until the Cholesky succeeds. Cheaper,
//            and the Levenberg idea: it interpolates between the Newton step
//            (tau = 0) and a scaled steepest-descent step (tau large).

class NewtonDirection : public Direction {
public:
  NewtonDirection(std::string mod, double floor_value)
    : mod_(mod), floor_(floor_value) {}

  bool needs_hessian() const { return true; }

  arma::vec compute(Objective& obj, const arma::vec& x, const arma::vec& g,
                    double, std::string& guard) {
    arma::mat H = obj.hess(x);

    if (!H.is_finite()) {
      guard = "hessian non-finite";
      return -g;
    }
    H = 0.5 * (H + H.t());

    arma::mat R;
    if (arma::chol(R, H)) {
      arma::vec d;
      if (arma::solve(d, H, -g, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) {
        if (arma::dot(g, d) < 0.0) return d;
      }
      // Positive definite and yet not a descent direction means the solve lost
      // accuracy; falling back is cheaper than pretending.
      guard = "hessian solve failed";
      return -g;
    }

    guard = "hessian modified";
    arma::mat Hm = (mod_ == "ridge") ? ridge(H) : eigen_floor(H);

    arma::vec d;
    if (!arma::solve(d, Hm, -g, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
        !d.is_finite() || arma::dot(g, d) >= 0.0) {
      guard = "hessian solve failed";
      return -g;
    }
    return d;
  }

private:
  arma::mat eigen_floor(const arma::mat& H) const {
    arma::vec val;
    arma::mat vec;
    if (!arma::eig_sym(val, vec, H)) return ridge(H);
    const double lo = std::max(floor_, floor_ * arma::abs(val).max());
    val.transform([&](double v) { return v < lo ? lo : v; });
    return vec * arma::diagmat(val) * vec.t();
  }

  arma::mat ridge(const arma::mat& H) const {
    const arma::uword p = H.n_rows;
    double tau = std::max(floor_, -H.diag().min() + floor_);
    arma::mat R;
    for (int k = 0; k < 60; ++k) {
      arma::mat Hm = H + tau * arma::eye(p, p);
      if (arma::chol(R, Hm)) return Hm;
      tau *= 2.0;
    }
    return arma::eye(p, p);
  }

  std::string mod_;
  double floor_;
};


// --- BFGS -------------------------------------------------------------------
//
// The inverse Hessian is approximated directly, so the direction is a
// matrix-vector product rather than a solve.
//
// The curvature guard is the part that matters. The update is only valid when
// s'y > 0, and when s'y is merely small the update is numerically poisonous:
// rho = 1/(s'y) blows up and the approximation is destroyed by one bad step.
// Skipping the update leaves a stale but sound matrix, which is far better than
// a fresh but corrupted one. A run that skips repeatedly has lost its curvature
// information anyway, so after enough skips the matrix is reset to the identity
// and the method restarts as steepest descent.

class BfgsDirection : public Direction {
public:
  BfgsDirection(arma::uword p, double curv_tol, int max_skip)
    : Binv_(arma::eye(p, p)), curv_tol_(curv_tol), max_skip_(max_skip) {}

  arma::vec compute(Objective&, const arma::vec&, const arma::vec& g, double,
                    std::string& guard) {
    arma::vec d = -(Binv_ * g);
    if (!d.is_finite() || arma::dot(g, d) >= 0.0) {
      guard = "bfgs reset";
      Binv_.eye();
      skips_ = 0;
      d = -g;
    }
    return d;
  }

  void update(const arma::vec& s, const arma::vec& y, std::string& guard) {
    const double sy = arma::dot(s, y);
    const double scale = arma::norm(s) * arma::norm(y);

    if (!(sy > curv_tol_ * scale) || !std::isfinite(sy)) {
      ++skips_;
      guard = "bfgs update skipped";
      if (skips_ >= max_skip_) {
        Binv_.eye();
        skips_ = 0;
        guard = "bfgs reset";
      }
      return;
    }
    skips_ = 0;

    // On the first accepted pair, scale the identity by s'y / y'y. Without it
    // the first Newton-like step is taken with a unit Hessian, which on a badly
    // scaled problem is a step of entirely the wrong magnitude.
    if (first_) {
      Binv_ *= sy / arma::dot(y, y);
      first_ = false;
    }

    const double rho = 1.0 / sy;
    const arma::uword p = s.n_elem;
    const arma::mat I = arma::eye(p, p);
    const arma::mat V = I - rho * s * y.t();
    Binv_ = V * Binv_ * V.t() + rho * s * s.t();
  }

private:
  arma::mat Binv_;
  double curv_tol_;
  int max_skip_;
  int skips_ = 0;
  bool first_ = true;
};


// --- L-BFGS -----------------------------------------------------------------
//
// The same approximation, never formed. Only the last `memory` secant pairs are
// kept and the product Binv * g is computed by the two-loop recursion, so the
// cost is O(m p) in time and memory instead of O(p^2). For a few dozen
// parameters that is no gain at all; for a few thousand it is the difference
// between possible and not.

class LbfgsDirection : public Direction {
public:
  LbfgsDirection(arma::uword m, double curv_tol)
    : m_(m), curv_tol_(curv_tol) {}

  arma::vec compute(Objective&, const arma::vec&, const arma::vec& g, double,
                    std::string& guard) {
    if (S_.empty()) return -g;

    arma::vec q = g;
    const std::size_t k = S_.size();
    std::vector<double> alpha(k);

    for (std::size_t i = k; i-- > 0; ) {
      alpha[i] = rho_[i] * arma::dot(S_[i], q);
      q -= alpha[i] * Y_[i];
    }

    // The initial scaling: the most recent pair's curvature, which is the
    // cheapest sensible estimate of the Hessian's scale and the reason L-BFGS
    // works without any matrix at all.
    const double gamma = arma::dot(S_.back(), Y_.back()) /
                         arma::dot(Y_.back(), Y_.back());
    arma::vec r = gamma * q;

    for (std::size_t i = 0; i < k; ++i) {
      const double beta = rho_[i] * arma::dot(Y_[i], r);
      r += S_[i] * (alpha[i] - beta);
    }

    arma::vec d = -r;
    if (!d.is_finite() || arma::dot(g, d) >= 0.0) {
      guard = "lbfgs reset";
      S_.clear(); Y_.clear(); rho_.clear();
      d = -g;
    }
    return d;
  }

  void update(const arma::vec& s, const arma::vec& y, std::string& guard) {
    const double sy = arma::dot(s, y);
    const double scale = arma::norm(s) * arma::norm(y);
    if (!(sy > curv_tol_ * scale) || !std::isfinite(sy)) {
      guard = "lbfgs update skipped";
      return;
    }
    if (S_.size() == m_) {
      S_.erase(S_.begin());
      Y_.erase(Y_.begin());
      rho_.erase(rho_.begin());
    }
    S_.push_back(s);
    Y_.push_back(y);
    rho_.push_back(1.0 / sy);
  }

private:
  arma::uword m_;
  double curv_tol_;
  std::vector<arma::vec> S_, Y_;
  std::vector<double> rho_;
};


inline Direction* make_direction(Rcpp::List method, arma::uword p) {
  std::string type = Rcpp::as<std::string>(method["type"]);
  if (type == "gradient_descent") {
    return new SteepestDescent();
  } else if (type == "newton") {
    return new NewtonDirection(Rcpp::as<std::string>(method["hessian_mod"]),
                               Rcpp::as<double>(method["floor"]));
  } else if (type == "bfgs") {
    return new BfgsDirection(p, Rcpp::as<double>(method["curv_tol"]),
                             Rcpp::as<int>(method["max_skip"]));
  } else if (type == "lbfgs") {
    return new LbfgsDirection(Rcpp::as<arma::uword>(method["memory"]),
                              Rcpp::as<double>(method["curv_tol"]));
  }
  Rcpp::stop("Unknown method '%s'.", type);
  return nullptr;
}

} // namespace optimizers7

#endif
