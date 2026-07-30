#ifndef OPTIMIZERS7_BOUNDED_H
#define OPTIMIZERS7_BOUNDED_H

// Box bounds, by reparametrisation.
//
// The toolkit's answer to a constraint is never to enforce it but to remove it:
// map the box onto the whole real line, optimise there, and map back. Every
// point the optimiser proposes is then admissible by construction, so there is
// no rejection step, no projection, and no boundary for a line search to trip
// over.
//
// THE TRANSFORMS ARE linkfunctions7's, WRITTEN OUT HERE.
//
// They are bounded_link()'s three cases and the identity, and they are
// implemented in C++ rather than called through R for one reason: the transform
// is evaluated on every objective evaluation, so calling back into R for it
// would destroy exactly the advantage a compiled objective exists to give. The
// duplication is real and is guarded by machine rather than by care: the test
// suite asserts that each of these agrees with linkinv(), dlinkinv() and
// d2linkinv() from linkfunctions7 to machine precision, so the two cannot drift
// apart unnoticed. That guard earns its keep -- linkfunctions7's exponential
// floor changed on 2026-07-30, and a copy made before then would now be wrong.

#include <RcppArmadillo.h>
#include "objective.h"

namespace optimizers7 {

// linkfunctions7's exp_floor: (24 / double.xmax)^(1/4), the smallest value
// keeping the log link's fourth forward derivative -6/theta^4 finite.
inline double exp_floor() {
  static const double v =
    std::pow(24.0 / std::numeric_limits<double>::max(), 0.25);
  return v;
}

inline double exp_floored(double eta) {
  const double e = std::exp(eta);
  return e < exp_floor() ? exp_floor() : e;
}

enum BoundKind { BOUND_NONE = 0, BOUND_LOWER, BOUND_UPPER, BOUND_BOTH };

struct Coord {
  BoundKind kind = BOUND_NONE;
  double lwr = R_NegInf;
  double upr = R_PosInf;
};

// theta = h(eta)
inline double h_value(const Coord& c, double eta) {
  switch (c.kind) {
    case BOUND_LOWER: return c.lwr + exp_floored(eta);
    case BOUND_UPPER: return c.upr - exp_floored(eta);
    case BOUND_BOTH: {
      const double W = c.upr - c.lwr;
      // plogis, written so that a large |eta| cannot overflow
      const double p = eta >= 0.0 ? 1.0 / (1.0 + std::exp(-eta))
                                  : std::exp(eta) / (1.0 + std::exp(eta));
      return c.lwr + W * p;
    }
    default: return eta;
  }
}

// h'(eta)
inline double h_deriv1(const Coord& c, double eta) {
  switch (c.kind) {
    case BOUND_LOWER: return exp_floored(eta);
    case BOUND_UPPER: return -exp_floored(eta);
    case BOUND_BOTH: {
      const double W = c.upr - c.lwr;
      const double p = eta >= 0.0 ? 1.0 / (1.0 + std::exp(-eta))
                                  : std::exp(eta) / (1.0 + std::exp(eta));
      return W * p * (1.0 - p);
    }
    default: return 1.0;
  }
}

// h''(eta)
inline double h_deriv2(const Coord& c, double eta) {
  switch (c.kind) {
    case BOUND_LOWER: return exp_floored(eta);
    case BOUND_UPPER: return -exp_floored(eta);
    case BOUND_BOTH: {
      const double W = c.upr - c.lwr;
      const double p = eta >= 0.0 ? 1.0 / (1.0 + std::exp(-eta))
                                  : std::exp(eta) / (1.0 + std::exp(eta));
      return W * p * (1.0 - p) * (1.0 - 2.0 * p);
    }
    default: return 0.0;
  }
}

// eta = g(theta), for the starting value
inline double g_value(const Coord& c, double theta) {
  switch (c.kind) {
    case BOUND_LOWER: return std::log(theta - c.lwr);
    case BOUND_UPPER: return std::log(c.upr - theta);
    case BOUND_BOTH: {
      const double p = (theta - c.lwr) / (c.upr - c.lwr);
      return std::log(p / (1.0 - p));
    }
    default: return theta;
  }
}

inline arma::vec to_theta(const std::vector<Coord>& cs, const arma::vec& eta) {
  arma::vec th(eta.n_elem);
  for (arma::uword i = 0; i < eta.n_elem; ++i) th[i] = h_value(cs[i], eta[i]);
  return th;
}

inline arma::vec to_eta(const std::vector<Coord>& cs, const arma::vec& theta) {
  arma::vec e(theta.n_elem);
  for (arma::uword i = 0; i < theta.n_elem; ++i) e[i] = g_value(cs[i], theta[i]);
  return e;
}


// The decorator. An algorithm sees an ordinary unconstrained objective in eta
// and never learns that a box was involved.
//
// The chain rule, with a DIAGONAL Jacobian because each coordinate carries its
// own transform and nothing else:
//
//   df/deta_i          = (df/dtheta_i) h_i'
//   d2f/deta_i deta_j  = (d2f/dtheta_i dtheta_j) h_i' h_j'
//                        + [i = j] (df/dtheta_i) h_i''
//
// The second term appears only on the diagonal, and only because h is not
// affine; it is the whole difference between transforming a Hessian and merely
// rescaling one.

class BoundedObjective : public Objective {
public:
  BoundedObjective(Objective* inner, std::vector<Coord> coords)
    : inner_(inner), coords_(std::move(coords)) {}

  double value(const arma::vec& eta) {
    const double v = inner_->value(to_theta(coords_, eta));
    sync();
    return v;
  }

  arma::vec gradient(const arma::vec& eta) {
    const arma::vec g = inner_->grad(to_theta(coords_, eta));
    sync();
    arma::vec out(eta.n_elem);
    for (arma::uword i = 0; i < eta.n_elem; ++i) {
      out[i] = g[i] * h_deriv1(coords_[i], eta[i]);
    }
    return out;
  }

  arma::mat hessian(const arma::vec& eta) {
    const arma::vec th = to_theta(coords_, eta);
    const arma::vec g = inner_->grad(th);
    const arma::mat H = inner_->hess(th);
    sync();

    const arma::uword p = eta.n_elem;
    arma::vec d1(p), d2(p);
    for (arma::uword i = 0; i < p; ++i) {
      d1[i] = h_deriv1(coords_[i], eta[i]);
      d2[i] = h_deriv2(coords_[i], eta[i]);
    }

    arma::mat out = arma::diagmat(d1) * H * arma::diagmat(d1);
    out.diag() += g % d2;
    return out;
  }

  // The transform costs nothing that the inner objective does not already do,
  // so what an algorithm may ask for is exactly what the inner objective offers.
  bool has_gradient() const { return true; }
  bool has_hessian() const { return true; }

  // Undo the transform for reporting: the user thinks in theta, and a result
  // reported in eta would be unreadable.
  arma::vec report_par(const arma::vec& eta) const {
    return to_theta(coords_, eta);
  }

  // The gradient of the ORIGINAL objective at the reported point, so that
  // `par` and `gradient` in the result refer to the same thing.
  arma::vec report_grad(const arma::vec& eta, const arma::vec& g_eta) const {
    arma::vec out(eta.n_elem);
    for (arma::uword i = 0; i < eta.n_elem; ++i) {
      const double d1 = h_deriv1(coords_[i], eta[i]);
      out[i] = d1 == 0.0 ? R_NaN : g_eta[i] / d1;
    }
    return out;
  }

private:
  // The counts the user sees must be the inner objective's, since those are the
  // evaluations that cost anything.
  void sync() {
    n_value = inner_->n_value;
    n_grad  = inner_->n_grad;
    n_hess  = inner_->n_hess;
  }

  Objective* inner_;
  std::vector<Coord> coords_;
};


inline std::vector<Coord> make_coords(Rcpp::List bounds) {
  const int p = bounds.size();
  std::vector<Coord> cs(p);
  for (int i = 0; i < p; ++i) {
    Rcpp::NumericVector b = bounds[i];
    Coord c;
    c.lwr = b[0];
    c.upr = b[1];
    const bool lo = R_finite(c.lwr);
    const bool hi = R_finite(c.upr);
    c.kind = lo && hi ? BOUND_BOTH : (lo ? BOUND_LOWER : (hi ? BOUND_UPPER
                                                             : BOUND_NONE));
    cs[i] = c;
  }
  return cs;
}

} // namespace optimizers7

#endif
