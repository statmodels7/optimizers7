#ifndef OPTIMIZERS7_RULES_H
#define OPTIMIZERS7_RULES_H

// The shipped stopping rules of the descent loop, compiled.
//
// Each class computes what its crit_met() method computes in R, with the same
// arithmetic, so that a loop asking a compiled rule stops at the same
// iteration as one asking the R rule: crit_grad(), crit_abs_obj(),
// crit_rel_obj(), crit_abs_par(), crit_rel_par(), crit_never(), crit_any() and
// crit_all(). R's comparisons return NA where an operand is NaN, and any() and
// all() combine NA by R's rules, so each rule returns one of three values and
// only TRUE stops the run, as ask_criterion() reads an R answer. A rule whose
// class is not one of these (a user's subclass) has no compiled form, and
// criterion_spec() on the R side returns NULL for it.

#include <RcppArmadillo.h>
#include <cmath>
#include <memory>
#include <vector>
#include "descent_loop.h"

namespace optimizers7 {

enum class Tri { False = 0, True = 1, NA = 2 };

inline Tri tri_less(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return Tri::NA;
  return a < b ? Tri::True : Tri::False;
}

// max() of R over doubles: NA where any element is NaN.
inline double r_max(const arma::vec& v) {
  double m = -std::numeric_limits<double>::infinity();
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    if (std::isnan(v[i])) return NA_REAL;
    if (v[i] > m) m = v[i];
  }
  return m;
}

struct RuleArgs {
  int iter;
  double f_new, f_old;
  const arma::vec& x_new;
  const arma::vec& x_old;
  const arma::vec& g;
  bool have_g, have_old;
};

class CompiledRule {
public:
  virtual ~CompiledRule() {}
  virtual Tri eval(const RuleArgs& a) const = 0;
};

// crit_grad(): FALSE without a gradient or with an NA in it; the max-norm or
// the 2-norm, the latter as sqrt(sum(g^2)) with sum() accumulating in long
// double as R does.
class RuleGrad : public CompiledRule {
public:
  RuleGrad(double tol, bool two) : tol_(tol), two_(two) {}
  Tri eval(const RuleArgs& a) const {
    if (!a.have_g || a.g.n_elem == 0) return Tri::False;
    for (arma::uword i = 0; i < a.g.n_elem; ++i) {
      if (std::isnan(a.g[i])) return Tri::False;
    }
    double v;
    if (two_) {
      long double s = 0.0;
      for (arma::uword i = 0; i < a.g.n_elem; ++i) {
        const double q = a.g[i] * a.g[i];
        s += q;
      }
      v = std::sqrt((double) s);
    } else {
      v = r_max(arma::abs(a.g));
    }
    return tri_less(v, tol_);
  }
private:
  double tol_;
  bool two_;
};

// crit_abs_obj() and crit_rel_obj(): FALSE without a previous value or with a
// non-finite one.
class RuleObj : public CompiledRule {
public:
  RuleObj(double tol, bool rel) : tol_(tol), rel_(rel) {}
  Tri eval(const RuleArgs& a) const {
    if (!a.have_old || !std::isfinite(a.f_old)) return Tri::False;
    const double d = std::fabs(a.f_new - a.f_old);
    if (rel_) return tri_less(d, tol_ * (std::fabs(a.f_old) + tol_));
    return tri_less(d, tol_);
  }
private:
  double tol_;
  bool rel_;
};

// crit_abs_par() and crit_rel_par(): FALSE without a previous iterate.
class RulePar : public CompiledRule {
public:
  RulePar(double tol, bool rel) : tol_(tol), rel_(rel) {}
  Tri eval(const RuleArgs& a) const {
    if (!a.have_old) return Tri::False;
    arma::vec d = arma::abs(a.x_new - a.x_old);
    if (rel_) d = d / (arma::abs(a.x_old) + tol_);
    return tri_less(r_max(d), tol_);
  }
private:
  double tol_;
  bool rel_;
};

class RuleNever : public CompiledRule {
public:
  Tri eval(const RuleArgs&) const { return Tri::False; }
};

// crit_any() and crit_all(): every member is evaluated, as vapply() does, and
// the answers combined by R's any() and all().
class RuleCombine : public CompiledRule {
public:
  RuleCombine(std::vector<std::unique_ptr<CompiledRule>> rs, bool all)
    : rs_(std::move(rs)), all_(all) {}
  Tri eval(const RuleArgs& a) const {
    bool any_true = false, any_false = false, any_na = false;
    for (const auto& r : rs_) {
      const Tri t = r->eval(a);
      if (t == Tri::True) any_true = true;
      else if (t == Tri::False) any_false = true;
      else any_na = true;
    }
    if (all_) {
      if (any_false) return Tri::False;
      return any_na ? Tri::NA : Tri::True;
    }
    if (any_true) return Tri::True;
    return any_na ? Tri::NA : Tri::False;
  }
private:
  std::vector<std::unique_ptr<CompiledRule>> rs_;
  bool all_;
};

// Built from the list criterion_spec() returns on the R side: `type` one of
// "grad", "abs_obj", "rel_obj", "abs_par", "rel_par", "never", "any", "all";
// `tol` and `norm` where the type has them; `criteria` for a combination.
inline std::unique_ptr<CompiledRule> make_compiled_rule(Rcpp::List spec) {
  const std::string type = Rcpp::as<std::string>(spec["type"]);
  if (type == "grad") {
    return std::unique_ptr<CompiledRule>(new RuleGrad(
      Rcpp::as<double>(spec["tol"]), Rcpp::as<std::string>(spec["norm"]) == "2"));
  }
  if (type == "abs_obj" || type == "rel_obj") {
    return std::unique_ptr<CompiledRule>(new RuleObj(
      Rcpp::as<double>(spec["tol"]), type == "rel_obj"));
  }
  if (type == "abs_par" || type == "rel_par") {
    return std::unique_ptr<CompiledRule>(new RulePar(
      Rcpp::as<double>(spec["tol"]), type == "rel_par"));
  }
  if (type == "never") return std::unique_ptr<CompiledRule>(new RuleNever());
  if (type == "any" || type == "all") {
    Rcpp::List cs = spec["criteria"];
    std::vector<std::unique_ptr<CompiledRule>> rs;
    for (R_xlen_t i = 0; i < cs.size(); ++i) {
      rs.push_back(make_compiled_rule(Rcpp::List(cs[i])));
    }
    return std::unique_ptr<CompiledRule>(new RuleCombine(std::move(rs), type == "all"));
  }
  Rcpp::stop("Unknown compiled rule '%s'.", type);
  return nullptr;
}

// A compiled rule as the loop asks it. A caller may also hold a rule that it
// replaces between runs (a resolution that moves): `rule` is then reset.
class CompiledStopRule : public StopRule {
public:
  explicit CompiledStopRule(std::unique_ptr<CompiledRule> r) : rule_(std::move(r)) {}
  bool met(int iter, double f_new, double f_old,
           const arma::vec& x_new, const arma::vec& x_old,
           const arma::vec& g, bool have_g, bool have_old) {
    const RuleArgs a{iter, f_new, f_old, x_new, x_old, g, have_g, have_old};
    return rule_->eval(a) == Tri::True;
  }
  const CompiledRule& rule() const { return *rule_; }
private:
  std::unique_ptr<CompiledRule> rule_;
};

} // namespace optimizers7

#endif
