#' @include optimizer_class.R
#' @include result_class.R
#' @include objective.R
NULL

#' @title Minimize a Function
#'
#' @description
#' Runs an optimizer on an objective and returns the point that it reached,
#' the value there, and the rule that stopped it. This is the single entry
#' point: every algorithm in the package is a method of this generic, so a
#' change of method is a change of the first argument. Box constraints are
#' accepted by every method except [prox_grad()], which takes its constraint
#' through the proximal operator, and [chebyshev()] when it is given `t3`.
#' Every method minimizes; [maximize()] negates the objective and restores
#' the sign of the answer.
#'
#' @details
#' The problem solved is
#'
#' \deqn{\min_{x \in \mathbb{R}^{p}} f(x)
#'   \qquad \text{subject to} \quad l \le x \le u,}
#'
#' with \eqn{f} the objective, \eqn{l} and \eqn{u} the bounds, and the
#' inequalities read coordinatewise. Every method reports the point where it
#' stopped together with the rule that stopped it. `converged` is `TRUE`
#' when the stopping rule is satisfied; it is also `TRUE` when a line search
#' with a declared resolution establishes that every step fails to decrease
#' the objective by more than that resolution, and when [brent()] reduces its
#' bracket to the resolution of the iterate, and the result then reports
#' that stop with its own label in `criterion_met`. It is never inferred from
#' the end of the run, so a run that exhausts its budget comes back with
#' `converged = FALSE`.
#'
#' Dispatch is on `optimizer` alone, so each algorithm is written once. The
#' objective is normalized separately, by [as_objective()], which dispatches
#' on `fn`; a caller with its own kind of objective registers one method
#' there and every algorithm accepts it.
#'
#' # Derivatives, supplied and differenced
#'
#' A gradient that is not supplied is computed by central finite differences,
#' at a cost of \eqn{2p} evaluations of the objective each time, and the
#' `message` of the result records that the gradient was obtained by finite
#' differences.
#'
#' A gradient that *is* supplied is checked before the run, with one central
#' difference along the gradient direction at `par`, which costs two
#' evaluations. If the difference disagrees with the gradient, it is repeated
#' at steps ten and a hundred times smaller, and a disagreement that persists
#' gives a warning that names both rates, as in
#'
#' ```
#' 'gr' does not appear to be the gradient of 'fn': along the gradient
#' direction at 'par', 'fn' changes at rate 4.47 where 'gr' predicts 22.4.
#' ```
#'
#' and the run proceeds. Without the check, a `gr` computed from a different
#' model than `fn` appears only as a line-search failure at the first
#' iteration, with no indication of its cause. Set
#' `options(optimizers7.check_gradient = FALSE)` to turn the check off.
#'
#' # Starters
#'
#' `par` may be a starter object instead of a vector: [start_zeros()] for
#' all zeros, [start_runif()] for a uniform draw from a chosen range. Both
#' work on the unconstrained scale and are mapped back through the bounds, so
#' one constant gives a sensible value for every kind of parameter (zero
#' becomes one for a variance, one half for a probability) and every draw
#' lies inside its box.
#'
#' The number of parameters is obtained in one of three ways. It can be given
#' with `start_zeros(npar = 3)`. Otherwise a `lower` or `upper` of more than
#' one element gives it, since bounds are one per parameter. Otherwise the
#' objective is probed by [infer_npar()], which tries the lengths 1 to 50
#' once, before the run begins, and succeeds when exactly one of them is
#' accepted.
#'
#' The probe works for an objective with a fixed width built into it, such
#' as `X %*% beta`, which raises an error for a parameter of the wrong
#' length. It does not work for a vectorized expression, since \R recycles a
#' shorter vector silently whenever its length divides, so
#' `sum((p - c(1, 2, 3))^2)` is a finite function of one parameter as well as
#' of three. The probe then signals an error that names the lengths found,
#' and `npar` has to be given.
#'
#' # Bounds by reparametrization
#'
#' Each bounded coordinate is reparametrized onto the whole real line, by a
#' shifted log for a one-sided bound and a scaled logit for a two-sided one,
#' and the optimizer runs unconstrained in the new variable. Every point that
#' it proposes is admissible by construction, so there is no rejection step
#' and no boundary for a line search to cross, and a method needs no code of
#' its own for the bounds. [bounded_transform()] is the map.
#'
#' `lower` and `upper` are two vectors, as in [stats::optim()] and
#' [stats::nlminb()], so `lower = 0` applies the lower bound zero to every
#' parameter without one value per coefficient. A length that is neither 1
#' nor `p` is rejected, and so is a `lower` that is at or above its `upper`.
#' The starting value must lie strictly inside its bounds: `par = 0` with
#' `lower = 0` is an error, because the transformed coordinate is infinite
#' there.
#'
#' # An optimum lying on a bound
#'
#' Reaching a bound exactly requires the transformed variable to run to
#' infinity, so a solution on a bound is approached and not attained. The
#' stopping rule ends the run, because the chain-rule factor
#' \eqn{\partial x / \partial \eta} decays to zero as the bound is
#' approached and takes the transformed gradient with it. The run therefore
#' reports `converged = TRUE` at a point near the bound, and how near is set
#' by the tolerance rather than by the budget.
#'
#' The bounds examples below show this on \eqn{\sum (x - (1, 2))^2} with
#' `upper = c(5, 1)`, whose minimum over the box lies on the upper bound of
#' the second coordinate: a smaller tolerance of the rule gives a point
#' closer to the bound, after more iterations.
#'
#' For the parameters for which the bounds are intended, such as a positive
#' variance or a probability inside the unit interval, the optimum is
#' usually interior and this does not arise. For a box-constrained problem
#' whose constraints are active at the solution, an active-set method is
#' more suitable.
#'
#' @param optimizer An [optimizer()] object carrying the algorithm and its
#'   settings. Dispatch reads this argument only.
#' @param fn The objective: a function of the parameter vector returning a
#'   single number, to be minimized. An object of another class works
#'   whenever [as_objective()] has a method for it.
#' @param par A numeric vector of starting values, or a starter object; see
#'   Starters below. With bounds it must lie strictly inside them.
#' @param gr The gradient, a function of the parameter vector. `NULL`, the
#'   default, has it differenced from `fn`. It is ignored when `fn` is an
#'   objective that carries its own gradient.
#' @param hess The Hessian, a function of the parameter vector. `NULL` is the
#'   default. [newton()] and [chebyshev()] read it; every other method
#'   accepts it and does not use it, so calling code need not branch on the
#'   algorithm.
#' @param lower,upper Box constraints, numeric of length 1 (applying to every
#'   parameter) or of length `p` (one each). Any other length is an error.
#'   The defaults `-Inf` and `Inf` are no constraint, and a
#'   coordinate whose pair is infinite on both sides is left untransformed.
#'   `lower` must be strictly below `upper` in every coordinate.
#' @param ... Passed to the method selected by dispatch. Among the shipped
#'   methods, only the method for [chebyshev()] reads from it, its argument
#'   `t3`.
#'
#' @return An [optimizer_result()] object: `par`, `value`, `gradient`,
#'   `counts`, `iterations`, `converged`, `criterion_met`, `message`,
#'   `trace`, `optimizer`, `elapsed` and `seed`.
#'
#' @examples
#' q <- function(p) sum((p - c(1, 2))^2)
#' qg <- function(p) 2 * (p - c(1, 2))
#'
#' # The first run supplies the gradient, and the second differences it.
#' # The two runs reach the same point and differ in the evaluation counts.
#' a <- minimize(bfgs(), q, par = c(0, 0), gr = qg)
#' b <- minimize(bfgs(), q, c(0, 0))
#' all.equal(a@par, b@par, tolerance = 1e-6)
#' rbind(supplied = unlist(a@counts), differenced = unlist(b@counts))
#'
#' # One bound applies to every parameter: a scale that must stay positive.
#' minimize(bfgs(), q, c(0.5, 0.5), lower = 0)@par
#'
#' # Bounds may also be given per parameter. The unconstrained minimum is at (1, 2), so the
#' # second coordinate approaches its upper bound of 1 and stops just short of
#' # it, at a distance set by the tolerance.
#' minimize(bfgs(), q, c(0.5, 0.5), lower = c(0, 0), upper = c(5, 1))@par
#' minimize(bfgs(criterion = crit_grad(1e-12)), q, c(0.5, 0.5),
#'          lower = c(0, 0), upper = c(5, 1))@par
#'
#' # Without a starting value, the bounds give the number of parameters.
#' minimize(bfgs(), q, start_zeros(), lower = c(0, 0), upper = c(5, 10))@par
#'
#' # A start on its own bound is rejected: the transform is infinite there.
#' try(minimize(bfgs(), q, c(0, 0), lower = 0))
#'
#' # A gradient that belongs to a different model gives a warning, before the
#' # run, that names both rates, and the run then proceeds.
#' wrong <- minimize(gd(maxit = 2), q, c(0, 0), gr = function(p) 5 * qg(p))
#' wrong@converged
#'
#' @seealso [maximize()] for the other direction, [optimizer()] for the
#'   algorithms, [criterion()] for the stopping rules, [start_zeros()] for
#'   the starters, [bounded_transform()] for the map that removes the bounds.
#' @export
minimize <- S7::new_generic("minimize", "optimizer",
  function(optimizer, fn, par, gr = NULL, hess = NULL,
           lower = -Inf, upper = Inf, ...) {
    # A starter becomes an ordinary numeric vector here, before dispatch, so
    # that every method -- including one a user wrote -- receives what it always
    # received and needs to know nothing about starters. A numeric `par` is
    # returned untouched, so the ordinary call pays nothing for this.
    par <- resolve_start(par, fn, gr, lower, upper)
    # One central difference against the supplied gradient, here for the same
    # reason: every method gets the guard without knowing it exists.
    if (!is.null(gr)) check_gradient_consistency(fn, gr, par)
    S7::S7_dispatch()
  })


#' @title Maximize a Function
#'
#' @description
#' Runs [minimize()] on the negated objective and hands back a result whose
#' value, gradient and traced objective are those of the objective as
#' written. Every algorithm in the package minimizes; this wrapper maximizes,
#' so that a log-likelihood need not be negated by hand and the answer
#' negated back.
#'
#' @details
#' \deqn{\arg\max_{x} f(x) = \arg\min_{x} \{-f(x)\},
#'   \qquad \max_{x} f(x) = -\min_{x}\{-f(x)\},}
#'
#' so `par` is exactly what [minimize()] returned, while `value`, `gradient`
#' and the `value` column of `trace` have their sign restored. The other
#' properties of the result pass through unchanged: `counts`, `iterations`,
#' `converged`, `criterion_met` and the optimizer.
#'
#' The rule named in `criterion_met` was evaluated on the negated
#' objective. Each shipped rule gives the same verdict on \eqn{f} and on
#' \eqn{-f}, since each reads a gradient norm, an absolute change or a
#' relative change.
#'
#' `fn` here must be a plain function. An objective of some other class,
#' reached through a method of [as_objective()], is rejected with an error:
#' this wrapper negates by composing a closure, and it cannot negate an
#' object whose evaluation it does not perform.
#'
#' @param optimizer An [optimizer()] object carrying the algorithm and its
#'   settings.
#' @param fn The objective, a plain function of the parameter vector
#'   returning a single number, to be maximized. Unlike [minimize()], an
#'   object of another class is rejected and not passed to
#'   [as_objective()]; such an objective has to be negated by the caller and
#'   passed to [minimize()].
#' @param par A numeric vector of starting values, or a starter object such
#'   as [start_zeros()]. With bounds it must lie strictly inside them.
#' @param gr The gradient of `fn` as written, a function of the parameter
#'   vector, or `NULL` for a central difference. It is negated here, so a
#'   caller supplies the gradient of the function being maximized.
#' @param hess The Hessian of `fn` as written, or `NULL`. It is negated here
#'   too, and only [newton()] and [chebyshev()] read it.
#' @param lower,upper Box constraints, numeric of length 1 or of length `p`,
#'   defaulting to `-Inf` and `Inf`. Passed through to [minimize()]
#'   unchanged; the sign convention does not touch them.
#' @param ... Passed to [minimize()], and from there to the method. A `t3`
#'   argument for [chebyshev()] is the contraction of the third derivative of
#'   `fn` as written, and is negated here like `gr` and `hess`.
#'
#' @return An [optimizer_result()] with the same twelve properties that
#'   [minimize()] returns, in which `value`, `gradient` and `trace$value`
#'   refer to `fn` as written.
#'
#' @examples
#' # A log-likelihood for the normal mean is maximized directly.
#' set.seed(1)
#' y <- rnorm(200, mean = 3)
#' ll <- function(p) sum(dnorm(y, mean = p, log = TRUE))
#' fit <- maximize(bfgs(), ll, par = 0)
#' c(fit@par, mean(y))
#'
#' # The reported value is the log-likelihood itself, not its negative.
#' all.equal(fit@value, ll(fit@par))
#'
#' # An objective that is not a plain function is rejected.
#' try(maximize(bfgs(), 1:3, c(0, 0)))
#'
#' @seealso [minimize()]
#' @export
maximize <- function(optimizer, fn, par, gr = NULL, hess = NULL,
                     lower = -Inf, upper = Inf, ...) {
  if (!is.function(fn)) {
    stop("maximize() takes a plain function; negate other objectives yourself.",
         call. = FALSE)
  }
  neg_fn <- function(p) -fn(p)
  neg_gr <- if (is.null(gr)) NULL else function(p) -gr(p)
  neg_he <- if (is.null(hess)) NULL else function(p) -hess(p)
  # the third-derivative contraction of chebyshev() arrives through `...` and
  # is a derivative of fn like the other two, so it is negated as well
  dots <- list(...)
  if (is.function(dots[["t3"]])) {
    t3 <- dots[["t3"]]
    dots[["t3"]] <- function(x, d) -t3(x, d)
  }

  res <- do.call(minimize, c(list(optimizer, neg_fn, par, gr = neg_gr,
                                  hess = neg_he, lower = lower, upper = upper),
                             dots))

  res@value <- -res@value
  if (!is.null(res@gradient)) res@gradient <- -res@gradient
  if (!is.null(res@trace)) res@trace$value <- -res@trace$value
  res
}
