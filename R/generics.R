#' @include optimizer_class.R
#' @include result_class.R
#' @include objective.R
NULL

#' @title Minimize a Function
#'
#' @description
#' Runs an optimizer on an objective and returns the point it reached, the
#' value there, and the rule that stopped it. This is the one entry point:
#' every algorithm in the package is a method of this generic, so changing
#' method means changing the first argument and nothing else. Box constraints
#' are accepted by every method. Everything here minimizes; [maximize()]
#' negates the objective and restores the sign of the answer.
#'
#' @details
#' The problem solved is
#'
#' \deqn{\min_{x \in \mathbb{R}^{p}} f(x)
#'   \qquad \text{subject to} \quad l \le x \le u,}
#'
#' with \eqn{f} the objective, \eqn{l} and \eqn{u} the bounds, and the
#' inequalities read coordinatewise. Every method reports the point where it
#' stopped together with the rule that stopped it. Convergence is what a
#' stopping rule confirmed and is never inferred from the run having ended,
#' so a run that exhausts its budget comes back with `converged = FALSE`.
#'
#' Dispatch is on `optimizer` alone, so each algorithm is written once. The
#' objective is normalized separately, by [as_objective()], which dispatches
#' on `fn`; a caller with its own kind of objective registers one method
#' there and every algorithm accepts it.
#'
#' # Derivatives, supplied and differenced
#'
#' A gradient that is not supplied is computed by central finite differences,
#' at a cost of \eqn{2p} evaluations of the objective each time. The result's
#' `message` records that, so a run is never silently less exact than it
#' looks.
#'
#' A gradient that *is* supplied is checked once, before the run: one central
#' difference along the gradient direction at `par`, two evaluations. A gross
#' disagreement draws a warning naming both rates, as in
#'
#' ```
#' 'gr' does not appear to be the gradient of 'fn': along the gradient
#' direction at 'par', 'fn' changes at rate 4.47 where 'gr' predicts 22.4.
#' ```
#'
#' and the run proceeds. A `gr` computed from a different model than `fn` is
#' otherwise very hard to see: it surfaces as a mute line-search failure at
#' the first iteration. Set `options(optimizers7.check_gradient = FALSE)` to
#' turn the check off.
#'
#' # Starters
#'
#' `par` may be a **starter** object instead of a vector: [start_zeros()] for
#' all zeros, [start_runif()] for a uniform draw from a chosen range. Both
#' work on the unconstrained scale and are mapped back through the bounds, so
#' one constant means something sensible for every kind of parameter (zero
#' becomes one for a variance, one half for a probability) and no draw can
#' land outside its box.
#'
#' A starter has to be told how many parameters there are, in one of three
#' ways. Say so with `start_zeros(npar = 3)`. Failing that, a `lower` or
#' `upper` of more than one element answers the question, bounds being one
#' per parameter. Failing that, the objective is probed by [infer_npar()],
#' which tries lengths 1 to 50 until one is accepted, once, before the run
#' begins.
#'
#' That last route settles any objective with a fixed width built into it,
#' which is most real ones: `X %*% beta` with a parameter of the wrong length
#' is an error. It cannot settle a vectorized toy, since \R recycles a shorter
#' vector silently whenever its length divides, so `sum((p - c(1, 2, 3))^2)`
#' is a perfectly finite function of one parameter as well as of three. It
#' then rejects, naming the lengths it found, and asks for `npar`.
#'
#' # Bounds are removed, not enforced
#'
#' Each bounded coordinate is reparametrized onto the whole real line, by a
#' shifted log for a one-sided bound and a scaled logit for a two-sided one,
#' and the optimizer runs unconstrained in the new variable. Every point it
#' proposes is admissible by construction, so there is no rejection step and
#' no boundary for a line search to trip over, and any method takes bounds
#' without knowing they exist. [bounded_transform()] is the map.
#'
#' `lower` and `upper` are two vectors, as in [stats::optim()] and
#' [stats::nlminb()], so `lower = 0` says that every parameter is positive
#' without writing out one pair per coefficient. A length that is neither 1
#' nor `p` is refused, and so is a `lower` at or above its `upper`. The
#' **starting value must lie strictly inside its bounds**: `par = 0` with
#' `lower = 0` is an error, the transformed coordinate being infinite there.
#'
#' # An optimum lying on a bound
#'
#' Reaching a bound exactly requires the transformed variable to run to
#' infinity, so a solution on a bound is approached and not attained. What
#' stops the run is the stopping rule, because the chain-rule factor
#' \eqn{\partial x / \partial \eta} decays to zero as the bound is
#' approached and takes the transformed gradient with it. The run therefore
#' reports `converged = TRUE` at a point near the bound, and how near is set
#' by the tolerance rather than by the budget.
#'
#' Measured on \eqn{\sum (x - (1, 2))^2} with `upper = c(5, 1)`, whose
#' second coordinate wants to sit at its ceiling: the gap to the bound is
#' `3.6e-07` at `crit_grad(1e-6)`, `1.1e-11` at `1e-10` and `2.7e-15` at
#' `1e-14`, and raising `maxit` from 50 to 5000 changes nothing, all three
#' runs stopping after 22 iterations.
#'
#' For the statistical use this exists to serve, a positive variance or a
#' probability inside the unit interval, the optimum is interior and none of
#' this arises. For a genuine box-constrained problem whose constraints are
#' active at the solution, an active-set method is the right tool.
#'
#' @param optimizer An [optimizer()] object carrying the algorithm and its
#'   settings. The only argument dispatch reads.
#' @param fn The objective: a function of the parameter vector returning a
#'   single number, to be minimized. An object of another class works
#'   whenever [as_objective()] has a method for it.
#' @param par A numeric vector of starting values, or a starter object; see
#'   Starters below. With bounds it must lie **strictly** inside them.
#' @param gr The gradient, a function of the parameter vector. `NULL`, the
#'   default, has it differenced from `fn`. Ignored when `fn` is an objective
#'   carrying its own gradient.
#' @param he The Hessian, a function of the parameter vector. `NULL` is the
#'   default and only [newton()] reads one; every other method accepts it and
#'   ignores it, so calling code need not branch on the algorithm.
#' @param lower,upper Box constraints, numeric of length 1 (applying to every
#'   parameter) or of length `p` (one each). Any other length is an error
#'   naming both. The defaults `-Inf` and `Inf` are no constraint, and a
#'   coordinate whose pair is infinite on both sides is left untransformed.
#'   `lower` must be strictly below `upper` in every coordinate.
#' @param ... Passed to the method dispatched on. No shipped method reads
#'   anything from it.
#'
#' @return An [optimizer_result()] object: `par`, `value`, `gradient`,
#'   `counts`, `iterations`, `converged`, `criterion_met`, `message`,
#'   `trace`, `optimizer`, `elapsed` and `seed`.
#'
#' @examples
#' q <- function(p) sum((p - c(1, 2))^2)
#' qg <- function(p) 2 * (p - c(1, 2))
#'
#' # With the gradient supplied, and without, so it is differenced. The two
#' # reach the same point; only the evaluation counts differ.
#' a <- minimize(bfgs(), q, par = c(0, 0), gr = qg)
#' b <- minimize(bfgs(), q, c(0, 0))
#' all.equal(a@par, b@par, tolerance = 1e-6)
#' rbind(supplied = unlist(a@counts), differenced = unlist(b@counts))
#'
#' # One bound for every parameter: a scale that must stay positive.
#' minimize(bfgs(), q, c(0.5, 0.5), lower = 0)@par
#'
#' # Or one per parameter. The unconstrained minimum is at (1, 2), so the
#' # second coordinate is pushed against its ceiling of 1 and stops just
#' # short of it, by an amount the tolerance sets.
#' minimize(bfgs(), q, c(0.5, 0.5), lower = c(0, 0), upper = c(5, 1))@par
#' minimize(bfgs(criterion = crit_grad(1e-12)), q, c(0.5, 0.5),
#'          lower = c(0, 0), upper = c(5, 1))@par
#'
#' # No starting value at all: the bounds say there are two parameters.
#' minimize(bfgs(), q, start_zeros(), lower = c(0, 0), upper = c(5, 10))@par
#'
#' # A start sitting on its own bound is refused: the transform is infinite
#' # there.
#' try(minimize(bfgs(), q, c(0, 0), lower = 0))
#'
#' # A gradient belonging to a different model draws a warning naming both
#' # rates, before the run, and the run then proceeds.
#' wrong <- minimize(gd(maxit = 2), q, c(0, 0), gr = function(p) 5 * qg(p))
#' wrong@converged
#'
#' @seealso [maximize()] for the other direction, [optimizer()] for the
#'   algorithms, [criterion()] for the stopping rules, [start_zeros()] for
#'   the starters, [bounded_transform()] for the map that removes the bounds.
#' @export
minimize <- S7::new_generic("minimize", "optimizer",
  function(optimizer, fn, par, gr = NULL, he = NULL,
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
#' written. Every algorithm in the package minimizes; this is the wrapper for
#' the other direction, so that a log-likelihood need not be negated by hand
#' and the answer negated back.
#'
#' @details
#' \deqn{\arg\max_{x} f(x) = \arg\min_{x} \{-f(x)\},
#'   \qquad \max_{x} f(x) = -\min_{x}\{-f(x)\},}
#'
#' so `par` is exactly what [minimize()] returned, while `value`, `gradient`
#' and the `value` column of `trace` have their sign restored. Everything
#' else on the result belongs to the run and passes through untouched:
#' `counts`, `iterations`, `converged`, `criterion_met` and the optimizer.
#'
#' A consequence worth knowing: `criterion_met` names a rule that was
#' evaluated on the **negated** objective. For a rule reading a gradient norm
#' or an absolute change that makes no difference, both being invariant under
#' the sign. For a rule reading the objective's own scale, such as
#' [crit_rel_obj()], it makes none either, the relative change being a ratio.
#'
#' `fn` here must be a plain function. An objective of some other class,
#' reached through a method of [as_objective()], is refused by name: this
#' wrapper negates by composing a closure, and it has no way to negate an
#' object whose evaluation it does not perform.
#'
#' @param optimizer An [optimizer()] object carrying the algorithm and its
#'   settings.
#' @param fn The objective, a **plain function** of the parameter vector
#'   returning a single number, to be maximized. Unlike [minimize()], an
#'   object of another class is refused rather than passed to
#'   [as_objective()]; negate such an objective yourself and call
#'   [minimize()].
#' @param par A numeric vector of starting values, or a starter object such
#'   as [start_zeros()]. With bounds it must lie strictly inside them.
#' @param gr The gradient of `fn` as written, a function of the parameter
#'   vector, or `NULL` for a central difference. It is negated here, so a
#'   caller supplies the gradient of the function being maximized.
#' @param he The Hessian of `fn` as written, or `NULL`. Negated here too, and
#'   read only by [newton()].
#' @param lower,upper Box constraints, numeric of length 1 or of length `p`,
#'   defaulting to `-Inf` and `Inf`. Passed through to [minimize()]
#'   unchanged; the sign convention does not touch them.
#' @param ... Passed to [minimize()], and from there to the method.
#'
#' @return An [optimizer_result()] with the same twelve properties
#'   [minimize()] returns, and with `value`, `gradient` and `trace$value`
#'   referring to `fn` as written.
#'
#' @examples
#' # A log-likelihood: the normal mean, maximized directly.
#' set.seed(1)
#' y <- rnorm(200, mean = 3)
#' ll <- function(p) sum(dnorm(y, mean = p, log = TRUE))
#' fit <- maximize(bfgs(), ll, par = 0)
#' c(fit@par, mean(y))
#'
#' # The reported value is the log-likelihood itself, not its negative.
#' all.equal(fit@value, ll(fit@par))
#'
#' # An objective that is not a plain function is refused by name.
#' try(maximize(bfgs(), 1:3, c(0, 0)))
#'
#' @seealso [minimize()]
#' @export
maximize <- function(optimizer, fn, par, gr = NULL, he = NULL,
                     lower = -Inf, upper = Inf, ...) {
  if (!is.function(fn)) {
    stop("maximize() takes a plain function; negate other objectives yourself.",
         call. = FALSE)
  }
  neg_fn <- function(p) -fn(p)
  neg_gr <- if (is.null(gr)) NULL else function(p) -gr(p)
  neg_he <- if (is.null(he)) NULL else function(p) -he(p)

  res <- minimize(optimizer, neg_fn, par, gr = neg_gr, he = neg_he,
                  lower = lower, upper = upper, ...)

  res@value <- -res@value
  if (!is.null(res@gradient)) res@gradient <- -res@gradient
  if (!is.null(res@trace)) res@trace$value <- -res@trace$value
  res
}
