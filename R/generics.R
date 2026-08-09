#' @include optimizer_class.R
#' @include result_class.R
#' @include objective.R
NULL

#' @title Minimize a Function
#'
#' @description
#' The entry point of the package. Everything here minimizes; see
#' \code{\link{maximize}} for the other direction.
#'
#' @param optimizer An \code{\link{optimizer}} object, carrying the algorithm and
#'   its settings.
#' @param fn The objective: a function of the parameter vector returning a
#'   single number.
#' @param par A numeric vector of starting values.
#' @param gr An optional gradient function. Ignored when \code{fn} carries its
#'   own gradient.
#' @param he An optional Hessian function. Methods that do not use one ignore
#'   it, so calling code need not branch on the algorithm.
#' @param lower,upper Box constraints. Numeric, of length one — applying to
#'   every parameter — or one value per parameter. Default \code{-Inf} and
#'   \code{Inf}, which is no constraint at all. See Details.
#' @param ... Passed to methods.
#'
#' @details
#' The problem solved is
#'
#' \deqn{\min_{x \in \mathbb{R}^{p}} f(x)
#'   \qquad \text{subject to} \quad l \le x \le u,}
#'
#' with \eqn{f} the objective, \eqn{l} and \eqn{u} the bounds, and the
#' inequalities read coordinatewise. Every method reports a point where it
#' stopped together with the rule that stopped it; convergence is what the
#' stopping rule says and is never inferred from the run having ended.
#'
#' The generic dispatches on \code{optimizer} alone, so each algorithm is
#' written once. The objective is normalized separately by
#' \code{\link{as_objective}}, which dispatches on \code{fn}, so a caller with
#' its own kind of objective registers one method there and every algorithm
#' accepts it.
#'
#' A gradient that is not supplied is computed by central finite differences.
#' The result records which derivatives were supplied and which were
#' differenced, so a run is never silently less exact than it appears. A
#' gradient that is supplied is checked once against the objective -- one
#' central difference along the gradient direction at \code{par}, two
#' evaluations -- and a gross disagreement draws a warning naming both rates,
#' since a \code{gr} computed from a different model than \code{fn} otherwise
#' surfaces as a mute line-search failure at the first iteration;
#' \code{options(optimizers7.check_gradient = FALSE)} disables the check.
#'
#' \subsection{Starters}{
#' \code{par} may be a \strong{starter} rather than a vector:
#' \code{\link{start_zeros}} for all zeros, \code{\link{start_runif}} for a
#' uniform draw from a chosen range. Both work on the \emph{unconstrained}
#' scale and are mapped back through the bounds, which is what makes a single
#' constant sensible for every kind of parameter — zero becomes one for a
#' variance, one half for a probability — and what guarantees the point is
#' admissible however tight the box.
#'
#' A starter needs to know how many parameters there are, and it is told in one
#' of three ways. Say so, with \code{start_zeros(npar = 3)}, and that is the end
#' of it. Otherwise a \code{lower} or \code{upper} with more than one element
#' answers the question, bounds being one per parameter. Otherwise the objective
#' is probed by \code{\link{infer_npar}}, which tries lengths until one is
#' accepted and costs at most fifty evaluations, once, before the run begins.
#'
#' That last route settles any objective with a fixed width built into it,
#' which is most real ones: \code{X \%*\% beta} with a parameter of the wrong
#' length is an error rather than a number. It cannot settle a vectorized toy,
#' since \R recycles a shorter vector silently whenever its length divides, so
#' \code{sum((p - c(1, 2, 3))^2)} is a perfectly finite function of one
#' parameter as well as of three. It then rejects, naming the lengths it found,
#' rather than optimizing a different problem from the one asked.
#' }
#'
#' \code{lower} and \code{upper} are two vectors rather than a list of pairs,
#' which is what \code{\link[stats]{optim}} and \code{\link[stats]{nlminb}}
#' take, and what lets \code{lower = 0} say "every parameter is positive"
#' without writing out one pair per coefficient. Recycling is length one or one
#' per parameter and nothing between, since anything else is far more likely to
#' be a mistake than a request.
#'
#' \strong{Bounds are removed, not enforced.} Each bounded coordinate is
#' reparametrized onto the whole real line — a shifted log for one-sided bounds,
#' a scaled logit for two — and the optimizer runs unconstrained in the new
#' variable. Every point it proposes is admissible by construction, so there is
#' no rejection step and no boundary for a line search to trip over, and any
#' method works with bounds without knowing about them.
#'
#' One limitation follows from the construction: \strong{an optimum
#' lying on a bound cannot be reached}. Getting there requires the transformed
#' variable to run to infinity, so the optimizer marches off, improves by less
#' and less, and stops on a budget at a point merely close to the bound. For the
#' statistical use this exists to serve — a positive variance, a probability
#' inside the unit interval — the optimum is interior and this never arises. For
#' a genuine box-constrained problem with active constraints at the solution,
#' an active-set method is the right tool and this is not one.
#'
#' @return An \code{\link{optimizer_result}}.
#'
#' @examples
#' # a quadratic, with the gradient supplied
#' minimize(gd(), function(p) sum((p - c(1, 2))^2),
#'          par = c(0, 0), gr = function(p) 2 * (p - c(1, 2)))
#'
#' # and without, so the gradient is differenced
#' minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0))
#'
#' # one bound for every parameter: a scale that must stay positive
#' minimize(bfgs(), function(p) sum((p - c(1, 2))^2), c(0.5, 0.5), lower = 0)
#'
#' # or one per parameter. The unconstrained minimum is at (1, 2), so the
#' # second coordinate is pushed against its ceiling of 1.
#' minimize(bfgs(), function(p) sum((p - c(1, 2))^2), c(0.5, 0.5),
#'          lower = c(0, 0), upper = c(5, 1))
#'
#' # no starting value at all: the bounds say there are two parameters
#' minimize(bfgs(), function(p) sum((p - c(1, 2))^2), start_zeros(),
#'          lower = c(0, 0), upper = c(5, 10))
#'
#' @seealso \code{\link{maximize}}, \code{\link{gd}},
#'   \code{\link{start_zeros}}, \code{\link{crit_any}}
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
#' Runs \code{\link{minimize}} on the negated objective and reports the value
#' with its sign restored.
#'
#' @inheritParams minimize
#'
#' @details
#' \deqn{\arg\max_{x} f(x) = \arg\min_{x} \{-f(x)\},
#'   \qquad \max_{x} f(x) = -\min_{x}\{-f(x)\},}
#'
#' so the point is the one \code{\link{minimize}} returns and the value, the
#' gradient and the traced objective have their sign restored.
#'
#' Every algorithm in the package minimizes, always; that is the convention and
#' the names say so. This is the thin wrapper for the other direction, and it
#' exists so that nobody has to remember to negate a log-likelihood by hand and
#' then negate the answer back.
#'
#' @return An \code{\link{optimizer_result}}, whose \code{value} and
#'   \code{gradient} refer to the original objective rather than the negated one.
#'
#' @examples
#' maximize(gd(), function(p) -sum((p - c(1, 2))^2), c(0, 0))
#'
#' @seealso \code{\link{minimize}}
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
