#' @include criterion.R
NULL

# The objective, and the two shapes it may arrive in.
#
# minimize() dispatches on the optimizer, so an algorithm is written once.
# as_objective() dispatches on the objective, so the two shapes are told apart
# once. Each generic dispatches where there is real variation; dispatching
# minimize() on both would need one method per algorithm per shape, and every
# algorithm would be written twice.
#
# There was a third, finite_sum(), an objective declaring itself a sum over
# observations so that a stochastic method could ask it for a subsample. It
# existed for exactly one caller, Adam, and went when Adam stopped drawing its
# own minibatches: an objective that resamples is a closure, and needs no class.


#' @title Normalize an Objective for the Optimizers
#'
#' @description
#' Packs a function to be minimized, together with whichever of its gradient
#' and Hessian the caller has, into the one handle every algorithm in the
#' package is written against. The result is a named list carrying the pieces
#' and two flags saying which derivatives were supplied. `as_objective()` is
#' a generic dispatching on `fn`, so a caller holding some other kind of
#' objective registers one method and every algorithm accepts it.
#'
#' @details
#' # Two generics, each dispatching once
#'
#' [minimize()] dispatches on the optimizer, so an algorithm is written once
#' whatever shape the objective arrived in. `as_objective()` dispatches on
#' the objective, so the shapes are told apart once whatever algorithm is
#' running. Dispatching [minimize()] on both would need one method per
#' algorithm per shape.
#'
#' # What happens when no gradient is supplied
#'
#' The compiled loop differences the objective, once per coordinate,
#' \deqn{\hat{g}_j = \frac{f(x + h e_j) - f(x - h e_j)}{2h}, \qquad
#'       h = \varepsilon^{1/3} \max(1, \lvert x_j \rvert),}
#' with \eqn{\varepsilon} the machine epsilon, so \eqn{h} is about
#' `6.06e-06` for a coordinate of size one. One gradient then costs
#' \eqn{2p} evaluations of the objective. The accuracy is what a central
#' difference gives: on the Rosenbrock function at \eqn{(0.5, 0.5)} the
#' differenced gradient is right to a relative `1.4e-10`.
#'
#' A run that differences its gradient says so in the result's `message`
#' field, so a fit is never quietly less exact than it looks.
#'
#' # The Hessian is the one place two differences compose
#'
#' Only [newton()] asks for a Hessian, and when none is supplied it is
#' obtained by differencing the gradient once. With an analytic gradient
#' that is a single differentiation and costs \eqn{2p} gradient evaluations.
#' With no gradient either, the gradient is itself a difference and the
#' Hessian is a difference of differences, which costs about \eqn{4p^2}
#' evaluations and loses accuracy accordingly: one Newton iteration on a
#' quadratic in 20 unknowns takes 1682 evaluations of the objective without
#' a gradient and 42 of the gradient with one. The Hessian is symmetrized
#' before it is returned, the two triangles differing by the differencing
#' error alone.
#'
#' # A point the objective refuses
#'
#' An objective that returns `NA`, `NaN` or an infinity at some point
#' declares that point unusable. The algorithms treat that as an ordinary
#' event, back away and continue, so a line search may probe outside the
#' region where the objective is defined without the run failing.
#'
#' @param fn The objective, a function of the parameter vector returning a
#'   single number to be **minimized**. An object of some other class is
#'   accepted whenever a method has been registered for it; without one,
#'   dispatch fails and reports the class it was handed.
#' @param gr The gradient, a function of the parameter vector returning a
#'   numeric vector of the same length. `NULL`, the default, has the
#'   gradient differenced from `fn` as above. Anything that is neither a
#'   function nor `NULL` is an error.
#' @param he The Hessian, a function of the parameter vector returning a
#'   `p x p` symmetric matrix. `NULL` is the default. Only [newton()]
#'   reads one; every other method accepts it and ignores it, so calling
#'   code can pass whatever it has without branching on the method.
#' @param ... Passed to the method dispatched on. The shipped method for an
#'   ordinary function reads nothing from it.
#'
#' @return A named list of six components describing the objective to the
#'   compiled side:
#'   \describe{
#'     \item{`kind`}{character, `"r"` for an ordinary R function.}
#'     \item{`fn`}{the objective as supplied.}
#'     \item{`gr`}{the gradient, or `NULL`.}
#'     \item{`he`}{the Hessian, or `NULL`.}
#'     \item{`has_gradient`}{logical, whether `gr` was supplied.}
#'     \item{`has_hessian`}{logical, whether `he` was supplied.}
#'   }
#'
#' @examples
#' # The flags record what was supplied, and they are what the loop branches on.
#' str(as_objective(function(p) sum(p^2)))
#' str(as_objective(function(p) sum(p^2), gr = function(p) 2 * p))
#'
#' # What the difference costs, on Rosenbrock from the origin. Both runs take
#' # 21 iterations and land in the same place; only the bill differs.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#'
#' a <- minimize(bfgs(), f, c(0, 0))
#' b <- minimize(bfgs(), f, c(0, 0), gr = g)
#' unlist(a@counts)          # 126 objective evaluations, no gradients
#' unlist(b@counts)          # 30 objective evaluations and 24 gradients
#' max(abs(a@par - b@par))   # the two answers agree to 1.5e-08
#'
#' # The run that differenced says so.
#' a@message
#'
#' # A class with no method is refused by name.
#' try(as_objective(1:3))
#'
#' @seealso [minimize()] for the run itself, [newton()] for the one method
#'   that reads a Hessian, [optimizer_provides()] for what a method offers a
#'   stopping rule in return.
#' @export
as_objective <- S7::new_generic("as_objective", "fn",
                                function(fn, gr = NULL, he = NULL, ...)
                                  S7::S7_dispatch())


#' @title An Ordinary R Function as an Objective
#' @name as_objective.function
#'
#' @description
#' The shipped method, and the case almost every caller is in: `fn(par)`
#' returns the number to be minimized, and `gr(par)` and `he(par)` return
#' its gradient and Hessian when the caller has them. Every evaluation is a
#' callback into R, which is the cost of accepting an arbitrary closure and
#' the reason a run reports its evaluation counts.
#'
#' @param fn A function of the parameter vector returning a single number.
#' @param gr A function of the parameter vector returning the gradient, or
#'   `NULL` for a central difference of `fn`.
#' @param he A function of the parameter vector returning the Hessian, or
#'   `NULL`. Read by [newton()] alone.
#' @param ... Unused.
#'
#' @return The six-component list described under [as_objective()], with
#'   `kind = "r"`.
#'
#' @section Validation:
#' `gr` and `he` are each checked to be a function or `NULL`, and anything
#' else raises immediately, before the run starts. `fn` itself is not
#' checked here: dispatch has already established that it is a function, and
#' whether it returns a single number is settled at the first evaluation,
#' where the offending value can be reported.
#'
#' @examples
#' obj <- as_objective(function(p) sum(p^2), gr = function(p) 2 * p)
#' obj$kind
#' c(obj$has_gradient, obj$has_hessian)
#'
#' # A gradient that is not a function is refused before any evaluation.
#' try(as_objective(function(p) sum(p^2), gr = 1))
#'
#' @keywords internal
S7::method(as_objective, S7::class_function) <-
  function(fn, gr = NULL, he = NULL, ...) {
    if (!is.null(gr) && !is.function(gr)) {
      stop("'gr' must be a function or NULL.", call. = FALSE)
    }
    if (!is.null(he) && !is.function(he)) {
      stop("'he' must be a function or NULL.", call. = FALSE)
    }
    list(kind = "r", fn = fn, gr = gr, he = he,
         has_gradient = !is.null(gr), has_hessian = !is.null(he))
  }
