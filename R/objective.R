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
#' and Hessian the caller has, into the handle on which every algorithm in
#' the package operates. The result is a named list that contains the pieces
#' and two flags recording which derivatives were supplied. `as_objective()`
#' is a generic that dispatches on `fn`, so a caller with another kind of
#' objective registers one method and every algorithm accepts that objective.
#'
#' @details
#' # Two generics, each dispatching once
#'
#' [minimize()] dispatches on the optimizer, so an algorithm is written once
#' for every kind of objective. `as_objective()` dispatches on the objective,
#' so each kind of objective is converted once for every algorithm.
#' Dispatching [minimize()] on both would need one method for each pair of an
#' algorithm and a kind of objective. The package ships one method, for an
#' ordinary R function.
#'
#' # What happens when no gradient is supplied
#'
#' The compiled loop differences the objective, once per coordinate,
#' \deqn{\hat{g}_j = \frac{f(x + h e_j) - f(x - h e_j)}{2h}, \qquad
#'       h = \varepsilon^{1/3} \max(1, \lvert x_j \rvert),}
#' with \eqn{\varepsilon} the machine epsilon, so \eqn{h} is about
#' `6.06e-06` for a coordinate of size one. One gradient then costs
#' \eqn{2p} evaluations of the objective, and its relative error is of the
#' order of \eqn{\varepsilon^{2/3}}.
#'
#' A run that differences its gradient records this in the `message` field
#' of the result.
#'
#' # A Hessian by differencing the gradient
#'
#' [newton()] and [chebyshev()] read a Hessian, and when none is supplied it
#' is obtained by differencing the gradient once. With an analytic gradient
#' this is a single differentiation and costs \eqn{2p} gradient evaluations.
#' Without a gradient, the gradient is itself a difference and the Hessian is
#' a difference of differences, which costs about \eqn{4p^2} evaluations of
#' the objective and is less accurate. The Hessian is symmetrized before it
#' is returned, since its two triangles differ by the error of the
#' differences.
#'
#' # Points where the objective is not finite
#'
#' If the objective returns `NA`, `NaN` or an infinity at a point other than
#' the starting value, the algorithms treat that point as unusable, reject
#' the trial step and continue, so a line search can try points outside the
#' region where the objective is defined. A value that is not finite at the
#' starting value is an error.
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
#'   `p x p` symmetric matrix. `NULL` is the default. [newton()] and
#'   [chebyshev()] read it; every other method accepts it and does not use
#'   it, so calling code can pass whatever it has without branching on the
#'   method.
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
#' # The flags record what was supplied, and the compiled loop branches on them.
#' str(as_objective(function(p) sum(p^2)))
#' str(as_objective(function(p) sum(p^2), gr = function(p) 2 * p))
#'
#' # The cost of the difference, on Rosenbrock from the origin. Both runs take
#' # 21 iterations and reach the same point; the evaluation counts differ.
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
#' # The message of the run that differenced records it.
#' a@message
#'
#' # A class with no method is rejected, and the error names the class.
#' try(as_objective(1:3))
#'
#' @seealso [minimize()] for the run itself, [newton()] and [chebyshev()]
#'   for the methods that read a Hessian, [optimizer_provides()] for the
#'   quantities that a method supplies to a stopping rule.
#' @export
as_objective <- S7::new_generic("as_objective", "fn",
                                function(fn, gr = NULL, he = NULL, ...)
                                  S7::S7_dispatch())


#' @title An Ordinary R Function as an Objective
#' @name as_objective.function
#'
#' @description
#' The shipped method: `fn(par)` returns the number to be minimized, and
#' `gr(par)` and `he(par)` return its gradient and Hessian when the caller
#' has them. Every evaluation is a callback into R. The compiled objective
#' counts the evaluations, and the result reports them in `counts`.
#'
#' @param fn A function of the parameter vector returning a single number.
#' @param gr A function of the parameter vector returning the gradient, or
#'   `NULL` for a central difference of `fn`.
#' @param he A function of the parameter vector returning the Hessian, or
#'   `NULL`. Read by [newton()] and [chebyshev()].
#' @param ... Unused.
#'
#' @return The six-component list described under [as_objective()], with
#'   `kind = "r"`.
#'
#' @section Validation:
#' `gr` and `he` are each checked to be a function or `NULL`, and any other
#' value raises an error before the run starts. `fn` is not checked here:
#' dispatch has already established that it is a function, and the compiled
#' objective checks at every evaluation that it returns a single number.
#'
#' @examples
#' obj <- as_objective(function(p) sum(p^2), gr = function(p) 2 * p)
#' obj$kind
#' c(obj$has_gradient, obj$has_hessian)
#'
#' # A gradient that is not a function is rejected before any evaluation.
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
