#' @include optimizer_class.R
#' @include generics.R
#' @include line_search.R
#' @include second_order.R
NULL


#' @title S7 Class for Chebyshev's Method
#'
#' @description
#' This class is a [Newton] optimizer whose direction carries Chebyshev's
#' third-order correction. [chebyshev()] builds it.
#'
#' @details
#' Beyond the properties of [Newton], a `Chebyshev` carries `ratio`, the
#' largest relative length of the correction at which the corrected
#' direction is used; above it the Newton step is used.
#'
#' @param ratio The largest \eqn{\lVert d - d_N\rVert / \lVert d_N\rVert} at
#'   which the corrected direction \eqn{d} is used.
#' @inheritParams Newton-class
#'
#' @return An S7 object of class `Chebyshev` inheriting from [Newton].
#'
#' @seealso [chebyshev()] for the constructor.
#' @name Chebyshev-class
#' @aliases Chebyshev
#' @keywords internal
Chebyshev <- S7::new_class("Chebyshev", parent = Newton,
  properties = list(ratio = S7::class_numeric))


#' @title Chebyshev's Third-Order Method
#'
#' @description
#' Newton's method with Chebyshev's correction: where the Hessian \eqn{H} is
#' positive definite the Newton step \eqn{d_N = -H^{-1}g} is corrected by the
#' third derivative of the objective,
#' \deqn{d = d_N - \tfrac12 H^{-1}\, T[d_N, d_N],}
#' where \eqn{T[u, v]_i = \sum_{jk} \partial^3 f/\partial x_i\partial x_j
#' \partial x_k\, u_j v_k}. Near a minimum the iteration converges cubically
#' where Newton's converges quadratically. Everything else (the repair of an
#' indefinite Hessian, the bound on the length of a step, the line search,
#' the stopping rule) is as in [newton()].
#'
#' @param ratio The largest relative change of the Newton step that the
#'   correction may make, \eqn{\lVert d - d_N\rVert \le
#'   \texttt{ratio}\,\lVert d_N\rVert}. A larger correction, or one that is
#'   not a descent direction, indicates that the cubic model does not
#'   describe the objective over the length of the step, and the Newton step
#'   is used instead. Defaults to `0.5`.
#' @param criterion,hessian_mod,floor,step,line_search,maxit,max_eval,verbose,refresh,keep_trace,max_length,typical
#'   As in [newton()].
#'
#' @details
#' # The third derivative
#'
#' The caller supplies to [minimize()] a function `t3(x, d)` that returns
#' the vector \eqn{T[d, d]}, the third derivative contracted twice with
#' \eqn{d}, and [minimize()] passes it to the method. The tensor itself is
#' never needed: for a sum over observations
#' \eqn{f = \sum_i \phi_i(x_i^\top\beta)} the contraction is
#' \eqn{T[d, d] = X^\top(\phi'''\,(Xd)^2)}, which costs as much as a
#' gradient. Without `t3` the contraction is one second difference of the
#' gradient along \eqn{d}, at two gradient evaluations in each iteration at
#' which the correction is attempted.
#'
#' # Where the correction is applied
#'
#' The correction is attempted only at an iterate where \eqn{H} is positive
#' definite and its Newton step is a descent direction; where the Hessian has
#' to be repaired, the step is the repaired Newton step, as in [newton()].
#' The correction costs one further solve with the same matrix. A correction
#' that is not used is recorded in the trace as `cubic correction refused`.
#'
#' Near the mode the correction is usually accepted, and the method needs
#' fewer iterations than [newton()]. Far from the mode the correction is
#' usually rejected, and the method then behaves as [newton()]. When the
#' starting point is already within one or two Newton steps of the mode, as
#' for an inner fit restarted from the previous mode, a faster local rate has
#' little to shorten, and the method gives no gain over [newton()].
#'
#' # Halley's method
#'
#' Halley's method solves with \eqn{H + \tfrac12 T[d_N]} in place of
#' \eqn{H}. This needs the contraction as a matrix and a second
#' factorization, of a matrix that need not be positive definite, and the
#' package does not provide it.
#'
#' # Box constraints
#'
#' With `lower` or `upper` the method runs on the unconstrained scale, as
#' every method that accepts bounds does. A `t3` written for the original
#' parameters is then the wrong derivative, so `t3` together with bounds is
#' rejected; without `t3` the second difference of the gradient is used.
#'
#' @return An S7 object of class [Chebyshev], inheriting from [newton()]'s
#'   class, to be handed to [minimize()].
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                           200 * (p[2] - p[1]^2))
#' rosen_he <- function(p) matrix(
#'   c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
#'     -400 * p[1], 200), 2, 2)
#' # The example computes T[d, d] of Rosenbrock: the only third derivatives are
#' # f_111 = 2400 p1 and f_112 = -400
#' rosen_t3 <- function(p, d) c(2400 * p[1] * d[1]^2 - 800 * d[1] * d[2],
#'                              -400 * d[1]^2)
#'
#' a <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
#' b <- minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he,
#'               t3 = rosen_t3)
#' c(newton = a@iterations, chebyshev = b@iterations)
#'
#' # Without t3 the contraction is a second difference of the gradient.
#' minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
#'
#' @seealso [newton()], whose settings this method shares.
#' @references
#' Gundersen, G. and Steihaug, T. (2010). On large-scale unconstrained
#' optimization problems and higher order methods. *Optimization Methods
#' and Software* **25**, 337--358.
#'
#' Gutierrez, J. M. and Hernandez, M. A. (1997). A family of
#' Chebyshev-Halley type methods in Banach spaces. *Bulletin of the
#' Australian Mathematical Society* **55**, 113--130.
#'
#' @export
chebyshev <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                      ratio = 0.5,
                      hessian_mod = c("eigen", "ridge"), floor = 1e-8,
                      step = 1, line_search = armijo(),
                      maxit = 200, max_eval = Inf,
                      verbose = FALSE, refresh = 10, keep_trace = FALSE,
                      max_length = Inf, typical = NULL) {
  base <- newton(criterion = criterion, hessian_mod = hessian_mod,
                 floor = floor, step = step, line_search = line_search,
                 maxit = maxit, max_eval = max_eval, verbose = verbose,
                 refresh = refresh, keep_trace = keep_trace,
                 max_length = max_length, typical = typical)
  if (length(ratio) != 1L || !is.numeric(ratio) || is.na(ratio) ||
      ratio <= 0) {
    stop("'ratio' must be a single positive number (Inf for no bound).",
         call. = FALSE)
  }
  Chebyshev(name = "Chebyshev", criterion = base@criterion,
            maxit = base@maxit, max_eval = base@max_eval,
            verbose = base@verbose, refresh = base@refresh,
            keep_trace = base@keep_trace, step = base@step,
            line_search = base@line_search, hessian_mod = base@hessian_mod,
            floor = base@floor, max_length = base@max_length,
            typical = base@typical, ratio = ratio)
}


#' @title Minimize by Chebyshev's Method
#' @name minimize.Chebyshev
#' @param optimizer A `Chebyshev` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @param t3 `NULL`, or a function of the parameters and a direction
#'   returning the third derivative of `fn` contracted twice with the
#'   direction, a vector as long as the parameters.
#' @return An [optimizer_result()].
#' @keywords internal
S7::method(minimize, Chebyshev) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, t3 = NULL,
           lower = -Inf, upper = Inf, ...) {
    if (!is.null(t3) && !is.function(t3)) {
      stop("'t3' must be a function of the parameters and a direction, or NULL.",
           call. = FALSE)
    }
    if (!is.null(t3) && length(check_bounds(lower, upper, par))) {
      stop(paste0("'t3' is the third derivative in the parameters, and with ",
                  "'lower' or 'upper' the method runs on the unconstrained ",
                  "scale, where it does not apply. Drop 't3' (a difference ",
                  "of the gradient is used) or the bounds."), call. = FALSE)
    }
    typ <- optimizer@typical
    if (!is.null(typ)) {
      if (length(typ) == 1L) typ <- rep(typ, length(par))
      if (length(typ) != length(par)) {
        stop(sprintf(paste0("'typical' has %d sizes and 'par' has %d ",
                            "values; give one size, or one per parameter."),
                     length(typ), length(par)), call. = FALSE)
      }
    }
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "chebyshev", hessian_mod = optimizer@hessian_mod,
                     floor = optimizer@floor,
                     max_length = optimizer@max_length,
                     typical = if (is.null(typ)) numeric(0) else typ,
                     ratio = optimizer@ratio, t3 = t3))
  }
