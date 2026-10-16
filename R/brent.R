#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for Brent's Method
#'
#' @description
#' A `Brent` object is an optimizer for a function of one variable. It holds
#' the length of the first step of the bracketing phase and the number of times
#' that phase may expand, and it is built by [brent()].
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Brent` carries two of
#' its own: `step`, the length of the first trial step, and `max_expand`, the
#' number of golden expansions that the bracketing phase may take.
#'
#' @param step The length of the first trial step.
#' @param max_expand The number of expansions allowed while bracketing.
#'
#' @return An S7 object of class `Brent` inheriting from [optimizer()], with
#'   the two properties above beside the seven shared ones.
#'
#' @seealso [brent()] for the constructor.
#' @name Brent-class
#' @aliases Brent
#' @keywords internal
Brent <- S7::new_class("Brent", parent = optimizer,
  properties = list(
    step       = S7::class_numeric,
    max_expand = S7::class_numeric
  ))


#' @title Brent's Method in One Dimension
#'
#' @description
#' Minimizes a function of one variable without derivatives: a bracket of the
#' minimum is found by golden expansion from the starting value, and the
#' bracket is then reduced by Brent's combination of parabolic interpolation
#' and golden-section steps.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_stationary(1e-5)`, which reads the half-width of the bracket, so
#'   the tolerance is on the scale of the parameter, or on the unconstrained
#'   scale when `lower` or `upper` is given. A rule that needs a gradient is
#'   rejected when the run starts.
#' @param step The length of the first trial step of the bracketing phase, a
#'   single positive number, default `1`.
#' @param max_expand How many golden expansions the bracketing phase may
#'   take, a single number of at least 1, default `50`. If the phase has not
#'   bracketed a minimum by then, the run stops at the lowest point evaluated
#'   and reports `converged = FALSE`.
#' @param maxit Maximum iterations of the reduction phase. Defaults to 100.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 1.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' The bracketing phase evaluates the starting value \eqn{x_0} and
#' \eqn{x_0 + h}, turns round if the second is the higher, and steps on in
#' the descending direction, where each step is \eqn{1.618} times the previous
#' one, until the value no longer falls. The three points then satisfy
#' \eqn{f(b) \le f(a)} and \eqn{f(b) \le f(c)} with \eqn{b} between \eqn{a}
#' and \eqn{c}, so a minimum lies in \eqn{[a, c]}. If the value is still
#' falling when `max_expand` expansions have been taken, or when the
#' evaluation budget runs out, the phase has not bracketed a minimum, and the
#' run stops at the lowest point evaluated, reports `converged = FALSE` and
#' records the reason in `message`.
#'
#' The reduction phase is the algorithm of Brent (1973, chapter 5). At each
#' iteration the parabola through the three best points is tried first; its
#' vertex is accepted when it lies inside the bracket and the step is less
#' than half the step before last, and a golden-section step into the larger
#' part of the bracket is taken otherwise. Parabolic steps converge
#' superlinearly on a smooth minimum, and the golden steps bound the number
#' of evaluations when the function is not well modeled by a parabola. The
#' stationarity measure is the half-width of the bracket.
#'
#' The rule receives `iter`, `f_new`, `f_old`, `x_new`, `x_old` and
#' `stationarity`; `f_old` and `x_old` are those of the previous iteration
#' and are `NULL` at the first. Besides the rule, the reduction stops on
#' Brent's own test, when the bracket has shrunk to about
#' \eqn{\sqrt{\varepsilon}\,\lvert x \rvert + 10^{-12}} around the best point
#' \eqn{x}, because at that width a further point cannot be distinguished from
#' \eqn{x}.
#' That stop reports `converged = TRUE` with `criterion_met` equal to
#' `bracket at machine resolution`, so a rule with a smaller tolerance than
#' this resolution, or [crit_never()], ends there.
#'
#' The method is for one parameter only; a starting value of any other
#' length is rejected. Bounds are removed by the same reparametrization as in
#' the other methods (see [bounded_transform()]).
#'
#' @return An S7 object of class [Brent], inheriting from [optimizer()], to
#'   be handed to [minimize()].
#'
#' @references
#' Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*.
#' Prentice-Hall, Englewood Cliffs.
#'
#' @examples
#' brent()
#'
#' # The method minimizes a smooth function of one variable.
#' f <- function(x) (x - 2)^2 + exp(-x)
#' r <- minimize(brent(), f, par = 0)
#' c(par = r@par, evaluations = r@counts[["f"]])
#'
#' # nelder_mead() is run on the same function for comparison.
#' minimize(nelder_mead(), f, par = 0)@counts[["f"]]
#'
#' @seealso [nelder_mead()] for the derivative-free method in several
#'   dimensions, [crit_stationary()] for the default rule of this method.
#' @export
brent <- function(criterion = crit_stationary(1e-5), step = 1,
                  max_expand = 50, maxit = 100, max_eval = Inf,
                  verbose = FALSE, refresh = 1, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  if (length(step) != 1L || !is.numeric(step) || !is.finite(step) ||
      step <= 0) {
    stop("'step' must be a single positive number.", call. = FALSE)
  }
  if (length(max_expand) != 1L || !is.numeric(max_expand) ||
      is.na(max_expand) || max_expand < 1) {
    stop("'max_expand' must be a single number of at least 1.", call. = FALSE)
  }
  Brent(
    name = "brent", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    step = step, max_expand = max_expand
  )
}


#' @rdname optimizer_provides.NelderMead
#' @name optimizer_provides.Brent
#' @keywords internal
S7::method(optimizer_provides, Brent) <- function(optimizer) "stationarity"


#' @title Minimize by Brent's Method
#' @name minimize.Brent
#'
#' @description
#' Runs [brent()] on a function of one variable: bracket the minimum by
#' golden expansion, then reduce the bracket by parabolic and golden-section
#' steps until the stopping rule fires or the bracket reaches the resolution
#' of the iterate.
#'
#' @param optimizer A `Brent` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` and `he` are
#'   accepted and ignored, since the method uses no derivative. Bounds are taken
#'   and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is `NULL` and whose trace
#'   has the columns `iteration`, `value`, `stationarity` (the half-width of
#'   the bracket) and `safeguard`, which names the step taken at each
#'   iteration: `parabolic` or `golden`.
#'
#' @keywords internal
S7::method(minimize, Brent) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
    if (length(par) != 1L) {
      stop("brent() minimizes a function of one variable; 'par' has length ",
           length(par), ".", call. = FALSE)
    }
    bounds <- check_bounds(lower, upper, par)
    t0 <- proc.time()[["elapsed"]]
    out <- brent_run(optimizer, spec, as.numeric(par), bounds)
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0)
  }


#' The Two Phases of Brent's Method
#'
#' @description
#' This is the body of [minimize()] for a [Brent] optimizer: the bracketing
#' phase and the reduction phase, on the unconstrained scale where bounds are
#' given.
#'
#' @param optimizer A `Brent` object.
#' @param spec The objective handle from [as_objective()].
#' @param par The starting value, a single number.
#' @param bounds What [check_bounds()] returned.
#'
#' @return A list in the shape that [build_result()] consumes.
#'
#' @keywords internal
brent_run <- function(optimizer, spec, par, bounds) {
  n_value <- 0L
  budget <- budget_int(optimizer@max_eval)
  boxed <- length(bounds) > 0L
  to_par <- if (boxed) function(e) bounded_transform(bounds[[1L]], e)$h else
    identity
  f <- function(e) {
    n_value <<- n_value + 1L
    v <- spec$fn(to_par(e))
    if (is.finite(v)) v else Inf
  }
  out_of_budget <- function() n_value >= budget
  gold <- (1 + sqrt(5)) / 2
  cgold <- (3 - sqrt(5)) / 2
  trace <- if (optimizer@keep_trace) list() else NULL

  # BRACKETING: a < b < cc (or the reverse) with f(b) below both ends
  a <- if (boxed) bounded_forward(bounds[[1L]], par) else par
  fa <- f(a)
  b <- a + optimizer@step
  fb <- f(b)
  if (fb > fa) {
    tmp <- a; a <- b; b <- tmp
    tmp <- fa; fa <- fb; fb <- tmp
  }
  cc <- b + gold * (b - a)
  fc <- f(cc)
  expansions <- 1L
  while (fc < fb && expansions < optimizer@max_expand && !out_of_budget()) {
    a <- b; fa <- fb
    b <- cc; fb <- fc
    cc <- b + gold * (b - a)
    fc <- f(cc)
    expansions <- expansions + 1L
  }
  # the loop also ends on max_expand or on the budget while the value is still
  # falling; (a, cc) then contains no minimum, and the run stops at the lowest
  # point it has seen instead of reducing a bracket that is not one
  if (fc < fb) {
    out_budget <- out_of_budget()
    return(list(
      par = to_par(cc), value = fc, gradient = NULL,
      n_value = n_value, n_grad = 0L, n_hess = 0L,
      iterations = 0L, converged = FALSE,
      stopped_by = if (out_budget) "max_eval" else "failed",
      message = if (out_budget)
        "the evaluation budget ran out before a minimum was bracketed" else
        paste0("the bracketing phase found no minimum in max_expand = ",
               optimizer@max_expand, " expansions"),
      trace = NULL))
  }

  # REDUCTION, after Brent (1973): x is the best point, w the second best,
  # v the previous value of w; (lo, hi) is the bracket
  lo <- min(a, cc)
  hi <- max(a, cc)
  x <- w <- v <- b
  fx <- fw <- fv <- fb
  d <- e <- 0
  converged <- FALSE
  stopped_by <- "maxit"
  it <- 0L
  f_prev <- NA_real_
  x_prev <- NULL
  for (it in seq_len(as.integer(optimizer@maxit))) {
    mid <- (lo + hi) / 2
    tol1 <- sqrt(.Machine$double.eps) * abs(x) + 1e-12
    tol2 <- 2 * tol1
    half <- (hi - lo) / 2
    state <- list(iter = it, f_new = fx,
                  f_old = if (is.finite(f_prev)) f_prev else NULL,
                  x_new = x, x_old = x_prev, gradient = NULL,
                  stationarity = half)
    if (isTRUE(crit_met(optimizer@criterion, state))) {
      converged <- TRUE
      stopped_by <- "criterion"
      break
    }
    # Brent's own test: the bracket has shrunk to the resolution of x, where no
    # further step can be distinguished from x, whatever the rule asks for
    if (abs(x - mid) <= tol2 - half) {
      converged <- TRUE
      stopped_by <- "bracket"
      break
    }
    if (out_of_budget()) {
      stopped_by <- "max_eval"
      break
    }
    kind <- "golden"
    if (abs(e) > tol1) {
      # the parabola through x, w and v
      r <- (x - w) * (fx - fv)
      q <- (x - v) * (fx - fw)
      p <- (x - v) * q - (x - w) * r
      q <- 2 * (q - r)
      if (q > 0) p <- -p
      q <- abs(q)
      etemp <- e
      e <- d
      if (abs(p) < abs(0.5 * q * etemp) && p > q * (lo - x) && p < q * (hi - x)) {
        d <- p / q
        u <- x + d
        if (u - lo < tol2 || hi - u < tol2) d <- if (mid >= x) tol1 else -tol1
        kind <- "parabolic"
      }
    }
    if (kind == "golden") {
      e <- if (x >= mid) lo - x else hi - x
      d <- cgold * e
    }
    u <- if (abs(d) >= tol1) x + d else x + (if (d >= 0) tol1 else -tol1)
    fu <- f(u)
    f_prev <- fx
    x_prev <- x
    if (fu <= fx) {
      if (u >= x) lo <- x else hi <- x
      v <- w; fv <- fw
      w <- x; fw <- fx
      x <- u; fx <- fu
    } else {
      if (u < x) lo <- u else hi <- u
      if (fu <= fw || w == x) {
        v <- w; fv <- fw
        w <- u; fw <- fu
      } else if (fu <= fv || v == x || v == w) {
        v <- u; fv <- fu
      }
    }
    if (optimizer@keep_trace) {
      trace[[length(trace) + 1L]] <- data.frame(
        iteration = it, value = fx, stationarity = (hi - lo) / 2,
        safeguard = kind, stringsAsFactors = FALSE)
    }
    if (optimizer@verbose && optimizer@refresh > 0 &&
        it %% optimizer@refresh == 0) {
      message(sprintf("  iter %4d  value %.10g  half-width %.3g  %s",
                      it, fx, (hi - lo) / 2, kind))
    }
  }
  list(
    par = to_par(x), value = fx, gradient = NULL,
    n_value = n_value, n_grad = 0L, n_hess = 0L,
    iterations = it, converged = converged, stopped_by = stopped_by,
    message = if (converged) "" else
      "the iteration ended before the stopping rule fired",
    trace = if (is.null(trace)) NULL else do.call(rbind, trace)
  )
}
