#' @include optimizer_class.R
NULL

#' @title S7 Class for the Proximal Gradient Method
#'
#' @description
#' This class is an optimizer that holds the two descriptions of the
#' non-smooth part of an objective, its proximal operator and its value,
#' together with the settings of the accelerated iteration. It is built by
#' [prox_grad()]. It is the only shipped class whose [optimizer_bounded()] is `FALSE`: its constraint
#' is passed inside `prox`, so box bounds are rejected.
#'
#' @details
#' Beyond the seven properties of every optimizer, a `ProxGrad` carries six
#' of its own. `prox` and `g` describe the same non-smooth term in two ways
#' and are both required. `accelerate`, `step`, `shrink` and `restart`
#' govern the iteration.
#'
#' @param prox The proximal operator of the non-smooth part, `prox(v, step)`.
#' @param g The value of the non-smooth part, `g(par)`.
#' @param accelerate Logical; whether the momentum extrapolation is applied.
#' @param step The initial step length offered to the backtracking search.
#' @param shrink The factor by which a rejected step is multiplied.
#' @param restart Logical; whether an increase in the objective resets the
#'   momentum.
#'
#' @return An S7 object of class `ProxGrad` inheriting from [optimizer()],
#'   with the six properties above beside the seven shared ones.
#'
#' @seealso [prox_grad()] for the constructor, [minimize.ProxGrad()] for the
#'   run.
#' @name ProxGrad-class
#' @aliases ProxGrad
#' @keywords internal
ProxGrad <- S7::new_class("ProxGrad", parent = optimizer,
  properties = list(
    prox       = S7::class_function,
    g          = S7::class_function,
    accelerate = S7::class_logical,
    step       = S7::class_numeric,
    shrink     = S7::class_numeric,
    restart    = S7::class_logical
  ))


#' Proximal Gradient Method
#'
#' @description
#' Minimizes \eqn{f(x) + g(x)} where \eqn{f} is smooth and \eqn{g} is
#' handled entirely through its proximal operator, so that a
#' non-differentiable term is minimized without being differenced. With
#' `accelerate = TRUE` the method is an accelerated iteration of the kind of
#' Beck and Teboulle, with momentum weight \eqn{(k-1)/(k+2)}, whose objective
#' gap falls as \eqn{O(1/k^2)} against the \eqn{O(1/k)} of the plain
#' iteration.
#'
#' @details
#' Each iteration takes a gradient step on the smooth part and applies the
#' proximal operator to the result,
#' \deqn{x_{k+1} = \mathrm{prox}_{t g}\big(y_k - t\nabla f(y_k)\big),}
#' with \eqn{y_k = x_k} for the plain method and
#' \eqn{y_k = x_k + \frac{k-1}{k+2}(x_k - x_{k-1})} for the accelerated
#' one. The step length is found by backtracking: \eqn{t} is multiplied by
#' `shrink` (one half by default) until the quadratic model built at
#' \eqn{y_k} dominates \eqn{f} at the new point, which is the condition used
#' by the convergence proof and needs no knowledge of the Lipschitz constant.
#' The step length is not increased again during the run.
#'
#' # What the stopping rule reads
#'
#' The gradient of \eqn{f + g} does not vanish at the solution, since \eqn{g}
#' is not differentiable there, so this method reports the proximal
#' gradient mapping
#' \deqn{G_t(x) = \frac{x - \mathrm{prox}_{tg}(x - t\nabla f(x))}{t}}
#' as its gradient, read at the iterate \eqn{x}, the point that is reported.
#' It vanishes exactly at a stationary point of \eqn{f + g} and reduces to
#' \eqn{\nabla f} when \eqn{g} is absent, so [crit_grad()] keeps its meaning
#' and its default tolerance. On the lasso of the examples the mapping agrees
#' with the KKT conditions, which the last example checks.
#'
#' Without momentum, the next iteration starts at the iterate, so the
#' gradient read there serves as the gradient of the next step, and the
#' plain method evaluates the gradient once per iteration. The accelerated
#' method takes its step at the extrapolated point, which is not the point
#' that it reports, and evaluates the gradient about twice per iteration.
#'
#' # Acceleration
#'
#' On an ill-conditioned problem the accelerated iteration needs far fewer
#' iterations than the plain one, as the \eqn{O(1/k^2)} and \eqn{O(1/k)}
#' rates indicate. On a well-conditioned problem with a tight tolerance the
#' plain iteration can need fewer: near the solution the objective changes
#' at the level of its rounding error, and each restart triggered there
#' discards the momentum built since the previous one. At the default
#' `crit_grad(1e-6)` this does not occur on the lasso of the examples.
#'
#' # Restarting
#'
#' Momentum makes the objective non-monotone, and an increase far from the
#' solution indicates momentum built in the wrong direction. With
#' `restart = TRUE` an increase resets the extrapolation to the current
#' point. This is the adaptive restart of O'Donoghue and Candes, and it costs
#' one comparison per iteration. On a badly conditioned problem it reduces
#' the number of iterations substantially.
#'
#' An increase is measured against an allowance for the rounding error of
#' the objective, \eqn{8\,\varepsilon\,\max(1, \lvert f \rvert)}, where
#' \eqn{f} is the current total objective, and not against zero. The
#' objective is a sum, so its rounding error grows with the number of terms,
#' and a bare `>` test fires on that error once the iteration is close
#' enough to the solution, restarting the momentum at nearly every
#' iteration. A larger allowance would suppress real restarts on an
#' ill-conditioned problem. The error of a sum of \eqn{n} terms can reach
#' about \eqn{n/2} units in the last place, so for a large sum the allowance
#' is smaller than the possible error.
#'
#' @param prox The proximal operator of the non-smooth part, a function of
#'   the point and the step length, `prox(v, step)`, returning the
#'   minimizer of \eqn{\lVert b - v \rVert^2/(2\,\mathrm{step}) + g(b)}.
#'   [penalties7::penalty_prox()] supplies one for every penalty
#'   that has it.
#' @param g The value of the non-smooth part, a function of the point.
#'   It is required alongside `prox`, because the two describe the same
#'   term, and without `g` the reported objective would include only the
#'   smooth part.
#' @param accelerate Apply the momentum extrapolation? Defaults to
#'   `TRUE`.
#' @param step The initial step length offered to the backtracking search.
#'   Defaults to `1`.
#' @param shrink The factor by which a rejected step is reduced. Defaults
#'   to `0.5`.
#' @param restart Reset the momentum when the objective increases by more
#'   than the allowance for its rounding error? Defaults to `TRUE`, and is
#'   ignored when `accelerate` is `FALSE`. See the details about a tight
#'   tolerance.
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   [crit_grad()], which here reads the proximal gradient mapping.
#' @param maxit Maximum iterations, a finite number at least 1. Defaults to
#'   1000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @return An S7 object of class [ProxGrad], inheriting from [optimizer()],
#'   to be handed to [minimize()]. Box bounds are rejected: pass the
#'   constraint through `prox` instead, composing the projection onto the box
#'   with the operator.
#'
#' @references
#' Beck, A. and Teboulle, M. (2009). A fast iterative shrinkage-thresholding
#' algorithm for linear inverse problems. *SIAM Journal on Imaging
#' Sciences*, 2(1), 183--202.
#'
#' O'Donoghue, B. and Candes, E. (2015). Adaptive restart for accelerated
#' gradient schemes. *Foundations of Computational Mathematics*, 15(3),
#' 715--732.
#'
#' @seealso [bundle()] for a non-smooth method that needs no
#'   proximal operator, [gd()] for the smooth case.
#'
#' @examples
#' # A lasso-penalized least squares problem is solved through the operator.
#' # Three of the eight coefficients are non-zero in the truth.
#' set.seed(1)
#' X <- matrix(rnorm(200 * 8), 200, 8)
#' b0 <- c(2, -1.5, 0, 0, 0.8, 0, 0, 0)
#' y <- as.numeric(X %*% b0 + rnorm(200))
#' lambda <- 0.4
#'
#' f  <- function(b) sum((y - X %*% b)^2) / (2 * nrow(X))
#' gf <- function(b) as.numeric(-crossprod(X, y - X %*% b) / nrow(X))
#'
#' fit <- minimize(
#'   prox_grad(prox = function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0),
#'             g = function(b) lambda * sum(abs(b))),
#'   f, gr = gf, par = rep(0, 8))
#' round(fit@par, 3)
#' sum(fit@par != 0)          # the operator sets coefficients exactly to zero
#'
#' # The KKT conditions confirm the answer, and share no arithmetic with the
#' # iteration: stationarity where a coefficient is not zero, and a gradient
#' # inside the interval opened by the kink where it is.
#' gr_at <- gf(fit@par)
#' active <- fit@par != 0
#' max(abs(gr_at[active] + lambda * sign(fit@par[active])))
#' max(abs(gr_at[!active])) < lambda
#'
#' # Box bounds are rejected, and the message names the argument that takes
#' # the constraint.
#' try(minimize(prox_grad(prox = function(v, t) v, g = function(b) 0),
#'              f, rep(0, 8), gr = gf, lower = 0))
#'
#' @export
prox_grad <- function(prox, g,
                      accelerate = TRUE, step = 1, shrink = 0.5,
                      restart = TRUE,
                      criterion = crit_grad(),
                      maxit = 1000, max_eval = Inf,
                      verbose = FALSE, refresh = 20, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  if (!is.function(prox)) {
    stop("'prox' must be a function of the point and the step length.",
         call. = FALSE)
  }
  if (missing(g) || !is.function(g)) {
    stop(paste0("'g' must be a function giving the value of the non-smooth ",
                "part.\n  Without it the reported objective would include ",
                "only the smooth part."), call. = FALSE)
  }
  for (nm in c("accelerate", "restart")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", nm), call. = FALSE)
    }
  }
  if (!is.numeric(shrink) || length(shrink) != 1L || is.na(shrink) ||
      shrink <= 0 || shrink >= 1) {
    stop("'shrink' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  ProxGrad(
    name = if (accelerate) "accelerated proximal gradient" else "proximal gradient",
    criterion = criterion, maxit = maxit, max_eval = max_eval,
    verbose = verbose, refresh = refresh, keep_trace = keep_trace,
    prox = prox, g = g, accelerate = accelerate,
    step = step, shrink = shrink, restart = restart
  )
}


#' @title The Proximal Gradient Method Does Not Take Box Bounds
#' @name optimizer_bounded.ProxGrad
#'
#' @description
#' Returns `FALSE`; it is the only shipped method that does. A box
#' constraint is itself a non-smooth term, expressed by the projection onto
#' the box, and this method already receives such a term through `prox`.
#' Bounds beside it would be a second, conflicting way of imposing a
#' constraint, and the two operators would have to be composed.
#'
#' @details
#' [check_optimizer()] reads this declaration before its bounds check, so a
#' method that returns `FALSE` is not failed for rejecting bounds. For a
#' penalty that is convex and separable coordinate by coordinate, such as
#' the lasso or the ridge, the composition
#' `prox(v, t) = pmin(pmax(prox_penalty(v, t), lower), upper)` is the
#' proximal operator of the penalty plus the box. For a nonconvex or a
#' non-separable penalty the composition is not the operator of the sum.
#'
#' @param optimizer A `ProxGrad` object.
#'
#' @return `FALSE`.
#'
#' @examples
#' pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
#' optimizer_bounded(pg)
#' optimizer_bounded(bfgs())
#'
#' @keywords internal
S7::method(optimizer_bounded, ProxGrad) <- function(optimizer) FALSE


#' @title Minimize by the Proximal Gradient Method
#' @name minimize.ProxGrad
#'
#' @description
#' Runs [prox_grad()] on the objective: a backtracked gradient step on the
#' smooth part, the proximal operator applied to the result, and the momentum
#' extrapolation with its restart. Box bounds are rejected here, with a
#' message that names `prox` as the place for the constraint.
#'
#' @param optimizer A `ProxGrad` object.
#' @param fn,par,gr,hess,lower,upper,... As in [minimize()]. A finite `lower`
#'   or `upper` raises an error; `hess` is accepted and ignored.
#'
#' @return An [optimizer_result()] whose `value` is the total objective
#'   \eqn{f(x) + g(x)}, which is why `g` is required alongside `prox`, and
#'   whose `gradient` is the proximal gradient mapping at `par` rather than
#'   \eqn{\nabla f}.
#'
#' @keywords internal
S7::method(minimize, ProxGrad) <-
  function(optimizer, fn, par, gr = NULL, hess = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, hess)
    if (any(is.finite(c(lower, upper)))) {
      stop(paste0("prox_grad() takes its constraint through the proximal ",
                  "operator, not through box bounds:\n  compose the ",
                  "projection onto the box into 'prox'."), call. = FALSE)
    }
    t0 <- proc.time()[["elapsed"]]
    out <- prox_grad_run(optimizer, spec, as.numeric(par))
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0)
  }


#' Run the Proximal Gradient Loop
#'
#' @description
#' The iteration behind [prox_grad()]: a backtracked gradient step on the
#' smooth part, the proximal operator applied to its result, and the momentum
#' extrapolation with its restart. It counts its own evaluations, including
#' those of a gradient obtained by differences.
#'
#' @details
#' The loop is written in R, as that of [brent()]. Every iteration calls
#' three R functions supplied by the caller: the objective, its gradient and
#' the proximal operator.
#'
#' The proximal gradient mapping is read at the iterate, the point that
#' is reported, and not at the extrapolated point. The trace records
#' `step shortened` at an iteration whose step was backtracked and
#' `momentum restart` at one whose momentum was reset.
#'
#' @param optimizer A `ProxGrad` object.
#' @param spec The objective handle from [as_objective()].
#' @param par The starting point, a numeric vector.
#'
#' @return A list in the shape that [build_result()] consumes: the point, the
#'   total objective there, the mapping, the counts, the iteration count, the
#'   verdict and the trace.
#'
#' @keywords internal
prox_grad_run <- function(optimizer, spec, par) {
  n_value <- 0L
  n_grad <- 0L
  f <- function(x) {
    n_value <<- n_value + 1L
    spec$fn(x)
  }
  gradf <- if (spec$has_gradient) {
    function(x) {
      n_grad <<- n_grad + 1L
      as.numeric(spec$gr(x))
    }
  } else {
    # differenced through f, so that the evaluations it costs are counted and
    # are charged to max_eval like any other
    function(x) {
      n_grad <<- n_grad + 1L
      fd_gradient(f, x)
    }
  }
  gval <- optimizer@g
  prox <- optimizer@prox
  total <- function(x, fx = NULL) (if (is.null(fx)) f(x) else fx) + gval(x)

  x <- par
  y <- x
  t <- optimizer@step
  fx_total <- total(x)
  f_prev <- NA_real_
  k <- 1L
  momentum <- 0
  trace <- if (optimizer@keep_trace) list() else NULL
  converged <- FALSE
  stopped_by <- "maxit"
  it <- 0L
  gmap <- rep(NA_real_, length(x))
  # the smooth value and gradient at the current iterate, read when the
  # mapping is computed there; the next iteration starts from the same point
  # whenever the extrapolation is zero (no momentum, the first iteration
  # after a restart, or accelerate = FALSE), and reuses them
  x_cache <- NULL
  f_cache <- NA_real_
  g_cache <- NULL

  for (it in seq_len(as.integer(optimizer@maxit))) {
    # A restart retries WITHIN the iteration. Spending an iteration on it
    # would make the accelerated method look slower than the plain one for
    # no reason but the bookkeeping, and would leave a gap in the trace.
    # The k > 1 guard makes the retry provably finite: at k = 1 there is no
    # momentum to discard, and the backtracked step decreases the objective.
    guard <- "none"
    repeat {
      if (!is.null(x_cache) && identical(y, x_cache)) {
        fy <- f_cache
        gy <- g_cache
      } else {
        fy <- f(y)
        gy <- gradf(y)
      }

      # backtracking: shrink until the quadratic model at y dominates f at
      # the candidate, which is the inequality the convergence proof uses
      repeat {
        x_new <- as.numeric(prox(y - t * gy, t))
        d <- x_new - y
        f_cand <- f(x_new)
        if (!is.finite(f_cand) ||
            f_cand > fy + sum(gy * d) + sum(d^2) / (2 * t)) {
          t <- t * optimizer@shrink
          if (guard == "none") guard <- "step shortened"
          if (t < .Machine$double.eps) break
        } else {
          break
        }
      }

      gmap <- (y - x_new) / t
      f_new_total <- f_cand + gval(x_new)

      # "Did the objective increase?" has to allow for the objective's own
      # rounding, or near the solution the test fires on it and the momentum
      # is discarded over and over. The objective is a sum, whose error grows
      # with the number of terms, so the allowance is a few ulps of its own
      # size rather than zero; measured, anything from 1 to 32 gives the same
      # run and 256 starts costing iterations on an ill-conditioned problem.
      slack <- 8 * .Machine$double.eps * max(1, abs(fx_total))
      if (optimizer@accelerate && optimizer@restart && k > 1L &&
          f_new_total > fx_total + slack) {
        y <- x
        k <- 1L
        guard <- "momentum restart"
      } else {
        break
      }
    }

    x_old <- x
    f_prev <- fx_total
    x <- x_new
    fx_total <- f_new_total

    # The mapping is read at the ITERATE x, the point that is reported. The
    # quantity computed above is the mapping at y, which is the previous
    # iterate without momentum and the extrapolated point with it. The
    # gradient read here is kept for the next iteration, which starts at x
    # whenever the extrapolation is zero.
    gx <- gradf(x)
    gmap <- (x - as.numeric(prox(x - t * gx, t))) / t
    x_cache <- x
    f_cache <- f_cand
    g_cache <- gx

    state <- list(iter = it, f_new = fx_total, f_old = f_prev,
                  x_new = x, x_old = x_old, gradient = gmap,
                  stationarity = NULL)
    met <- crit_met(optimizer@criterion, state)

    if (optimizer@keep_trace) {
      trace[[length(trace) + 1L]] <- data.frame(
        iteration = it, value = fx_total, gnorm = max(abs(gmap)), step = t,
        safeguard = guard, stringsAsFactors = FALSE)
    }
    if (optimizer@verbose && optimizer@refresh > 0 &&
        it %% optimizer@refresh == 0) {
      message(sprintf("  iter %4d  value %.6g  |G| %.3g  step %.3g",
                      it, fx_total, max(abs(gmap)), t))
    }

    if (isTRUE(met)) {
      converged <- TRUE
      stopped_by <- "criterion"
      break
    }
    if (n_value >= budget_int(optimizer@max_eval)) {
      stopped_by <- "max_eval"
      break
    }

    momentum <- if (optimizer@accelerate) (k - 1) / (k + 2) else 0
    y <- x + momentum * (x - x_old)
    k <- k + 1L
  }

  list(
    par = x, value = fx_total, gradient = gmap,
    n_value = n_value, n_grad = n_grad, n_hess = 0L,
    iterations = it, converged = converged,
    stopped_by = stopped_by,
    message = if (converged) "" else
      "the iteration ended before the stopping rule fired",
    trace = if (is.null(trace)) NULL else do.call(rbind, trace)
  )
}


#' A Central-Difference Gradient
#'
#' @description
#' Computes the finite-difference gradient that the proximal loop uses when
#' the caller supplies none. The step is scaled to each coordinate.
#'
#' @param fn A function of the parameter vector.
#' @param x The point.
#'
#' @return A numeric vector.
#'
#' @keywords internal
fd_gradient <- function(fn, x) {
  h <- .Machine$double.eps^(1 / 3) * pmax(1, abs(x))
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- x[j] + h[j]
    xm[j] <- x[j] - h[j]
    (fn(xp) - fn(xm)) / (2 * h[j])
  }, numeric(1))
}
