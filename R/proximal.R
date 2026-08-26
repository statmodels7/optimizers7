#' @include optimizer_class.R
NULL

#' @title S7 Class for the Proximal Gradient Method
#'
#' @description
#' An optimizer holding the two descriptions of the non-smooth part of an
#' objective, its proximal operator and its value, together with the three
#' settings of the accelerated iteration. Built by [prox_grad()]. It is the
#' one shipped class whose [optimizer_bounded()] is `FALSE`: its constraint
#' travels inside `prox`, so box bounds are refused.
#'
#' @details
#' Beyond the seven properties every optimizer has, a `ProxGrad` carries six
#' of its own. `prox` and `g` describe the same non-smooth term from two
#' sides and are both required. `accelerate`, `step`, `shrink` and `restart`
#' govern the iteration.
#'
#' @param prox The proximal operator of the non-smooth part, `prox(v, step)`.
#' @param g The value of the non-smooth part, `g(par)`.
#' @param accelerate Logical; whether the momentum extrapolation is applied.
#' @param step The initial step length offered to the backtracking search.
#' @param shrink The factor a rejected step is multiplied by.
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
#' non-differentiable term is minimized without ever differencing it. With
#' `accelerate = TRUE` the method is the accelerated one of Beck and
#' Teboulle, whose objective gap falls as \eqn{O(1/k^2)} against the
#' \eqn{O(1/k)} of the plain iteration.
#'
#' @details
#' Each iteration takes a gradient step on the smooth part and applies the
#' proximal operator to the result,
#' \deqn{x_{k+1} = \mathrm{prox}_{t g}\big(y_k - t\nabla f(y_k)\big),}
#' with \eqn{y_k = x_k} for the plain method and
#' \eqn{y_k = x_k + \frac{k-1}{k+2}(x_k - x_{k-1})} for the accelerated
#' one. The step length is found by backtracking: \eqn{t} is halved until
#' the quadratic model built at \eqn{y_k} dominates \eqn{f} at the new
#' point, which is the condition the convergence proof uses and which
#' needs no knowledge of the Lipschitz constant.
#'
#' # What the stopping rule reads
#'
#' The gradient of \eqn{f + g} does not vanish at the solution, \eqn{g} being
#' non-differentiable there, so this method reports the **proximal gradient
#' mapping**
#' \deqn{G_t(x) = \frac{x - \mathrm{prox}_{tg}(x - t\nabla f(x))}{t}}
#' as its gradient, read at the iterate. It vanishes exactly at a stationary
#' point of \eqn{f + g} and reduces to \eqn{\nabla f} when \eqn{g} is absent,
#' so [crit_grad()] keeps its meaning and its default tolerance. Measured on
#' the lasso of the examples below, the mapping agrees with the KKT
#' violation: at the reported point the active coordinates satisfy
#' \eqn{\lvert \nabla f + \lambda \operatorname{sign}(\beta) \rvert \le}
#' `4.5e-07` and the inactive ones have \eqn{\lvert \nabla f \rvert \le
#' 0.129} against a \eqn{\lambda} of 0.4.
#'
#' Acceleration pays 2.1 gradient evaluations per iteration against the plain
#' method's 1.0, the extrapolated point at which it takes its step not being
#' the point it reports.
#'
#' # Acceleration, and where it earns its keep
#'
#' On an ill-conditioned problem the difference is the one the theory
#' predicts. Measured on a smooth quadratic in eight unknowns at
#' `crit_grad(1e-8)`, plain against accelerated: 64 iterations against 31 at
#' a condition number near 3, 1154 against 100 at 55, and 9321 against 268 at
#' 480.
#'
#' On a well-conditioned problem asked for a tight tolerance the plain
#' iteration is still the better one. The lasso below converges in 12
#' iterations at `crit_grad(1e-9)` with `accelerate = FALSE`, in 89 with
#' `restart = FALSE`, and in 782 with the defaults, all three reaching the
#' same support, the same objective to the last bit and the same
#' coefficients to `8e-09`. What costs the iterations is the restart, not the
#' momentum: near the solution the objective moves at the rounding level, and
#' each spurious reset discards the momentum built since the last one. At the default `crit_grad(1e-6)` none of this appears, both settings
#' taking about ten iterations.
#'
#' # Restarting
#'
#' Momentum makes the objective non-monotone, and an increase far from the
#' solution is a symptom of momentum built in the wrong direction. With
#' `restart = TRUE` an increase resets the extrapolation to the current
#' point. This is the adaptive restart of O'Donoghue and Candes and costs one
#' comparison per iteration. It is worth a great deal where the problem is
#' badly conditioned: on a smooth quadratic in eight unknowns at
#' `crit_grad(1e-8)`, with the restart against without, 150 iterations
#' against 858 at a condition number of 55, 573 against 5418 at 480, and
#' 1042 against 15104 at 2400.
#'
#' An increase is measured against the objective's own rounding, not against
#' zero. The objective is a sum, so its error grows with the number
#' of terms, and a test reading a bare `>` fires on that error once the
#' iteration is near enough to the solution: the lasso below took 21646
#' iterations at `crit_grad(1e-9)` before the allowance and takes 782 after
#' it, at the same answer. The allowance is eight units in the last place of
#' the current objective, measured to give the same run anywhere between one
#' and thirty-two while 256 begins costing iterations on the ill-conditioned
#' quadratic, where it suppresses restarts that are real. What it costs
#' there is about a tenth: 149, 571 and 934 iterations at the three
#' condition numbers before, against the 150, 573 and 1042 above.
#'
#' @param prox The proximal operator of the non-smooth part, a function of
#'   the point and the step length, `prox(v, step)`, returning the
#'   minimizer of \eqn{\lVert b - v \rVert^2/(2\,\mathrm{step}) + g(b)}.
#'   [penalties7::penalty_prox()] supplies one for every penalty
#'   that has it.
#' @param g The value of the non-smooth part, a function of the point.
#'   Required alongside `prox`: the two describe the same term, and
#'   without `g` the reported objective would be the smooth part
#'   alone.
#' @param accelerate Apply the momentum extrapolation? Defaults to
#'   `TRUE`.
#' @param step The initial step length offered to the backtracking search.
#'   Defaults to `1`.
#' @param shrink The factor by which a rejected step is reduced. Defaults
#'   to `0.5`.
#' @param restart Reset the momentum when the objective increases by more
#'   than the objective's own rounding? Defaults to `TRUE`, and is ignored
#'   when `accelerate` is `FALSE`. See the measured cost at a tight
#'   tolerance above.
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
#'   to be handed to [minimize()]. Box bounds are refused: pass the
#'   constraint through `prox` instead, composing the projection onto the box
#'   into it.
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
#' # A lasso-penalized least squares problem, solved through the operator.
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
#' # iteration: stationary where a coefficient survives, and inside the
#' # interval the kink opens where one does not.
#' gr_at <- gf(fit@par)
#' active <- fit@par != 0
#' max(abs(gr_at[active] + lambda * sign(fit@par[active])))
#' max(abs(gr_at[!active])) < lambda
#'
#' # Box bounds are refused, and the message says where the constraint goes.
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
                "part.\n  Without it the reported objective would be the ",
                "smooth part alone."), call. = FALSE)
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
#' Returns `FALSE`, the only shipped method that does. A box constraint is
#' itself a non-smooth term, expressed by the projection onto the box, and
#' this method already has a slot for such a term: `prox`. Offering bounds
#' beside it would be a second and conflicting route to the same thing, and
#' the two operators would have to be composed by somebody.
#'
#' @details
#' [check_optimizer()] consults this before running its bounds check, so a
#' method that answers `FALSE` is not failed for refusing a box it never
#' promised. The composition a caller writes instead is
#' `prox(v, t) = pmin(pmax(prox_penalty(v, t), lower), upper)`, which imposes
#' both terms exactly.
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
#' extrapolation with its restart. Box bounds are refused here, with a
#' message naming `prox` as where the constraint belongs.
#'
#' @param optimizer A `ProxGrad` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. A finite `lower`
#'   or `upper` raises an error; `he` is accepted and ignored.
#'
#' @return An [optimizer_result()] whose `value` is the **total** objective
#'   \eqn{f(x) + g(x)}, which is why `g` is required alongside `prox`, and
#'   whose `gradient` is the proximal gradient mapping at `par` rather than
#'   \eqn{\nabla f}.
#'
#' @keywords internal
S7::method(minimize, ProxGrad) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
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
#' extrapolation with its restart. Counts its own evaluations, so the caller
#' need not.
#'
#' @details
#' This is the one method in the package written in R instead of compiled.
#' Every iteration calls the objective, its gradient and the proximal
#' operator, all three R functions supplied by the caller, and the loop
#' around them costs a fraction of a microsecond against those. Compiling it
#' would move the callbacks and change nothing else.
#'
#' The stationarity measure is read **at the iterate** and not at the
#' extrapolated point. With momentum the two differ, and reading the
#' extrapolated one leaves a mapping that never vanishes.
#'
#' @param optimizer A `ProxGrad` object.
#' @param spec The objective handle from [as_objective()].
#' @param par The starting point, a numeric vector.
#'
#' @return A list in the shape [build_result()] consumes: the point, the
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
    function(x) {
      n_grad <<- n_grad + 1L
      fd_gradient(spec$fn, x)
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

  for (it in seq_len(as.integer(optimizer@maxit))) {
    # A restart retries WITHIN the iteration. Spending an iteration on it
    # would make the accelerated method look slower than the plain one for
    # no reason but the bookkeeping, and would leave a gap in the trace.
    # The k > 1 guard makes the retry provably finite: at k = 1 there is no
    # momentum to discard, and the backtracked step decreases the objective.
    repeat {
      fy <- f(y)
      gy <- gradf(y)

      # backtracking: shrink until the quadratic model at y dominates f at
      # the candidate, which is the inequality the convergence proof uses
      repeat {
        x_new <- as.numeric(prox(y - t * gy, t))
        d <- x_new - y
        f_cand <- f(x_new)
        if (!is.finite(f_cand) ||
            f_cand > fy + sum(gy * d) + sum(d^2) / (2 * t)) {
          t <- t * optimizer@shrink
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
      } else {
        break
      }
    }

    x_old <- x
    f_prev <- fx_total
    x <- x_new
    fx_total <- f_new_total

    # The mapping must be read at the ITERATE. Without momentum y equals x
    # and the quantity above already is it; with momentum the two differ by
    # the extrapolation, and a mapping read at y does not vanish at the
    # solution -- the run would then circle the answer reporting failure.
    # Re-reading it at x costs one gradient and one prox per iteration, which
    # the accelerated method repays in iterations saved.
    if (optimizer@accelerate) {
      gx <- gradf(x)
      gmap <- (x - as.numeric(prox(x - t * gx, t))) / t
    }

    state <- list(x_new = x, x_old = x_old, f_new = fx_total,
                  f_old = f_prev, gradient = gmap)
    met <- crit_met(optimizer@criterion, state)

    if (optimizer@keep_trace) {
      trace[[length(trace) + 1L]] <- data.frame(
        iteration = it, value = fx_total, step = t,
        gradient = max(abs(gmap)))
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
#' The finite-difference gradient the R-level proximal loop uses when the
#' caller supplies none, with the step scaled to each coordinate.
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
