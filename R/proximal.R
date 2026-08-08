#' @include optimizer_class.R
NULL

#' @title S7 Class for the Proximal Gradient Method
#' @description The class \code{\link{prox_grad}} instantiates.
#' @param prox The proximal operator of the non-smooth part.
#' @param g The value of the non-smooth part.
#' @param accelerate Whether the momentum extrapolation is applied.
#' @param step The initial step length.
#' @param shrink The backtracking factor.
#' @param restart Whether an increase in the objective resets the momentum.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{prox_grad}}
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
#' \code{accelerate = TRUE} the method is the accelerated one of Beck and
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
#' \subsection{What the stopping rule reads}{
#' The gradient of the total objective does not vanish at the solution --
#' that is what non-differentiable means -- so this method reports the
#' \strong{proximal gradient mapping}
#' \deqn{G_t(x) = \frac{x - \mathrm{prox}_{tg}(x - t\nabla f(x))}{t}}
#' as its gradient, read at the iterate. It vanishes exactly at a
#' stationary point of \eqn{f + g} and reduces to \eqn{\nabla f} when
#' \eqn{g} is absent, so \code{\link{crit_grad}} keeps its meaning and
#' its default tolerance. The accelerated variant pays one extra gradient
#' per iteration for it, the extrapolated point at which it takes its step
#' not being the point it reports.
#' }
#'
#' \subsection{Restarting}{
#' Momentum makes the objective non-monotone, and an increase far from the
#' solution is a symptom of momentum built in the wrong direction. With
#' \code{restart = TRUE} an increase resets the extrapolation to the
#' current point, which is the adaptive restart of O'Donoghue and Candes
#' and costs one comparison per iteration.
#' }
#'
#' @param prox The proximal operator of the non-smooth part, a function of
#'   the point and the step length, \code{prox(v, step)}, returning the
#'   minimizer of \eqn{\lVert b - v \rVert^2/(2\,\mathrm{step}) + g(b)}.
#'   \code{\link[penalties7]{penalty_prox}} supplies one for every penalty
#'   that has it.
#' @param g The value of the non-smooth part, a function of the point.
#'   Required alongside \code{prox}: the two describe the same term, and
#'   without \code{g} the reported objective would be the smooth part
#'   alone.
#' @param accelerate Apply the momentum extrapolation? Defaults to
#'   \code{TRUE}.
#' @param step The initial step length offered to the backtracking search.
#'   Defaults to \code{1}.
#' @param shrink The factor by which a rejected step is reduced. Defaults
#'   to \code{0.5}.
#' @param restart Reset the momentum when the objective increases?
#'   Defaults to \code{TRUE}, and is ignored when \code{accelerate} is
#'   \code{FALSE}.
#' @param criterion The stopping rule; see \code{\link{crit_any}}.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to \code{Inf}.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @return An S7 object of class \code{\link{ProxGrad}}, to be handed to
#'   \code{\link{minimize}}.
#'
#' @references
#' Beck, A. and Teboulle, M. (2009). A fast iterative shrinkage-thresholding
#' algorithm for linear inverse problems. \emph{SIAM Journal on Imaging
#' Sciences}, 2(1), 183--202.
#'
#' O'Donoghue, B. and Candes, E. (2015). Adaptive restart for accelerated
#' gradient schemes. \emph{Foundations of Computational Mathematics}, 15(3),
#' 715--732.
#'
#' @seealso \code{\link{bundle}} for a non-smooth method that needs no
#'   proximal operator, \code{\link{gd}} for the smooth case.
#'
#' @examples
#' # a lasso-penalized least squares problem, solved through the operator
#' set.seed(1)
#' X <- matrix(rnorm(200 * 8), 200, 8)
#' b0 <- c(2, -1.5, 0, 0, 0.8, 0, 0, 0)
#' y <- as.numeric(X %*% b0 + rnorm(200))
#' lambda <- 0.4
#'
#' fit <- minimize(
#'   prox_grad(prox = function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0),
#'             g = function(b) lambda * sum(abs(b))),
#'   fn = function(b) sum((y - X %*% b)^2) / (2 * nrow(X)),
#'   gr = function(b) -crossprod(X, y - X %*% b) / nrow(X),
#'   par = rep(0, 8))
#' round(fit@par, 3)
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


S7::method(optimizer_bounded, ProxGrad) <- function(optimizer) FALSE


#' @title Minimize by the Proximal Gradient Method
#' @name minimize.ProxGrad
#' @description Runs \code{\link{prox_grad}} on the objective.
#' @param optimizer A \code{ProxGrad} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
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
#' The iteration behind \code{\link{prox_grad}}: a backtracked gradient
#' step on the smooth part, the proximal operator applied to its result,
#' and the momentum extrapolation with its restart.
#'
#' @details
#' Written in R rather than compiled, because every iteration calls the
#' objective, its gradient and the proximal operator, all of which are R
#' functions supplied by the caller; the loop around them costs a fraction
#' of a microsecond against those.
#'
#' @param optimizer A \code{ProxGrad} object.
#' @param spec The objective handle from \code{\link{as_objective}}.
#' @param par The starting point.
#'
#' @return A list in the shape \code{\link{build_result}} consumes.
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

      if (optimizer@accelerate && optimizer@restart && k > 1L &&
          f_new_total > fx_total) {
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
