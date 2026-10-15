#' @include optimizer_class.R
NULL

#' @title S7 Class for the Result of an Optimization
#'
#' @description
#' The value of [minimize()] and [maximize()]: the point reached and the
#' value there, the evaluation counts, whether the run converged and which
#' rule stopped it, and enough of the run to diagnose it otherwise. The
#' optimizer itself is kept on the object, so a run can be repeated from its
#' result.
#'
#' @details
#' # The meaning of `converged`
#'
#' `converged` is `TRUE` when the stopping rule was satisfied, and also when
#' a line search with a declared resolution established that no step
#' decreases the objective by more than that resolution (`criterion_met` is
#' then `no decrease above the objective's resolution`) or when [brent()]
#' reduced its bracket to the resolution of the iterate (`bracket at machine
#' resolution`). Running out of iterations or evaluations leaves it `FALSE`
#' and puts `iteration budget reached` or `evaluation budget reached` into
#' `criterion_met`; `message` then reads `evaluation budget exhausted` for
#' the second. A line search that finds no acceptable step also leaves it
#' `FALSE`, with `criterion_met` equal to `stopped without converging`. In
#' the methods with a line search, an accepted step that moves the iterate
#' by no more than \eqn{4\epsilon\max(1, \lVert x\rVert_\infty)}, with
#' \eqn{\epsilon} the machine precision, is treated in the same way: the
#' stopping rule is then evaluated without the previous iterate, so only a
#' rule on the gradient can report convergence, and otherwise `message`
#' reads `the step did not move the iterate beyond rounding` and the trace
#' records `step below rounding`. A change in the parameters or in the
#' objective over such a step is rounding, and it is not read as
#' convergence. For the same reason, when the stopping rule is met on a step
#' \eqn{\Delta x} along which the first-order change at the new point,
#' \eqn{\lvert g^\top \Delta x\rvert}, exceeds \eqn{\max(1, \lvert f\rvert)},
#' the step is undone and the run is handled as above, with the message
#' `the step that met the stopping rule is not resolved by the objective`
#' and the safeguard `step not resolved`. A run that ends on a budget is therefore reported as not converged, so
#' that code reading the result can distinguish it from a run that
#' satisfied a rule.
#'
#' # The trace
#'
#' With `keep_trace = TRUE` on the optimizer, `trace` is a data frame with
#' one row per iteration. For a single method it has the columns
#' `iteration`, `value`, a third column holding the quantity read by the
#' stopping rule, `step` and `safeguard`; [brent()] has no `step` column. The
#' third column is `gnorm` for a method that computes a gradient and
#' `stationarity` for one that does not, so code that reads a trace should
#' select the column by position or check `names(trace)`. A [multistart()]
#' trace has one row per start, with the columns `start`, `value`,
#' `converged` and `iterations`, and a [chain()] trace has an additional
#' `stage` column.
#'
#' `safeguard` names the repair that the algorithm applied to its own step,
#' or `"none"`. The names differ by method; examples are `step shortened`
#' and `step adjusted` for the line-search methods, `cg restart` when
#' conjugacy is lost, `bb curvature reset` when a secant pair has no usable
#' curvature, and `reflect`, `expand`, `contract in` and `contract out` for
#' the simplex. `summary()` tabulates them.
#'
#' # Repeating a stochastic run
#'
#' `seed` holds `.Random.seed` as it was when the run began, and is filled in
#' only by the methods that draw random numbers themselves: [sa()],
#' [compass()] with its default `directions = "mads"`, and [multistart()].
#' A [chain()] carries the seed of its first stage. [adam()] draws no random
#' numbers, so its `seed` is `NULL`; a resampling objective draws in the
#' session of the caller. Assigning the seed back reproduces the run
#' exactly:
#'
#' ```r
#' a <- minimize(sa(maxit = 300), f, c(-1.2, 1))
#' assign(".Random.seed", a@seed, globalenv())
#' b <- minimize(sa(maxit = 300), f, c(-1.2, 1))
#' identical(a@par, b@par)   # TRUE
#' ```
#'
#' For a deterministic method `seed` is `NULL`.
#'
#' @param par The minimizer reached, a numeric vector of the same length as
#'   the starting point.
#' @param value The objective at `par`, a single number. For a run started by
#'   [maximize()] this is the value of the original objective, with the sign
#'   already undone.
#' @param gradient The gradient at `par`, a numeric vector, or `NULL` for a
#'   method that computes none. [prox_grad()] reports the proximal gradient
#'   mapping here.
#' @param counts A named integer vector of length three, `f`, `g` and `h`:
#'   evaluations of the objective, the gradient and the Hessian. A gradient
#'   obtained by differencing is counted in `f`, one per coordinate per
#'   side.
#' @param iterations The number of iterations performed, a single integer.
#' @param converged A single logical; see Details.
#' @param criterion_met A single string: the label of the rule that fired,
#'   `no decrease above the objective's resolution` or `bracket at machine
#'   resolution` for the two stops described in Details, `iteration budget
#'   reached` or `evaluation budget reached` for a budget, and `stopped
#'   without converging` for a run that stopped for another reason.
#' @param message A single string, empty when there is nothing to report.
#'   Carries notes such as `gradient obtained by finite differences`.
#' @param trace The iteration path as a data frame, or `NULL` when the
#'   optimizer was built with `keep_trace = FALSE`, which is the default.
#' @param optimizer The [optimizer()] that produced the result, kept whole so
#'   the run can be repeated or restarted from where it stopped.
#' @param elapsed Wall-clock seconds, a single number. Measured to the
#'   platform's clock resolution, so a fast run can read exactly `0`.
#' @param seed The state of the random number generator when the run began,
#'   an integer vector, or `NULL` for a deterministic method.
#'
#' @return An S7 object of class `optimizer_result` carrying the twelve
#'   properties above. Built by the methods of [minimize()]; a caller
#'   receives one and does not construct it.
#'
#' @examples
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#'
#' res <- minimize(bfgs(keep_trace = TRUE), f, c(-1.2, 1), gr = g)
#' res@par
#' res@converged
#' res@criterion_met
#' res@counts
#'
#' # A budget stops the run and leaves converged FALSE, with the reason in
#' # criterion_met.
#' short <- minimize(bfgs(maxit = 3), f, c(-1.2, 1), gr = g)
#' c(short@converged, short@criterion_met)
#'
#' # The name of the third column of the trace is the quantity read by the
#' # rule.
#' names(res@trace)
#' names(minimize(nelder_mead(keep_trace = TRUE), f, c(-1.2, 1))@trace)
#'
#' # The optimizer travels with the result, so the run can be repeated.
#' again <- minimize(res@optimizer, f, c(-1.2, 1), gr = g)
#' all.equal(again@par, res@par)
#'
#' @seealso [minimize()], [print.optimizer_result()] and
#'   [summary.optimizer_result()] for the two views of it,
#'   [plot.optimizer_result()] for the trace.
#' @export
optimizer_result <- S7::new_class(
  "optimizer_result",
  properties = list(
    par           = S7::class_numeric,
    value         = S7::class_numeric,
    gradient      = S7::class_any,
    counts        = S7::class_any,
    iterations    = S7::class_numeric,
    converged     = S7::class_logical,
    criterion_met = S7::class_character,
    message       = S7::class_character,
    trace         = S7::class_any,
    optimizer     = S7::class_any,
    elapsed       = S7::class_numeric,
    seed          = S7::class_any
  )
)


#' Format a Duration With a Unit Matched to Its Size
#'
#' @description
#' Renders a time in seconds in a unit matched to its size: microseconds
#' below a millisecond, milliseconds below a second and seconds below a
#' minute, each to three significant figures, then whole minutes and seconds
#' below an hour, and whole hours and minutes above. The unit is chosen on
#' the rounded value, so that no field reaches its carry (`"59 min 60 s"`).
#' Used by [print.optimizer_result()] for the `elapsed` line.
#'
#' @param sec A single number of seconds.
#'
#' @return A character string such as `"250 ms"`, `"1 min 30 s"` or
#'   `"1 h 7 min"`. `NA_character_` when `sec` is empty, missing or not
#'   finite, so an unmeasured duration is reported as unmeasured.
#'
#' @keywords internal
format_elapsed <- function(sec) {
  if (!length(sec) || !is.finite(sec)) return(NA_character_)
  # each unit is chosen on the value as it will be printed, so that rounding
  # cannot produce "1e+03 ms", "60 s" or "59 min 60 s"
  if (signif(sec * 1e6, 3) < 1000) return(sprintf("%.3g us", sec * 1e6))
  if (signif(sec * 1e3, 3) < 1000) return(sprintf("%.3g ms", sec * 1e3))
  if (signif(sec, 3) < 60) return(sprintf("%.3g s", sec))
  s <- round(sec)
  if (s < 3600) return(sprintf("%.0f min %.0f s", s %/% 60, s %% 60))
  m <- round(sec / 60)
  sprintf("%.0f h %.0f min", m %/% 60, m %% 60)
}


#' @title Print Method for an Optimization Result
#' @name print.optimizer_result
#'
#' @description
#' Shows the run in six lines: the name of the method, the objective value,
#' the leading parameters, the iteration and evaluation counts, the elapsed
#' time and the convergence status with the rule that fired. The elapsed line
#' is omitted when the time is missing, and a non-empty `message` adds a
#' seventh line, `note`. A run that did not converge prints `NO` in
#' capitals.
#'
#' @param x An [optimizer_result()].
#' @param digits Decimal places the parameters are rounded to. A single
#'   non-negative whole number, default 4. Anything else raises an error.
#'   The objective value is not affected: it always prints to six
#'   significant figures.
#' @param max_par How many parameters to show. A single positive whole
#'   number, default 6; the remainder is reported as
#'   `... (6 of 40 shown)`.
#' @param ... Unused.
#'
#' @return `x`, invisibly. Called for the output.
#'
#' @examples
#' res <- minimize(gd(), function(p) sum((p - 1:2)^2), c(0, 0))
#' res
#' print(res, digits = 2, max_par = 1)
#'
#' # The converged line of a run stopped by its budget names the budget.
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' print(minimize(gd(maxit = 5), rosen, c(-1.2, 1)))
#'
#' @keywords internal
S7::method(print, optimizer_result) <- function(x, digits = 4, max_par = 6,
                                                ...) {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits < 0 || digits != round(digits)) {
    stop("'digits' must be a single non-negative whole number.", call. = FALSE)
  }
  if (!is.numeric(max_par) || length(max_par) != 1L || is.na(max_par) ||
      max_par < 1 || max_par != round(max_par)) {
    stop("'max_par' must be a single positive whole number.", call. = FALSE)
  }
  p <- length(x@par)
  shown <- round(x@par[seq_len(min(p, max_par))], digits)
  tail_note <- if (p > max_par) {
    sprintf(" ... (%d of %d shown)", as.integer(max_par), p)
  } else ""

  cat("<optimizer_result> ", x@optimizer@name, "\n", sep = "")
  cat("  value      : ", format(x@value, digits = 6), "\n", sep = "")
  cat("  par        : ", paste(format(shown), collapse = " "), tail_note,
      "\n", sep = "")
  cat("  iterations : ", x@iterations,
      "   evaluations: f ", x@counts[["f"]], ", g ", x@counts[["g"]],
      "\n", sep = "")
  el <- format_elapsed(x@elapsed)
  if (!is.na(el)) cat("  elapsed    : ", el, "\n", sep = "")
  if (x@converged) {
    cat("  converged  : yes (", x@criterion_met, ")\n", sep = "")
  } else {
    cat("  converged  : NO (", x@criterion_met, ")\n", sep = "")
  }
  if (nzchar(x@message)) cat("  note       : ", x@message, "\n", sep = "")
  invisible(x)
}


#' @title Summary Method for an Optimization Result
#' @name summary.optimizer_result
#'
#' @description
#' Prints what [print.optimizer_result()] shows, followed by a count of each
#' safeguard that the run applied to its own steps. The counts show how the
#' run proceeded: a Newton run that reports `step shortened` had steps
#' shortened by the line search, one that reports `hessian modified` was
#' repairing an indefinite Hessian, and a Barzilai-Borwein run that reports
#' `bb curvature reset` had secant pairs without usable curvature.
#'
#' @details
#' The safeguard table needs a trace, so the optimizer must have been built
#' with `keep_trace = TRUE`. Without a trace the method prints what
#' [print.optimizer_result()] prints. With a trace in which no safeguard
#' fired, the output reads `safeguards : none fired`.
#'
#' @param object An [optimizer_result()].
#' @param ... Unused.
#'
#' @return `object`, invisibly. Called for the printed summary.
#'
#' @examples
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#'
#' # Newton on the curved valley: the line search shortens a few steps.
#' summary(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g))
#'
#' # On a quadratic no safeguard fires.
#' summary(minimize(newton(keep_trace = TRUE),
#'                  function(p) sum((p - 1:2)^2), c(0, 0)))
#'
#' @keywords internal
S7::method(summary, optimizer_result) <- function(object, ...) {
  print(object)
  tr <- object@trace
  if (!is.null(tr) && nrow(tr)) {
    fired <- tr$safeguard[tr$safeguard != "none"]
    if (length(fired)) {
      cat("  safeguards :\n")
      tab <- table(fired)
      for (nm in names(tab)) {
        cat("    ", nm, ": ", tab[[nm]], "\n", sep = "")
      }
    } else {
      cat("  safeguards : none fired\n")
    }
  }
  invisible(object)
}


#' @title Plot Method for an Optimization Result
#' @name plot.optimizer_result
#'
#' @description
#' Draws the objective against the iteration number as a line, with a filled
#' point at every iteration where a safeguard fired and a legend for the
#' points. The title is the name of the method. The marks show the
#' iterations at which the algorithm repaired its step and how the objective
#' moved there.
#'
#' @details
#' The plot is built from the `trace`, so the optimizer must have been
#' created with `keep_trace = TRUE`. Without a trace the method stops with an
#' error that names that argument.
#'
#' The vertical axis is the objective on its own scale; for a run whose
#' objective falls over several orders of magnitude, `log = "y"` is more
#' readable, and it is passed to [graphics::plot()] like any other
#' argument.
#'
#' @param x An [optimizer_result()].
#' @param ... Passed to [graphics::plot()]. `type`, `lwd`, `las`, `xlab`,
#'   `ylab` and `main` are already set and passing them again is an error, as
#'   it is for any duplicated argument.
#'
#' @return `NULL`, invisibly. Called for the plot.
#'
#' @examples
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#'
#' # The marked iterations are those at which a step was repaired.
#' plot(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g), log = "y")
#'
#' # Without a trace there is nothing to draw, and the error names the
#' # missing argument.
#' try(plot(minimize(newton(), f, c(-1.2, 1), gr = g)))
#'
#' @importFrom graphics plot points grid legend
#' @keywords internal
S7::method(plot, optimizer_result) <- function(x, ...) {
  tr <- x@trace
  if (is.null(tr) || !nrow(tr)) {
    stop("No trace to plot; build the optimizer with keep_trace = TRUE.",
         call. = FALSE)
  }
  graphics::plot(tr$iteration, tr$value, type = "l", lwd = 2, las = 1,
                 xlab = "iteration", ylab = "objective",
                 main = x@optimizer@name, ...)
  graphics::grid()
  hit <- tr$safeguard != "none"
  if (any(hit)) {
    graphics::points(tr$iteration[hit], tr$value[hit], pch = 19, col = "#9C3E11")
    graphics::legend("topright", pch = 19, col = "#9C3E11", bty = "n",
                     legend = "safeguard fired")
  }
  invisible(NULL)
}
