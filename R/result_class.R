#' @include optimizer_class.R
NULL

#' @title S7 Class for the Result of an Optimization
#'
#' @description
#' What [minimize()] and [maximize()] return: the point reached and the value
#' there, the evaluation counts, whether a stopping rule was satisfied and
#' which one, and enough of the run to diagnose it when none was. The
#' optimizer itself is kept on the object, so a run can be repeated from what
#' it returned.
#'
#' @details
#' # `converged` means a rule fired
#'
#' `converged` is `TRUE` only when the stopping rule was satisfied. Running
#' out of iterations or of evaluations leaves it `FALSE` and puts
#' `iteration budget reached` or `evaluation budget exhausted` into
#' `criterion_met`. Reporting a budget as a success is the commonest defect
#' in a hand-written optimization loop; it turns a failure into a wrong
#' answer that nothing downstream can detect.
#'
#' # The trace, and the column that changes name
#'
#' With `keep_trace = TRUE` on the optimizer, `trace` is a data frame with
#' one row per iteration and five columns: `iteration`, `value`, the
#' quantity the stopping rule is watching, `step` and `safeguard`. The third
#' column is `gnorm` for a method that computes a gradient and
#' `stationarity` for one that does not, so code reading a trace should ask
#' for the column by position or check `names(trace)`.
#'
#' `safeguard` names the repair the algorithm applied to its own step, or
#' `"none"`. The names differ by method and are worth reading: `step
#' shortened` and `step adjusted` for the line-search methods, `cg restart`
#' when conjugacy is lost, `bb curvature reset` when a secant pair reports
#' none, and `reflect`, `expand`, `contract in` and `contract out` for the
#' simplex. `summary()` tabulates them.
#'
#' # Repeating a stochastic run
#'
#' `seed` holds `.Random.seed` as it stood when the run began, and is filled
#' in only by the methods that draw: [sa()], a resampling [adam()], a
#' [multistart()] generating its own starts. Assigning it back reproduces
#' the run exactly:
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
#' @param criterion_met A single string: the label of the rule that fired, or
#'   `iteration budget reached` or `evaluation budget exhausted` when the run
#'   ended without one.
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
#'   receives one rather than constructing it.
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
#' # A budget stops the run and leaves converged FALSE, with the reason.
#' short <- minimize(bfgs(maxit = 3), f, c(-1.2, 1), gr = g)
#' c(short@converged, short@criterion_met)
#'
#' # The trace names its third column after what the rule watches.
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
#' Renders a time in seconds using the unit its magnitude calls for, to three
#' significant figures: microseconds below a millisecond, milliseconds below
#' a second, seconds below a minute, whole minutes and seconds below an hour,
#' whole hours and minutes above. Used by
#' [print.optimizer_result()] for the `elapsed` line.
#'
#' @param sec A single number of seconds.
#'
#' @return A character string such as `"250 ms"`, `"1 min 30 s"` or
#'   `"1 h 7 min"`. `NA_character_` when `sec` is empty, missing or not
#'   finite, so an unmeasured duration is reported as unmeasured.
#'
#' @examples
#' vapply(c(1e-5, 5e-4, 0.25, 12, 90, 4000), format_elapsed, "")
#' format_elapsed(NA)
#'
#' @keywords internal
format_elapsed <- function(sec) {
  if (!length(sec) || !is.finite(sec)) return(NA_character_)
  if (sec < 1e-3) return(sprintf("%.3g us", sec * 1e6))
  if (sec < 1)    return(sprintf("%.3g ms", sec * 1e3))
  if (sec < 60)   return(sprintf("%.3g s", sec))
  if (sec < 3600) {
    m <- floor(sec / 60)
    return(sprintf("%d min %.0f s", m, sec - 60 * m))
  }
  h <- floor(sec / 3600)
  sprintf("%d h %.0f min", h, (sec - 3600 * h) / 60)
}


#' @title Print Method for an Optimization Result
#' @name print.optimizer_result
#'
#' @description
#' Shows the run in six lines at most: the method's name, the objective
#' value, the leading parameters, the iteration and evaluation counts, the
#' elapsed time and the convergence status with the rule that fired. A
#' non-empty `message` adds a `note` line. A failure prints `NO` in capitals,
#' so a run that did not converge cannot be skimmed past.
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
#' # A run stopped by its budget says so on the converged line.
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
#' Everything [print.optimizer_result()] shows, followed by a count of each
#' safeguard the run applied to its own steps. This is the view that answers
#' *why* a run behaved as it did: a Newton run that shortened four steps was
#' repairing an indefinite Hessian, and a Barzilai-Borwein run resetting its
#' curvature had secant pairs carrying none.
#'
#' @details
#' The safeguard table needs a trace, so the optimizer must have been built
#' with `keep_trace = TRUE`. Without one the method prints what
#' [print.optimizer_result()] prints and stops there. With a trace in which
#' nothing fired it says `safeguards : none fired`, which is information
#' rather than silence.
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
#' # Newton on the curved valley shortens a few steps on the way in.
#' summary(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g))
#'
#' # A quadratic gives it no trouble at all.
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
#' point at every iteration where a safeguard fired and a legend saying so.
#' The method's name is the title. Reading the marks against the curve shows
#' where the algorithm was in trouble and whether the objective moved when it
#' was.
#'
#' @details
#' The plot is built from the `trace`, so the optimizer must have been
#' created with `keep_trace = TRUE`. Without one the method stops with a
#' message naming that argument, since an empty plot would say nothing.
#'
#' The vertical axis is the objective on its own scale; a run whose objective
#' falls over several orders of magnitude is better read with `log = "y"`,
#' which passes through to [graphics::plot()] like any other argument.
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
#' # The marked iterations are the ones where a step had to be repaired.
#' plot(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g), log = "y")
#'
#' # Without a trace there is nothing to draw, and the method says which
#' # argument was missing.
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
