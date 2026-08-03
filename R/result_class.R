#' @include optimizer_class.R
NULL

#' @title S7 Class for the Result of an Optimisation
#'
#' @description
#' What \code{\link{minimize}} returns: the answer, how it was reached, and
#' enough of the run to diagnose it when it was not reached.
#'
#' @param par The minimiser.
#' @param value The objective there.
#' @param gradient The gradient there, or \code{NULL} if the method computes none.
#' @param counts Evaluations of the objective and of the gradient.
#' @param iterations Iterations performed.
#' @param converged Logical; see Details.
#' @param criterion_met Which rule ended the run.
#' @param message A human-readable account.
#' @param trace The iteration path, or \code{NULL}.
#' @param optimizer The \code{\link{optimizer}} that produced this, kept so the
#'   run can be repeated exactly.
#' @param elapsed Seconds.
#' @param seed The state of the random number generator when the run began, for
#'   a method that draws any, and \code{NULL} otherwise.
#'
#' @details
#' \code{converged} is \code{TRUE} only when the stopping rule was satisfied. It
#' is \strong{never} \code{TRUE} because the iteration budget ran out — that is
#' the commonest defect in hand-written optimisation loops, and it turns a
#' failure into a silently wrong answer.
#'
#' \code{trace}, present when the optimiser was built with
#' \code{keep_trace = TRUE}, records for each iteration the objective, the
#' quantity the criterion is watching, the step actually taken, and the name of
#' any safeguard that fired. That last column is what turns "it did not
#' converge" into a diagnosis.
#'
#' \code{seed} is filled in by the methods that draw random numbers -- a
#' \code{"mads"} poll, a subsampling \code{\link{adam}}, a
#' \code{\link{multistart}} generating its own starts. Assigning it back with
#' \code{assign(".Random.seed", res@seed, globalenv())} reproduces the run
#' exactly. A stochastic method that cannot be repeated is very hard to debug,
#' and remembering to call \code{set.seed()} beforehand is not something anyone
#' does until the second time they need it.
#'
#' @return An S7 object of class \code{optimizer_result}.
#'
#' @examples
#' res <- minimize(bfgs(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'                 gr = function(p) 2 * (p - c(1, 2)))
#' res@par
#' res@converged
#' res@criterion_met
#' res@counts
#'
#' @seealso \code{\link{minimize}}
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
#' Renders a time in seconds using the unit its magnitude calls for:
#' microseconds below a millisecond, milliseconds below a second, seconds below
#' a minute, minutes and seconds below an hour, hours and minutes above.
#'
#' @param sec A single non-negative number of seconds.
#' @return A character string, or \code{NA_character_} when \code{sec} is
#'   missing or not finite.
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


#' @title Print Method for an Optimisation Result
#' @name print.optimizer_result
#' @description
#' Prints the objective value, the leading parameters, the evaluation counts,
#' the elapsed time and the convergence status.
#' @param x An \code{\link{optimizer_result}}.
#' @param digits Decimal places the parameters are rounded to. Defaults to 4.
#' @param max_par How many parameters to show; any remainder is summarised as
#'   a count. Defaults to 6.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' res <- minimize(gd(), function(p) sum((p - 1:2)^2), c(0, 0))
#' print(res)
#' print(res, digits = 2, max_par = 1)
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


#' @title Summary Method for an Optimisation Result
#' @name summary.optimizer_result
#' @param object An \code{\link{optimizer_result}}.
#' @param ... Unused.
#' @return \code{object}, invisibly. Called for the printed summary.
#' @examples
#' res <- minimize(gd(keep_trace = TRUE),
#'                 function(p) sum((p - 1:2)^2), c(0, 0))
#' summary(res)
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


#' @title Plot Method for an Optimisation Result
#' @name plot.optimizer_result
#'
#' @description
#' The objective against iteration, with any iteration at which a safeguard
#' fired marked. Requires \code{keep_trace = TRUE}.
#'
#' @param x An \code{\link{optimizer_result}}.
#' @param ... Passed to \code{\link[graphics]{plot}}.
#'
#' @return No return value; called for the plot.
#'
#' @examples
#' res <- minimize(gd(keep_trace = TRUE),
#'                 function(p) sum((p - 1:2)^2), c(0, 0))
#' plot(res)
#'
#' @importFrom graphics plot points grid legend
#' @keywords internal
S7::method(plot, optimizer_result) <- function(x, ...) {
  tr <- x@trace
  if (is.null(tr) || !nrow(tr)) {
    stop("No trace to plot; build the optimiser with keep_trace = TRUE.",
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
