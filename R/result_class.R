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


#' @title Print Method for an Optimisation Result
#' @name print.optimizer_result
#' @param x An \code{\link{optimizer_result}}.
#' @param digits Significant digits.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' res <- minimize(gd(), function(p) sum((p - 1:2)^2), c(0, 0))
#' print(res)
#' @keywords internal
S7::method(print, optimizer_result) <- function(x, digits = 6, ...) {
  cat("<optimizer_result> ", x@optimizer@name, "\n", sep = "")
  cat("  value      : ", format(x@value, digits = digits), "\n", sep = "")
  cat("  par        : ", paste(format(x@par, digits = digits), collapse = " "),
      "\n", sep = "")
  cat("  iterations : ", x@iterations,
      "   evaluations: f ", x@counts[["f"]], ", g ", x@counts[["g"]],
      "\n", sep = "")
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
