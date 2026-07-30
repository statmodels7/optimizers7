#' @include criterion.R
NULL

#' @title S7 Class for Optimisation Algorithms
#'
#' @description
#' An optimiser is an object carrying an algorithm and every setting that
#' algorithm obeys: its stopping rule, its budgets, what it reports, and what it
#' keeps. It is created once and may be reused, inspected, stored beside a result
#' and passed around; nothing about a run is hidden in the call that started it.
#'
#' @details
#' The class is abstract. Each algorithm is a subclass with its own constructor —
#' \code{\link{gradient_descent}}, and in time \code{newton()}, \code{bfgs()} and
#' the rest — adding only the settings that are genuinely its own.
#'
#' These properties are shared by every algorithm without exception, which is
#' what lets a caller swap one optimiser for another without changing anything
#' else.
#'
#' @param name A short character name, used when printing and reporting.
#' @param criterion The stopping rule, a \code{\link{criterion}} object.
#' @param maxit Maximum number of iterations.
#' @param max_eval Maximum number of objective evaluations. A budget on work
#'   rather than on progress: a line search can spend many evaluations in one
#'   iteration, and a runaway one is invisible to \code{maxit}.
#' @param verbose Logical; whether to report progress.
#' @param refresh Report every \code{refresh} iterations. \code{0} reports only
#'   the final summary.
#' @param keep_trace Logical; whether to store the iteration path in the result.
#'
#' @return An S7 object of class \code{optimizer}.
#'
#' @seealso \code{\link{minimize}}, \code{\link{gradient_descent}}
#' @export
optimizer <- S7::new_class(
  "optimizer",
  properties = list(
    name       = S7::class_character,
    criterion  = S7::class_any,
    maxit      = S7::class_numeric,
    max_eval   = S7::class_numeric,
    verbose    = S7::class_logical,
    refresh    = S7::class_numeric,
    keep_trace = S7::class_logical
  ),
  abstract = TRUE,
  validator = function(self) {
    if (!S7::S7_inherits(self@criterion, criterion)) {
      return("Property 'criterion' must be a 'criterion' object.")
    }
    if (length(self@maxit) != 1L || is.na(self@maxit) || self@maxit < 1) {
      return("Property 'maxit' must be a single positive number.")
    }
    if (length(self@refresh) != 1L || is.na(self@refresh) || self@refresh < 0) {
      return("Property 'refresh' must be a single non-negative number.")
    }
    NULL
  }
)


#' Validate the Settings Every Optimiser Shares
#'
#' @description
#' The checks each constructor would otherwise repeat, in one place and in one
#' wording.
#'
#' @param criterion The stopping rule.
#' @param maxit,max_eval,refresh Numeric budgets.
#' @param verbose,keep_trace Logical flags.
#'
#' @return Invisibly \code{TRUE}; raises an error otherwise.
#'
#' @keywords internal
check_optimizer_args <- function(criterion, maxit, max_eval, verbose, refresh,
                                 keep_trace) {
  if (!S7::S7_inherits(criterion, criterion_class())) {
    stop("'criterion' must be a criterion object, e.g. crit_grad().",
         call. = FALSE)
  }
  pos <- function(v, nm) {
    if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 1) {
      stop("'", nm, "' must be a single positive number.", call. = FALSE)
    }
  }
  pos(maxit, "maxit")
  pos(max_eval, "max_eval")
  if (length(refresh) != 1L || !is.numeric(refresh) || is.na(refresh) ||
      refresh < 0) {
    stop("'refresh' must be a single non-negative number.", call. = FALSE)
  }
  flag <- function(v, nm) {
    if (length(v) != 1L || !is.logical(v) || is.na(v)) {
      stop("'", nm, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }
  flag(verbose, "verbose")
  flag(keep_trace, "keep_trace")
  invisible(TRUE)
}

# The `criterion` class object, fetched rather than captured, so that the check
# above cannot be fooled by the class being re-created -- the trap recorded in
# linkfunctions7, where comparing S7 classes with identical() silently degraded
# a numerical fallback.
criterion_class <- function() criterion


#' @title Print Method for Optimisers
#' @name print.optimizer
#' @param x An \code{\link{optimizer}} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' print(gradient_descent())
#' @keywords internal
S7::method(print, optimizer) <- function(x, ...) {
  cat("<optimizer> ", x@name, "\n", sep = "")
  cat("  stop when : ", x@criterion@label, "\n", sep = "")
  cat("  budgets   : maxit ", x@maxit, ", evaluations ", x@max_eval, "\n", sep = "")
  extra <- setdiff(names(S7::props(x)),
                   c("name", "criterion", "maxit", "max_eval", "verbose",
                     "refresh", "keep_trace"))
  if (length(extra)) {
    # A setting need not be a number: a line search is an object, and format()
    # on one is an error rather than a string. Anything carrying a label shows
    # its label, anything atomic shows its value, and anything else shows what
    # it is -- so adding a new kind of setting can never break printing.
    describe <- function(p) {
      v <- S7::prop(x, p)
      if (S7::S7_inherits(v) && "label" %in% names(S7::props(v))) return(v@label)
      if (S7::S7_inherits(v) && "name" %in% names(S7::props(v))) return(v@name)
      if (is.atomic(v) && length(v) == 1L) return(format(v))
      paste0("<", class(v)[1], ">")
    }
    vals <- vapply(extra, describe, character(1))
    cat("  settings  : ", paste(extra, vals, sep = " = ", collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}
