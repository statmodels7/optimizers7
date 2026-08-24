#' @include criterion.R
NULL

#' @title S7 Class for Optimization Algorithms
#'
#' @description
#' The abstract parent of every algorithm in the package. An optimizer object
#' carries the algorithm and every setting the algorithm obeys: its stopping
#' rule, its two budgets, what it reports as it goes and what it keeps. It is
#' built once and can be reused, printed, stored beside a result and passed
#' around, so nothing about a run is hidden inside the call that started it.
#'
#' @details
#' # The seven properties every algorithm has
#'
#' Each algorithm is a subclass adding the settings that are its own:
#' [newton()] a Hessian repair, [bfgs()] a curvature tolerance, [adam()] its
#' three decay rates. The seven below are common to all of them, so a caller
#' can replace one optimizer with another and leave every other line alone.
#' `print()` shows the shared settings first and the algorithm's own on a
#' `settings` line.
#'
#' # The two budgets are different budgets
#'
#' `maxit` bounds progress and `max_eval` bounds work, and they diverge
#' whenever a line search is expensive: one iteration may spend many
#' evaluations, and a search that backtracks thirty times is invisible to
#' `maxit`. On Rosenbrock from the customary start, [bfgs()] with no cap
#' takes 35 iterations and 210 evaluations; capped at `max_eval = 40` it
#' stops after 7 iterations and 42 evaluations, reporting `converged = FALSE`
#' and `evaluation budget exhausted`.
#'
#' `maxit` must be finite, because it is the stop of last resort; `max_eval`
#' may be `Inf`, which is its default.
#'
#' @param name A short character name, shown when the object is printed and
#'   carried into the result. Set by each constructor; a caller has no reason
#'   to supply it.
#' @param criterion The stopping rule, a [criterion()] object. Anything else
#'   raises an error naming [crit_grad()] as an example. An algorithm also
#'   rejects a rule it cannot evaluate, through [check_criterion()].
#' @param maxit Maximum iterations, a single number at least 1 and **finite**.
#'   `Inf` is refused: a run with no last resort can hang on an objective
#'   that never satisfies its rule.
#' @param max_eval Maximum evaluations of the objective, a single number at
#'   least 1. `Inf` is the default and is allowed, so the budget is off
#'   unless it is asked for. A run stopped here reports `converged = FALSE`,
#'   the budget having ended it rather than the rule.
#' @param verbose `TRUE` or `FALSE`; whether to report progress as the run
#'   goes. `NA` is refused.
#' @param refresh Report every `refresh` iterations when `verbose` is `TRUE`.
#'   A single non-negative number; `0` is allowed and reports only the final
#'   summary.
#' @param keep_trace `TRUE` or `FALSE`; whether to store the iteration path in
#'   the result's `trace` field. Off by default, the path costing one row per
#'   iteration.
#'
#' @return An S7 object of class `optimizer`, with properties `name`
#'   (character), `criterion` (a [criterion()]), `maxit`, `max_eval` and
#'   `refresh` (numeric) and `verbose` and `keep_trace` (logical). The class
#'   is abstract, so it is never returned directly: every value is an object
#'   of one of its subclasses.
#'
#' @examples
#' # Abstract, so it cannot be built: use one of the constructors.
#' try(optimizer(name = "mine"))
#'
#' # The same seven properties on every algorithm, plus each one's own.
#' shared <- c("name", "criterion", "maxit", "max_eval", "verbose",
#'             "refresh", "keep_trace")
#' stopifnot(all(shared %in% names(S7::props(nelder_mead()))))
#' setdiff(names(S7::props(newton())), shared)
#' setdiff(names(S7::props(adam())), shared)
#'
#' # Printing shows the rule and the budgets before anything algorithm-specific.
#' bfgs()
#'
#' # The two budgets bound different things.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' a <- minimize(bfgs(maxit = 1000), f, c(-1.2, 1))
#' b <- minimize(bfgs(maxit = 1000, max_eval = 40), f, c(-1.2, 1))
#' c(iterations = a@iterations, evaluations = a@counts[["f"]])
#' c(iterations = b@iterations, evaluations = b@counts[["f"]])
#' b@message
#'
#' @seealso [minimize()] for the run, [criterion()] for the stopping rule,
#'   [check_optimizer()] for what a subclass of your own has to promise.
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


#' Validate the Settings Every Optimizer Shares
#'
#' @description
#' Checks the seven arguments common to every constructor, so that all
#' fifteen reject the same nonsense in the same words. Called at the top of
#' each constructor, before any algorithm-specific validation.
#'
#' @details
#' The rules, and they differ from one another:
#' - `criterion` must inherit from [criterion()].
#' - `maxit` must be a single number at least 1 **and finite**.
#' - `max_eval` must be a single number at least 1; `Inf` passes.
#' - `refresh` must be a single number at least 0.
#' - `verbose` and `keep_trace` must each be `TRUE` or `FALSE`, `NA` refused.
#'
#' @param criterion The stopping rule.
#' @param maxit,max_eval,refresh Numeric budgets.
#' @param verbose,keep_trace Logical flags.
#'
#' @return Invisibly `TRUE`. Raises an error naming the offending argument
#'   otherwise.
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
  if (!is.finite(maxit)) {
    stop("'maxit' must be finite: it is the stop of last resort.",
         call. = FALSE)
  }
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

#' The criterion Class Object
#'
#' @description
#' Returns the [criterion()] class itself, for code that has to test whether
#' a value inherits from it. The indirection exists because the functions
#' asking that question take an argument named `criterion`, which shadows the
#' class inside their own body.
#'
#' @details
#' Written directly, `S7_inherits(criterion, criterion)` inside a function
#' whose formal is `criterion` compares the value with itself, and S7 stops
#' with ``` `class` must be an <S7_class> or NULL ```. Reaching the class
#' through a function of no arguments looks it up in the namespace, past the
#' formal. [check_optimizer_args()] is the caller.
#'
#' @return The [criterion()] class object, an `S7_class`.
#'
#' @keywords internal
criterion_class <- function() criterion


#' The optimizer Class Object
#'
#' @description
#' Returns the [optimizer()] class itself, for the same reason
#' [criterion_class()] exists: the three functions that accept an arbitrary
#' optimizer and must check what they were given all take an argument named
#' `optimizer`, which shadows the class inside their body.
#'
#' @details
#' The callers are [multistart()], which wraps an optimizer, [chain()],
#' which holds a list of them, and [check_optimizer()], which tests one.
#' Each rejects a non-optimizer with a message of its own; without the check
#' the failure would surface later, inside [minimize()] dispatch.
#'
#' @return The [optimizer()] class object, an `S7_class`.
#'
#' @keywords internal
optimizer_class <- function() optimizer


#' @title Print Method for Optimizers
#' @name print.optimizer
#'
#' @description
#' Shows an optimizer in four lines: its name, the label of its stopping
#' rule, its two budgets, and its algorithm-specific settings. Everything the
#' run will obey is on the object, so printing it is how a caller checks what
#' a fit was configured to do.
#'
#' @details
#' The settings line describes each value by what it is. A number prints as
#' itself; an object carrying a `label`, such as a [criterion()], prints its
#' label; an object carrying a `name`, such as a [line_search()] or a nested
#' optimizer, prints its name; anything else prints as its class in angle
#' brackets. A setting of a kind the package does not yet have therefore
#' prints something sensible instead of stopping the method.
#'
#' @param x An [optimizer()] object.
#' @param ... Unused.
#'
#' @return `x`, invisibly. Called for the output.
#'
#' @examples
#' # The line search is an object and shows its own name.
#' bfgs()
#'
#' # A derivative-free method has a different rule and different settings.
#' nelder_mead()
#'
#' # A wrapper names the optimizer it wraps.
#' multistart(bfgs(), n = 4)
#'
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
