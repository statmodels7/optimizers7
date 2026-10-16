#' @include criterion.R
NULL

#' @title S7 Class for Optimization Algorithms
#'
#' @description
#' This class is the abstract parent of every algorithm in the package. An
#' optimizer object carries the algorithm and every setting that the algorithm
#' uses: its stopping rule, its two budgets, what it reports during the run and
#' what it keeps. It is built once and can be reused, printed, stored beside a
#' result and passed around, so that all the settings of a run are stored in
#' the object and not in the call that started it.
#'
#' @details
#' # The seven properties of every algorithm
#'
#' Each algorithm is a subclass that adds its own settings: a Hessian repair
#' for [newton()], a curvature tolerance for [bfgs()], the two moment decay
#' rates for [adam()]. The seven properties below are common to all of them,
#' so a caller can replace one optimizer with another and leave the rest of
#' the call unchanged. `print()` shows the stopping rule and the two budgets,
#' and then the settings of the algorithm on a `settings` line; `verbose`,
#' `refresh` and `keep_trace` are not shown.
#'
#' # The iteration and evaluation budgets
#'
#' `maxit` bounds the iterations and `max_eval` bounds the evaluations of the
#' objective, and the two differ whenever a line search is expensive: one
#' iteration may spend many evaluations, and the number of backtracks does
#' not count towards `maxit`. A run stopped by `max_eval` reports
#' `converged = FALSE`, with `criterion_met` equal to
#' `evaluation budget reached` (see the examples).
#'
#' `maxit` must be finite, because it is the stop of last resort; `max_eval`
#' may be `Inf`, which is its default.
#'
#' @param name A short character name, shown when the object is printed and
#'   carried into the result. Set by each constructor; a caller does not
#'   supply it.
#' @param criterion The stopping rule, a [criterion()] object. Anything else
#'   raises an error naming [crit_grad()] as an example. An algorithm also
#'   rejects a rule that it cannot evaluate, through [check_criterion()].
#' @param maxit Maximum iterations, a single number at least 1 and finite.
#'   `Inf` is rejected, because a run without a last resort can run forever
#'   on an objective that never satisfies its rule.
#' @param max_eval Maximum evaluations of the objective, a single number at
#'   least 1. `Inf` is the default and is allowed, so the budget is off
#'   unless a finite value is given. A run stopped here reports
#'   `converged = FALSE`, because the budget ended it and not the rule.
#' @param verbose `TRUE` or `FALSE`; whether to report progress during the
#'   run. `NA` is rejected.
#' @param refresh Report every `refresh` iterations when `verbose` is `TRUE`.
#'   It is a single non-negative number; `0` is allowed and reports only the
#'   final summary.
#' @param keep_trace `TRUE` or `FALSE`; whether to store the iteration path in
#'   the `trace` field of the result. It is off by default, since the path
#'   costs one row per iteration.
#'
#' @return An S7 object of class `optimizer`, with properties `name`
#'   (character), `criterion` (a [criterion()]), `maxit`, `max_eval` and
#'   `refresh` (numeric) and `verbose` and `keep_trace` (logical). The class
#'   is abstract, so it is never returned directly: every value is an object
#'   of one of its subclasses.
#'
#' @examples
#' # The class is abstract, so it cannot be built: use one of the constructors.
#' try(optimizer(name = "mine"))
#'
#' # Every algorithm has the same seven properties, and each one has properties
#' # of its own.
#' shared <- c("name", "criterion", "maxit", "max_eval", "verbose",
#'             "refresh", "keep_trace")
#' stopifnot(all(shared %in% names(S7::props(nelder_mead()))))
#' setdiff(names(S7::props(newton())), shared)
#' setdiff(names(S7::props(adam())), shared)
#'
#' # Printing shows the rule and the budgets before the settings of the
#' # algorithm.
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
#'   [check_optimizer()] for the checks that a new subclass has to pass.
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


#' Validate the Settings That Every Optimizer Shares
#'
#' @description
#' Checks the six arguments common to the constructors (the stopping rule,
#' `maxit`, `max_eval`, `verbose`, `refresh` and `keep_trace`), so that every
#' constructor that calls it rejects invalid settings with the same messages.
#' It is called at the top of each such constructor, before any validation
#' specific to the algorithm.
#'
#' @details
#' The rules are:
#' - `criterion` must inherit from [criterion()].
#' - `maxit` must be a single number at least 1 and finite.
#' - `max_eval` must be a single number at least 1; `Inf` passes.
#' - `refresh` must be a single number at least 0.
#' - `verbose` and `keep_trace` must each be `TRUE` or `FALSE`; `NA` is
#'   rejected.
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
#' that make this test take an argument named `criterion`, which shadows the
#' class inside their own body.
#'
#' @details
#' Written directly, `S7_inherits(criterion, criterion)` inside a function
#' whose formal is `criterion` passes the value of the argument as the class,
#' and S7 stops with ``` `class` must be an <S7_class> or NULL ```. A
#' function of no arguments looks the class up in the namespace, where the
#' formal does not shadow it. [check_optimizer_args()] is the caller.
#'
#' @return The [criterion()] class object, an `S7_class`.
#'
#' @keywords internal
criterion_class <- function() criterion


#' The optimizer Class Object
#'
#' @description
#' Returns the [optimizer()] class itself, for the same reason as
#' [criterion_class()]: the functions that accept an arbitrary
#' optimizer and check their argument all take an argument named
#' `optimizer`, which shadows the class inside their body.
#'
#' @details
#' The callers are [multistart()], which wraps an optimizer, [chain()],
#' which holds a list of them, and [check_optimizer()], which tests one.
#' Each rejects a value that is not an optimizer with a message of its own;
#' without the check the failure would appear later, in the dispatch of
#' [minimize()].
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
#' rule, its two budgets, and the settings specific to the algorithm.
#' `verbose`, `refresh` and `keep_trace` are not printed.
#'
#' @details
#' The settings line describes each value according to its type. A number
#' prints as itself; an object with a `label`, such as a [criterion()] or a
#' [line_search()], prints its label; an object with a `name`, such as a
#' nested optimizer, prints its name; any other value prints as its class in
#' angle brackets. A setting of a new kind therefore prints a readable value
#' and does not stop the method.
#'
#' @param x An [optimizer()] object.
#' @param ... Unused.
#'
#' @return `x`, invisibly. Called for the output.
#'
#' @examples
#' # The line search is an object and shows its own label.
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
