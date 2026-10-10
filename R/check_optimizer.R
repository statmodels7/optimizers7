#' @include test_problems.R
#' @include generics.R
#' @include methods.R
NULL

#' @title Check What an Optimizer Reports
#'
#' @description
#' Puts an optimizer through twelve checks on what it *reports*, prints the
#' verdict on each, and then runs it over the standard test problems and
#' reports the gap that it reached on each one. It is intended for a caller
#' who adds a method. Every shipped method except [brent()] passes all
#' twelve checks; [brent()] minimizes a function of one variable and cannot
#' run the battery, whose problems have more parameters.
#'
#' @details
#' # Requirements and performance
#'
#' The twelve numbered checks are statements that an optimizer must satisfy
#' whatever its performance. The performance comes afterwards, as a table of
#' gaps without a verdict. The two are kept apart because a slow method is
#' not a defective one: gradient descent does not solve Rosenbrock in five
#' hundred iterations, and that is a documented property of the method.
#'
#' The only performance requirement among the numbered checks is check 12,
#' the minimization of a quadratic, which every correct method passes.
#'
#' # The twelve checks
#'
#' 1. `value` is the objective at `par`. A method that reports the value at a
#'    point that it has since left passes every test written in terms of the
#'    value alone.
#' 2. the reported gradient is the gradient at `par`. Checked only for an
#'    optimizer that declares `"gradient"` through [optimizer_provides()].
#'    [bundle()] reports an aggregate subgradient and does not declare a
#'    gradient, so the check does not apply to it.
#' 3. `converged` follows the stopping rule and is never inferred from the
#'    run having ended. Checked by starving the optimizer of iterations: a
#'    run cut off after one must not report success.
#' 4. budgets are respected: `iterations` never exceeds `maxit`.
#' 5. evaluations are counted: the reported number of evaluations is not
#'    zero.
#' 6. the trace, when kept, is a data frame whose key column (`iteration`, or
#'    `start` for a [multistart()] trace) starts at one and increases.
#' 7. bounds are respected **strictly**. A probability of exactly 1 is not a
#'    probability inside \eqn{(0, 1)}, and the next computation of the caller
#'    often divides by it.
#' 8. the run repeats. A deterministic method gives the same answer twice; a
#'    stochastic one gives it again from the seed it recorded, which tests
#'    the recording as well as the repeatability.
#' 9. [maximize()] is [minimize()] of the negative.
#' 10. a stopping rule that the optimizer cannot evaluate is rejected when
#'     the run starts, and is not accepted and left never to fire.
#' 11. a starting point where the objective is not finite raises an error,
#'     so that a run does not return `NaN` silently.
#' 12. it minimizes a quadratic.
#'
#' An optimizer that returns an invented point with `converged = TRUE` fails
#' six of the twelve checks: 1, 2, 3, 10, 11 and 12. The six that it passes
#' concern bookkeeping.
#'
#' # The problem battery
#'
#' The table reports the gap between the value reached and the known minimum,
#' as information. A large gap on `rastrigin` or `himmelblau` means the
#' method found a different local minimum, which for a local method is
#' correct behavior, and the `note` column labels those two. A large gap on
#' `abs_sum` means that the method stopped at a kink, the case for which
#' [bundle()] and the derivative-free methods are intended.
#'
#' @param optimizer The [optimizer()] to check. Anything else raises an error
#'   naming `bfgs()` as an example.
#' @param problems A list of problems in the shape [test_problems()] returns.
#'   Defaults to all eight. Pass a subset to keep the report short:
#'   `test_problems("sphere")`.
#' @param verbose `TRUE` (default) prints the twelve verdicts and the battery
#'   table; `FALSE` returns the same information silently.
#' @param tol Tolerance for the checks that compare numbers, a single
#'   positive number, default `1e-6`. It applies to check 1 and to the
#'   comparison of the two values in check 9. Check 2 uses a fixed `1e-3` on
#'   the gradient, check 9 compares the two points with a fixed relative
#'   tolerance of `1e-5`, check 12 uses a fixed `1e-3` on the distance to the
#'   solution, and check 7 tests the bounds strictly, without a tolerance.
#'
#' @return Invisibly, a named list of two:
#'   \describe{
#'     \item{`checks`}{a named logical vector of length 12, one entry per
#'       numbered check, named by what the check asserts.}
#'     \item{`battery`}{a data frame with one row per problem and the columns
#'       `problem`, `value`, `gap`, `converged`, `evaluations` and `note`. A
#'       problem that the optimizer could not run keeps its name in
#'       `problem`, has `NA` in the other numeric and logical columns, and has
#'       the error message in `note`.}
#'   }
#'
#' @examples
#' check_optimizer(bfgs())
#'
#' # For a method that computes no gradient, the gradient check does not
#' # apply, and the other checks are the same.
#' res <- check_optimizer(nelder_mead(), problems = test_problems("sphere"))
#' all(res$checks)
#'
#' # The battery is returned as a table.
#' check_optimizer(cg(), verbose = FALSE)$battery
#'
#' @seealso [test_problems()] for the battery, [check_criterion()] and
#'   [check_bounds()] for the two functions that a new method has to call,
#'   [optimizer_provides()] for the declaration read by check 2.
#' @export
check_optimizer <- function(optimizer, problems = test_problems(),
                            verbose = TRUE, tol = 1e-6) {
  if (!S7::S7_inherits(optimizer, optimizer_class())) {
    stop("'optimizer' must be an optimizer object, e.g. bfgs().", call. = FALSE)
  }
  sph <- test_problems("sphere")[[1]]
  ros <- test_problems("rosenbrock")[[1]]
  provides <- optimizer_provides(optimizer)

  run <- function(o, p, ...) {
    tryCatch(minimize(o, p$fn, p$par, gr = p$gr, he = p$he, ...),
             error = function(e) e)
  }
  failed <- function(r) inherits(r, "error")

  ok <- stats::setNames(rep(NA, 12), c(
    "value agrees with par", "gradient agrees with par",
    "convergence is not assumed", "budgets are respected",
    "evaluations are counted", "trace is well formed",
    "bounds are respected strictly", "the run repeats",
    "maximize mirrors minimize", "an unevaluable rule is rejected",
    "a bad starting point is an error", "it minimizes a quadratic"))

  base <- run(optimizer, sph)
  if (failed(base)) {
    stop("The optimizer could not run the simplest problem in the battery: ",
         conditionMessage(base), call. = FALSE)
  }

  # [1] and [12]
  ok[1] <- abs(base@value - sph$fn(base@par)) <= tol * (1 + abs(base@value))
  ok[12] <- max(abs(base@par - sph$solution)) <= 1e-3

  # [2] only for an optimizer that offers its gradient to a stopping rule, and
  # is therefore claiming the thing it reports IS the gradient at par.
  ok[2] <- if (!("gradient" %in% provides)) TRUE
           else !is.null(base@gradient) &&
                max(abs(base@gradient - sph$gr(base@par))) <= 1e-3

  # [3] Starved of iterations on a hard problem, a run must not claim success.
  starved <- run(with_maxit(optimizer, 1), ros)
  ok[3] <- !failed(starved) && !starved@converged &&
    !identical(starved@criterion_met, optimizer@criterion@label)

  # [4] Against the budget the object itself carries, not against the number
  # passed in. For an ordinary optimizer they are the same and this asserts
  # `<= 3`; for a wrapper they are not, since multistart()'s own maxit counts
  # STARTS while the budget being varied belongs to the optimizer inside.
  # Comparing the outer iteration count against the inner budget was the first
  # version of this check, and it failed a correct optimizer.
  capped_o <- with_maxit(optimizer, 3)
  capped <- run(capped_o, ros)
  ok[4] <- !failed(capped) && capped@iterations <= capped_o@maxit

  # [5]
  ok[5] <- is.numeric(base@counts) && base@counts[["f"]] > 0

  # [6] A trace is one row per STEP for a descent method and one row per START
  # for multistart(), so the key is `iteration` in the first case and `start`
  # in the second; whichever it is must begin at one and not go backwards.
  #
  # ⚠️ This read `traced@trace$iteration`, and `$` on a data frame PARTIALLY
  # MATCHES: on a multistart trace, whose columns are start/value/converged/
  # iterations, it silently returned the `iterations` column and compared the
  # first start's iteration COUNT against one. That passed for as long as the
  # inner method happened to solve the sphere in a single iteration and broke
  # the moment it took two, which is how it was found. `[[` does not partial
  # match, and asking for the key the trace actually has is the check.
  traced <- run(with_trace(optimizer), sph)
  ok[6] <- !failed(traced) &&
    (is.null(traced@trace) ||
       local({
         tr <- traced@trace
         if (!is.data.frame(tr) || !nrow(tr)) return(FALSE)
         key <- if ("iteration" %in% names(tr)) "iteration" else
                if ("start" %in% names(tr)) "start" else NA_character_
         if (is.na(key)) return(FALSE)
         v <- tr[[key]]
         is.numeric(v) && v[1] == 1 && !is.unsorted(v)
       }))

  # [7] A box whose ceiling binds, so the transform is pushed hard. Vacuous
  # for a method that declares it takes its constraint another way.
  far <- list(fn = function(p) sum((p - 50)^2), gr = function(p) 2 * (p - 50),
              he = function(p) diag(2, length(p)), par = c(0, 0, 0))
  ok[7] <- if (!optimizer_bounded(optimizer)) TRUE else {
    bnd <- run(optimizer, far, lower = -5, upper = 5)
    !failed(bnd) && all(bnd@par > -5) && all(bnd@par < 5)
  }

  # [8] Deterministic twice over, or stochastic and repeatable from its own
  # recorded seed -- which tests the recording at the same time.
  again <- if (is.null(base@seed)) {
    run(optimizer, sph)
  } else {
    assign(".Random.seed", base@seed, envir = globalenv())
    run(optimizer, sph)
  }
  ok[8] <- !failed(again) && isTRUE(all.equal(base@par, again@par))

  # [9] From the SAME random stream as the run it is compared against, for the
  # reason check [8] re-seeds: the property under test is that maximize mirrors
  # minimize, and for a stochastic method that statement is only meaningful
  # given the same draws. Without it the check silently required the optimizer
  # to converge tightly enough that two independent runs agree to 1e-5 --
  # which the mesh-shrinking methods happen to do and a global search does not,
  # so the check passed for the wrong reason everywhere and failed for the
  # wrong reason on simulated annealing.
  if (!is.null(base@seed)) {
    assign(".Random.seed", base@seed, envir = globalenv())
  }
  mx <- tryCatch(maximize(optimizer, function(p) -sph$fn(p), sph$par,
                          gr = function(p) -sph$gr(p),
                          he = function(p) -sph$he(p)),
                 error = function(e) e)
  ok[9] <- !failed(mx) &&
    isTRUE(all.equal(mx@par, base@par, tolerance = 1e-5)) &&
    abs(mx@value + base@value) <= tol * (1 + abs(base@value))

  # [10] Every state component the optimizer does NOT provide must have its
  # criterion refused. An optimizer providing everything passes vacuously.
  rules <- list(gradient = crit_grad(), stationarity = crit_stationary())
  unmet <- setdiff(names(rules), provides)
  ok[10] <- all(vapply(unmet, function(nm) {
    r <- run(with_criterion(optimizer, rules[[nm]]), sph)
    failed(r) && grepl("does not provide", conditionMessage(r))
  }, logical(1)))

  # [11]
  bad <- tryCatch(minimize(optimizer, function(p) NaN, sph$par),
                  error = function(e) e)
  ok[11] <- failed(bad)

  battery <- run_battery(optimizer, problems)

  if (verbose) print_optimizer_check(optimizer, ok, battery)
  invisible(list(checks = ok, battery = battery))
}


#' Run an Optimizer Over the Battery
#'
#' @description
#' Runs the optimizer on each problem from that problem's own starting point
#' and records the gap between the value reached and the known minimum. The
#' result is information: no row is a pass or a failure.
#'
#' @details
#' An error on a problem is caught and reported in the `note` column, so
#' that it does not discard the results on the other problems. Such a row
#' keeps the name of the problem and has `NA` in the other numeric and
#' logical columns.
#'
#' @param optimizer The [optimizer()] to run.
#' @param problems A list in the shape [test_problems()] returns.
#'
#' @return A data frame with one row per problem and the columns `problem`
#'   (character), `value` and `gap` (numeric), `converged` (logical),
#'   `evaluations` (integer) and `note` (character, `"multimodal"`,
#'   `"non-smooth"`, an error message, or empty).
#'
#' @keywords internal
run_battery <- function(optimizer, problems) {
  rows <- lapply(problems, function(p) {
    r <- tryCatch(minimize(optimizer, p$fn, p$par, gr = p$gr, he = p$he),
                  error = function(e) e)
    if (inherits(r, "error")) {
      data.frame(problem = p$name, value = NA_real_, gap = NA_real_,
                 converged = NA, evaluations = NA_integer_,
                 note = conditionMessage(r), stringsAsFactors = FALSE)
    } else {
      data.frame(problem = p$name, value = r@value, gap = r@value - p$value,
                 converged = r@converged,
                 evaluations = as.integer(r@counts[["f"]]),
                 note = if (isTRUE(p$multimodal)) "multimodal"
                        else if (isFALSE(p$smooth)) "non-smooth" else "",
                 stringsAsFactors = FALSE)
    }
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Print the Report of check_optimizer
#'
#' @description
#' Writes the twelve verdicts, one per line as `[PASSED]` or `[FAILED]`, then
#' a summary naming every failing check, then the battery as one line per
#' problem with its gap, its convergence flag, its evaluation count and its
#' note.
#'
#' @param optimizer The [optimizer()] checked, read for its name.
#' @param ok The named logical vector of twelve checks.
#' @param battery The data frame [run_battery()] returned.
#'
#' @return Invisibly `NULL`. Called for the output.
#'
#' @keywords internal
print_optimizer_check <- function(optimizer, ok, battery) {
  cat("Checking optimizer: ", optimizer@name, "\n", sep = "")
  w <- max(nchar(names(ok)))
  for (i in seq_along(ok)) {
    verdict <- if (isTRUE(ok[i])) "[PASSED]" else "[FAILED]"
    cat(sprintf("  [%2d] %-*s %s\n", i, w + 1, paste0(names(ok)[i], ":"),
                verdict))
  }
  bad <- names(ok)[!vapply(ok, isTRUE, logical(1))]
  cat("\n  ", if (!length(bad)) "All checks passed."
              else paste0(length(bad), " check(s) FAILED: ",
                          paste(bad, collapse = "; ")), "\n", sep = "")

  cat("\n  battery (gap from the known minimum; information, not a verdict)\n")
  for (i in seq_len(nrow(battery))) {
    cat(sprintf("    %-12s gap %9.2e  %-5s  %6d evals  %s\n",
                battery$problem[i], battery$gap[i],
                if (isTRUE(battery$converged[i])) "conv" else "-",
                if (is.na(battery$evaluations[i])) 0L else battery$evaluations[i],
                battery$note[i]))
  }
  invisible(NULL)
}


# --- rebuilding an optimizer with one setting changed -----------------------
#
# check_optimizer has to vary maxit, the trace and the criterion on whatever
# optimizer it was handed, without knowing its class. S7::set_props does that
# for an ordinary one; a wrapper has to pass the change inwards as well, which
# is what the MultiStart methods below are for.

#' Rebuild an Optimizer With a Different Iteration Budget
#'
#' @description
#' Returns a copy of the optimizer with `maxit` replaced, keeping its class
#' and every other setting. [check_optimizer()] uses it for check 3, which
#' gives an optimizer a single iteration and tests that it does not report
#' convergence, and which receives an optimizer of any class.
#'
#' @details
#' The default method sets the property with `S7::set_props()`, which is
#' correct for an optimizer whose own budget governs the run. A wrapper needs
#' a method of its own so that the change reaches the optimizer inside;
#' [multistart()] and [chain()] have one.
#'
#' @param optimizer The [optimizer()] to copy.
#' @param maxit The new budget, a single finite number at least 1.
#'
#' @return An optimizer of the same class as `optimizer`.
#'
#' @aliases with_maxit.optimizer
#' @keywords internal
with_maxit <- S7::new_generic("with_maxit", "optimizer",
                              function(optimizer, maxit) S7::S7_dispatch())

S7::method(with_maxit, optimizer) <- function(optimizer, maxit)
  S7::set_props(optimizer, maxit = maxit)


#' Rebuild an Optimizer With the Trace Switched On
#'
#' @description
#' Returns a copy of the optimizer with `keep_trace = TRUE`, keeping its
#' class and every other setting. [check_optimizer()] uses it for check 6,
#' which tests that the trace is well formed and therefore needs a run that
#' records one.
#'
#' @details
#' The default method sets the property with `S7::set_props()`. The method
#' for [chain()] also switches on the trace of every stage, because the trace
#' of a chain is assembled from the traces of its stages. A [multistart()]
#' needs no method of its own: its trace is the table of starts, which the
#' wrapper builds from the results of the inner runs.
#'
#' @param optimizer The [optimizer()] to copy.
#'
#' @return An optimizer of the same class as `optimizer`, with
#'   `keep_trace` `TRUE`.
#'
#' @aliases with_trace.optimizer
#' @keywords internal
with_trace <- S7::new_generic("with_trace", "optimizer",
                              function(optimizer) S7::S7_dispatch())

S7::method(with_trace, optimizer) <- function(optimizer)
  S7::set_props(optimizer, keep_trace = TRUE)


#' Rebuild an Optimizer With a Different Stopping Rule
#'
#' @description
#' Returns a copy of the optimizer with its criterion replaced, keeping its
#' class and every other setting. For a wrapper it replaces the rule that is
#' actually consulted, which is not always the one on the outside.
#'
#' @details
#' [multistart()] carries a criterion property that is used for printing;
#' the rule that the run evaluates belongs to the optimizer inside. Setting
#' the outer one alone changes only the printing, and a check that relied on
#' it would test the default rule:
#'
#' ```r
#' ms <- multistart(bfgs(), n = 3)
#' with_criterion(ms, crit_abs_obj(1e-4))@optimizer@criterion@label
#' # "|df| < 1e-04"  (the inner rule changed too)
#'
#' S7::set_props(ms, criterion = crit_abs_obj(1e-4))@optimizer@criterion@label
#' # "gradient (max-norm) < 1e-06 or ..."  (unchanged)
#' ```
#'
#' [chain()] has a method of its own for the same reason, since its rule is
#' the rule of the last stage.
#'
#' @param optimizer The [optimizer()] to copy.
#' @param criterion The new rule, a [criterion()] object.
#'
#' @return An optimizer of the same class as `optimizer`.
#'
#' @aliases with_criterion.optimizer
#' @keywords internal
with_criterion <- S7::new_generic("with_criterion", "optimizer",
                                  function(optimizer, criterion)
                                    S7::S7_dispatch())

S7::method(with_criterion, optimizer) <- function(optimizer, criterion)
  S7::set_props(optimizer, criterion = criterion)
