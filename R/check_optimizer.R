#' @include test_problems.R
#' @include generics.R
#' @include methods.R
NULL

#' @title Check That an Optimizer Keeps Its Promises
#'
#' @description
#' Puts an optimizer through twelve checks on what it *reports*, prints the
#' verdict on each, and then runs it over the standard test problems and
#' reports the gap it reached on every one. Written for whoever adds a method
#' of their own; all twelve shipped methods pass all twelve checks.
#'
#' @details
#' # The contract and the power are different questions
#'
#' The twelve numbered checks are statements an optimizer must satisfy
#' however weak it is. How strong it is comes afterwards, as a table of gaps
#' with no verdict attached. Conflating the two would make the function
#' useless: gradient descent does not solve Rosenbrock in five hundred
#' iterations, and it is slow rather than broken, which is a documented
#' property of the method.
#'
#' The one performance requirement among the numbered checks is check 12,
#' minimizing a quadratic. That is a floor no correct method can fail.
#'
#' # The twelve checks
#'
#' 1. `value` is the objective at `par`. A method reporting a value from a
#'    point it has since left is a defect that survives every test written in
#'    terms of the value alone.
#' 2. the reported gradient is the gradient at `par`. Checked only for an
#'    optimizer that offers `"gradient"` through [optimizer_provides()], and
#'    so is claiming that what it reports is one. [bundle()] reports an
#'    aggregate subgradient and makes no such claim, so it is exempt.
#' 3. `converged` follows the stopping rule and is never inferred from the
#'    run having ended. Checked by starving the optimizer of iterations: a
#'    run cut off after one must not report success.
#' 4. budgets are respected: `iterations` never exceeds `maxit`.
#' 5. evaluations are counted. A method reporting zero of them evaluated
#'    nothing.
#' 6. the trace, when kept, is a data frame whose iteration numbers start at
#'    one and increase.
#' 7. bounds are respected **strictly**. A probability of exactly 1 is not a
#'    probability inside \eqn{(0, 1)}, and the caller's next act is usually
#'    to divide by it.
#' 8. the run repeats. A deterministic method gives the same answer twice; a
#'    stochastic one gives it again from the seed it recorded, which tests
#'    the recording as well as the repeatability.
#' 9. [maximize()] is [minimize()] of the negative.
#' 10. a stopping rule the optimizer cannot evaluate is rejected at
#'     construction, not accepted and left never to fire.
#' 11. a starting point where the objective is not finite raises an error, so
#'     no run quietly returns `NaN`.
#' 12. it minimizes a quadratic.
#'
#' A deliberately lying optimizer, one returning a made-up point with
#' `converged = TRUE`, fails six of the twelve: 1, 2, 3, 10, 11 and 12. The
#' six it passes are the bookkeeping ones, which is why the battery is not
#' the whole check.
#'
#' # The problem battery
#'
#' The table reports the gap between the value reached and the known minimum,
#' as information. A large gap on `rastrigin` or `himmelblau` means the
#' method found a different local minimum, which for a local method is
#' correct behavior, and the `note` column labels those two. A large gap on
#' `abs_sum` means the method was defeated by a kink, the case [bundle()] and
#' the derivative-free methods exist for.
#'
#' @param optimizer The [optimizer()] to check. Anything else raises an error
#'   naming `bfgs()` as an example.
#' @param problems A list of problems in the shape [test_problems()] returns.
#'   Defaults to all eight. Pass a subset to keep the report short:
#'   `test_problems("sphere")`.
#' @param verbose `TRUE` (default) prints the twelve verdicts and the battery
#'   table; `FALSE` returns the same information silently.
#' @param tol Tolerance for the checks that compare numbers, a single
#'   positive number, default `1e-6`. It governs checks 1, 2 and 9; check 12
#'   uses a fixed `1e-3` on the distance to the solution, and check 7 tests
#'   the bounds strictly and reads no tolerance at all.
#'
#' @return Invisibly, a named list of two:
#'   \describe{
#'     \item{`checks`}{a named logical vector of length 12, one entry per
#'       numbered check, named by what the check asserts.}
#'     \item{`battery`}{a data frame with one row per problem and the columns
#'       `problem`, `value`, `gap`, `converged`, `evaluations` and `note`. A
#'       problem the optimizer could not run at all has `NA` throughout and
#'       the error message in `note`.}
#'   }
#'
#' @examples
#' check_optimizer(bfgs())
#'
#' # A method that computes no gradient is held to fewer claims, and to the
#' # same standard on the ones it does make.
#' res <- check_optimizer(nelder_mead(), problems = test_problems("sphere"))
#' all(res$checks)
#'
#' # The battery is a table, so it can be read rather than printed.
#' check_optimizer(cg(), verbose = FALSE)$battery
#'
#' @seealso [test_problems()] for the battery, [check_criterion()] and
#'   [check_bounds()] for the two pieces a method of your own has to use,
#'   [optimizer_provides()] for the claim check 2 rests on.
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
    tryCatch(minimize(o, p$fn, p$par, gr = p$gr, ...),
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
              par = c(0, 0, 0))
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
                          gr = function(p) -sph$gr(p)),
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
#' A problem the optimizer cannot run at all is caught rather than
#' propagated, so one method that refuses one problem does not lose the other
#' seven. Such a row carries `NA` in every numeric column and the error
#' message in `note`.
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
    r <- tryCatch(minimize(optimizer, p$fn, p$par, gr = p$gr),
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
#' and every other setting. [check_optimizer()] needs it for check 3, which
#' starves an optimizer of iterations to see whether it still claims
#' convergence, and it cannot name the class it was handed.
#'
#' @details
#' The default method sets the property with `S7::set_props()`, which is
#' right for any optimizer whose own budget is the one the run obeys. A
#' wrapper needs a method of its own so that the change reaches the optimizer
#' inside; [multistart()] has one.
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
#' class and every other setting. [check_optimizer()] needs it for check 6,
#' which asks whether the trace is well formed, and a trace has to be asked
#' for.
#'
#' @details
#' The default method sets the property with `S7::set_props()`. A wrapper
#' needs a method of its own so that the inner optimizer records a path too;
#' [multistart()] has one.
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
#' [multistart()] carries a criterion so that printing it tells the truth;
#' the rule the run evaluates belongs to the optimizer inside. Setting the
#' outer one alone changes the printing and nothing else, which is exactly
#' the sort of thing that makes a check pass while testing nothing:
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
#' [chain()] has a method of its own for the same reason, its reported rule
#' being the last stage's.
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
