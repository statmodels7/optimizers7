#' @include test_problems.R
#' @include generics.R
#' @include methods.R
NULL

#' @title Check That an Optimiser Keeps Its Promises
#'
#' @description
#' Runs an optimiser through a series of checks on what it \emph{reports}, and
#' then through the standard test problems. Written for whoever adds a method of
#' their own, and run against every method here.
#'
#' @param optimizer The \code{\link{optimizer}} to check.
#' @param problems The battery; defaults to \code{\link{test_problems}()}.
#' @param verbose Print the report? Defaults to \code{TRUE}.
#' @param tol Tolerance for the checks that compare numbers. Defaults to
#'   \code{1e-6}.
#'
#' @details
#' \strong{What this checks is the contract, not the power.} Those are different
#' questions and conflating them would make the function useless: gradient
#' descent does not solve Rosenbrock in five hundred iterations, and it is not
#' broken — it is slow, which is a documented property of the method and not a
#' defect for a validator to report. So the numbered checks are all statements an
#' optimiser must satisfy however weak it is, and how strong it is comes
#' afterwards, as a table of gaps rather than as a verdict.
#'
#' The one performance requirement among the numbered checks is that the
#' optimiser minimises a quadratic. That is a floor no correct method can fail.
#'
#' \subsection{The checks}{
#' \enumerate{
#'   \item \code{value} is the objective at \code{par}. A method that reports a
#'     value from a point it has since left is the kind of defect that survives
#'     every test written in terms of the value alone.
#'   \item the reported gradient is the gradient at \code{par}, checked only for
#'     optimisers that offer \code{"gradient"} to a stopping rule and so are
#'     claiming it is one. \code{\link{bundle}} reports an aggregate subgradient
#'     and does not make that claim, so it is not held to it.
#'   \item \code{converged} follows the stopping rule and is never inferred from
#'     the run having ended. Checked by starving the optimiser of iterations: a
#'     run cut off after one must not report success.
#'   \item budgets are respected — \code{iterations} never exceeds \code{maxit}.
#'   \item evaluations are counted; a method reporting zero of them did not
#'     evaluate anything.
#'   \item the trace, when kept, is a data frame whose iteration numbers run
#'     from one and increase.
#'   \item bounds are respected \strong{strictly}: a probability of exactly 1 is
#'     not a probability inside \eqn{(0, 1)}, and the caller's next act is
#'     usually to divide by it.
#'   \item the run repeats. A deterministic method must give the same answer
#'     twice; a stochastic one must give it again from the seed it recorded,
#'     which tests the recording as well as the repeatability.
#'   \item \code{\link{maximize}} is \code{\link{minimize}} of the negative.
#'   \item a stopping rule the optimiser cannot evaluate is refused, rather than
#'     accepted and left never to fire.
#'   \item a starting point where the objective is not finite is an error, not a
#'     run that quietly returns \code{NaN}.
#'   \item it minimises a quadratic.
#' }
#' }
#'
#' \subsection{Reading the battery}{
#' The table reports the gap between the value reached and the known minimum,
#' and it is information rather than judgement. A large gap on
#' \code{rastrigin} or \code{himmelblau} means the method found a different
#' local minimum, which for a local method is correct behaviour; a large gap on
#' \code{abs_sum} means it was defeated by a kink, which is what
#' \code{\link{bundle}} and the derivative-free methods are for.
#' }
#'
#' @return Invisibly, a named list: \code{checks}, a logical vector with one
#'   entry per numbered check, and \code{battery}, a data frame of gaps.
#'
#' @examples
#' check_optimizer(bfgs())
#'
#' # a method that does not compute a gradient is held to fewer claims, and to
#' # the same standard on the ones it does make
#' check_optimizer(nelder_mead(), problems = test_problems("sphere"))
#'
#' @seealso \code{\link{test_problems}}, \code{\link{minimize}}
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
    "maximize mirrors minimize", "an unevaluable rule is refused",
    "a bad starting point is an error", "it minimises a quadratic"))

  base <- run(optimizer, sph)
  if (failed(base)) {
    stop("The optimiser could not run the simplest problem in the battery: ",
         conditionMessage(base), call. = FALSE)
  }

  # [1] and [12]
  ok[1] <- abs(base@value - sph$fn(base@par)) <= tol * (1 + abs(base@value))
  ok[12] <- max(abs(base@par - sph$solution)) <= 1e-3

  # [2] only for an optimiser that offers its gradient to a stopping rule, and
  # is therefore claiming the thing it reports IS the gradient at par.
  ok[2] <- if (!("gradient" %in% provides)) TRUE
           else !is.null(base@gradient) &&
                max(abs(base@gradient - sph$gr(base@par))) <= 1e-3

  # [3] Starved of iterations on a hard problem, a run must not claim success.
  starved <- run(with_maxit(optimizer, 1), ros)
  ok[3] <- !failed(starved) && !starved@converged &&
    !identical(starved@criterion_met, optimizer@criterion@label)

  # [4] Against the budget the object itself carries, not against the number
  # passed in. For an ordinary optimiser they are the same and this asserts
  # `<= 3`; for a wrapper they are not, since multistart()'s own maxit counts
  # STARTS while the budget being varied belongs to the optimiser inside.
  # Comparing the outer iteration count against the inner budget was the first
  # version of this check, and it failed a correct optimiser.
  capped_o <- with_maxit(optimizer, 3)
  capped <- run(capped_o, ros)
  ok[4] <- !failed(capped) && capped@iterations <= capped_o@maxit

  # [5]
  ok[5] <- is.numeric(base@counts) && base@counts[["f"]] > 0

  # [6]
  traced <- run(with_trace(optimizer), sph)
  ok[6] <- !failed(traced) &&
    (is.null(traced@trace) ||
       (is.data.frame(traced@trace) && nrow(traced@trace) > 0 &&
          !is.null(traced@trace$iteration) &&
          traced@trace$iteration[1] == 1 &&
          !is.unsorted(traced@trace$iteration)))

  # [7] A box whose ceiling binds, so the transform is pushed hard.
  bx <- list(c(-5, 5), c(-5, 5), c(-5, 5))
  far <- list(fn = function(p) sum((p - 50)^2), gr = function(p) 2 * (p - 50),
              par = c(0, 0, 0))
  bnd <- run(optimizer, far, bounds = bx)
  ok[7] <- !failed(bnd) && all(bnd@par > -5) && all(bnd@par < 5)

  # [8] Deterministic twice over, or stochastic and repeatable from its own
  # recorded seed -- which tests the recording at the same time.
  again <- if (is.null(base@seed)) {
    run(optimizer, sph)
  } else {
    assign(".Random.seed", base@seed, envir = globalenv())
    run(optimizer, sph)
  }
  ok[8] <- !failed(again) && isTRUE(all.equal(base@par, again@par))

  # [9]
  mx <- tryCatch(maximize(optimizer, function(p) -sph$fn(p), sph$par,
                          gr = function(p) -sph$gr(p)),
                 error = function(e) e)
  ok[9] <- !failed(mx) &&
    isTRUE(all.equal(mx@par, base@par, tolerance = 1e-5)) &&
    abs(mx@value + base@value) <= tol * (1 + abs(base@value))

  # [10] Every state component the optimiser does NOT provide must have its
  # criterion refused. An optimiser providing everything passes vacuously.
  rules <- list(gradient = crit_grad(), objective = crit_abs_obj(),
                stationarity = crit_stationary())
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


#' Run an Optimiser Over the Battery
#'
#' @description
#' The gap from each known minimum, as information rather than as a verdict.
#'
#' @param optimizer The \code{\link{optimizer}}.
#' @param problems A list in the shape \code{\link{test_problems}} returns.
#'
#' @return A data frame with one row per problem.
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
#' @param optimizer The optimiser checked.
#' @param ok The logical vector of checks.
#' @param battery The data frame of gaps.
#'
#' @return Invisibly \code{NULL}.
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


# --- rebuilding an optimiser with one setting changed -----------------------
#
# check_optimizer has to vary maxit, the trace and the criterion on whatever
# optimiser it was handed, without knowing its class. S7::set_props does that
# for an ordinary one; a wrapper has to pass the change inwards as well, which
# is what the MultiStart methods below are for.

#' Rebuild an Optimiser With a Different Iteration Budget
#' @param optimizer The \code{\link{optimizer}}.
#' @param maxit The new budget.
#' @return An optimiser of the same class.
#' @keywords internal
with_maxit <- S7::new_generic("with_maxit", "optimizer",
                              function(optimizer, maxit) S7::S7_dispatch())

S7::method(with_maxit, optimizer) <- function(optimizer, maxit)
  S7::set_props(optimizer, maxit = maxit)


#' Rebuild an Optimiser With the Trace Switched On
#' @param optimizer The \code{\link{optimizer}}.
#' @return An optimiser of the same class.
#' @keywords internal
with_trace <- S7::new_generic("with_trace", "optimizer",
                              function(optimizer) S7::S7_dispatch())

S7::method(with_trace, optimizer) <- function(optimizer)
  S7::set_props(optimizer, keep_trace = TRUE)


#' Rebuild an Optimiser With a Different Stopping Rule
#'
#' @description
#' Replaces the criterion, and for a wrapper replaces the one that will actually
#' be consulted.
#'
#' @details
#' The distinction matters. \code{\link{multistart}} carries a criterion only so
#' that printing it tells the truth; the rule that is evaluated belongs to the
#' optimiser inside. Setting the outer one and expecting a different run is the
#' sort of thing that makes a check pass while testing nothing.
#'
#' @param optimizer The \code{\link{optimizer}}.
#' @param criterion The new rule.
#'
#' @return An optimiser of the same class.
#'
#' @keywords internal
with_criterion <- S7::new_generic("with_criterion", "optimizer",
                                  function(optimizer, criterion)
                                    S7::S7_dispatch())

S7::method(with_criterion, optimizer) <- function(optimizer, criterion)
  S7::set_props(optimizer, criterion = criterion)
