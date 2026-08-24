#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
#' @include check_optimizer.R
NULL

#' @title S7 Class for a Sequence of Optimizers
#'
#' @description
#' An optimizer holding a list of optimizers to be run in order, each from
#' the point the previous one reached. Built by [chain()]. Because it
#' inherits from [optimizer()] it is accepted anywhere a single method is,
#' including inside [multistart()].
#'
#' @details
#' Beyond the seven properties every optimizer has, a `Chain` carries
#' `stages`, the list of optimizers in running order. The seven it inherits
#' describe the **last** stage: its `criterion`, `maxit` and `max_eval` are
#' copied from there, because the last stage is the one whose rule ends the
#' run and whose result is reported. `refresh` is fixed at 1, the chain's own
#' progress being one line per stage.
#'
#' @param stages A list of [optimizer()] objects, in the order they run.
#'
#' @return An S7 object of class `Chain` inheriting from [optimizer()], with
#'   the property `stages` beside the seven shared ones.
#'
#' @seealso [chain()] for the constructor, [minimize.Chain()] for the run.
#' @name Chain-class
#' @aliases Chain
#' @keywords internal
Chain <- S7::new_class("Chain", parent = optimizer,
  properties = list(stages = S7::class_list))


#' @title Run Optimizers One After Another
#'
#' @description
#' Builds one optimizer out of several, run in order, each starting from the
#' point the previous one reached. The composition a global search needs:
#' `chain(sa(), lbfgs())` explores and then descends from wherever the
#' exploration left off, and neither method has to know about the other. The
#' result is itself an [optimizer()], so it goes anywhere a single method
#' goes, `multistart(chain(sa(), lbfgs()))` included.
#'
#' @details
#' # Each stage keeps its own settings
#'
#' A stage carries its own stopping rule and its own budgets, which is the
#' reason a chain is built out of optimizers instead of taking arguments: a
#' coarse rule and a small budget for the exploration, a tight rule for the
#' descent, written as `chain(sa(maxit = 200), newton(criterion =
#' crit_grad(1e-12)))`.
#'
#' # What the result reports
#'
#' The point, the value, the gradient and `iterations` are the last stage's,
#' that being where the run ended. So is `converged`: a chain has converged
#' when the method that finished it says so. An earlier stage exhausting its
#' own budget is the ordinary way a global search ends and leaves the flag
#' alone.
#'
#' The evaluation counts are **summed over the stages**, so a chain can be
#' compared against a single optimizer on total work. `message` is the
#' stages' messages joined, each prefixed with its stage number.
#'
#' A trace is assembled when both the chain and the stage were built with
#' `keep_trace = TRUE`; the chain's flag decides whether one is kept at all
#' and each stage's decides whether that stage contributes rows. It carries a
#' `stage` column. Stages need not report the same columns, a derivative-free
#' method having a `stationarity` where a descent has a `gnorm`; when they
#' disagree only the last stage's trace is kept.
#'
#' A stage that raises propagates: a method that cannot run on the objective
#' is a fact about the objective. A stage that runs without converging simply
#' passes its point on, and the first stage of a chain usually does exactly
#' that. Every stage's rule is checked by [check_criterion()] before any of
#' them runs, so a chain whose third stage cannot evaluate its criterion
#' fails before spending the first two.
#'
#' `chain(x)` with a single stage is that stage's run reported through the
#' chain, identical to the stage's own down to the evaluation counts.
#'
#' @param ... Optimizers, in the order they should run. At least one is
#'   required; one is accepted. Anything that is not an [optimizer()] raises
#'   an error naming `sa()` and `bfgs()` as examples.
#' @param verbose `TRUE` or `FALSE`, default `FALSE`: report which stage is
#'   starting, one line each. Independent of the stages' own `verbose`.
#' @param keep_trace `TRUE` or `FALSE`, default `FALSE`: assemble a trace
#'   across the stages. The stages must also have been built with
#'   `keep_trace = TRUE` to contribute anything.
#'
#' @return An S7 object of class [Chain], inheriting from [optimizer()]. Its
#'   `name` is the stages' names joined by `then`, and its `criterion`,
#'   `maxit` and `max_eval` are copied from the last stage.
#'
#' @examples
#' chain(sa(maxit = 10), bfgs())
#'
#' # Rastrigin has 121 local minima, so which one a descent method finds is
#' # decided by where it began. The search picks the basin; the descent
#' # finishes the job inside it.
#' rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
#' rg <- function(p) 2 * p + 20 * pi * sin(2 * pi * p)
#'
#' set.seed(3)
#' minimize(bfgs(), rastrigin, c(3.5, -2.5), gr = rg)@value   # 12.93
#' set.seed(3)
#' minimize(sa(maxit = 20), rastrigin, c(3.5, -2.5))@value    # 1.48
#' set.seed(3)
#' minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))@value
#'
#' # The counts are the whole chain's work; the iteration count is the last
#' # stage's.
#' set.seed(3)
#' res <- minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))
#' unlist(res@counts)
#' res@iterations
#'
#' # One stage is not a special case: it reports exactly what the stage does.
#' identical(minimize(chain(bfgs()), rastrigin, c(3.5, -2.5), gr = rg)@par,
#'           minimize(bfgs(), rastrigin, c(3.5, -2.5), gr = rg)@par)
#'
#' # A stage that cannot evaluate its own rule is caught before anything runs.
#' try(minimize(chain(bfgs(), nelder_mead(criterion = crit_grad())),
#'              rastrigin, c(1, 1), gr = rg))
#'
#' @seealso [multistart()] for the other wrapper, [sa()] for the search a
#'   chain usually opens with, [check_criterion()] for the rule check.
#' @export
chain <- function(..., verbose = FALSE, keep_trace = FALSE) {
  stages <- list(...)
  if (!length(stages)) {
    stop("'chain' needs at least one optimizer.", call. = FALSE)
  }
  ok <- vapply(stages, function(o) S7::S7_inherits(o, optimizer_class()),
               logical(1))
  if (!all(ok)) {
    stop("every stage must be an optimizer object, e.g. sa() or bfgs().",
         call. = FALSE)
  }
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(keep_trace) != 1L || !is.logical(keep_trace) ||
      is.na(keep_trace)) {
    stop("'keep_trace' must be TRUE or FALSE.", call. = FALSE)
  }
  last <- stages[[length(stages)]]
  Chain(
    name = paste(vapply(stages, function(o) o@name, character(1)),
                 collapse = " then "),
    # The rule and the budgets shown are the last stage's, because the result
    # is reported from there: check_criterion() then asks the question that
    # decides whether the reported run could evaluate its own rule.
    criterion = last@criterion, maxit = last@maxit, max_eval = last@max_eval,
    verbose = verbose, refresh = 1, keep_trace = keep_trace,
    stages = stages
  )
}


#' @title What a Chain Can Offer a Stopping Rule
#' @name optimizer_provides.Chain
#'
#' @description
#' Reports whatever the last stage offers. The chain's stopping rule is the
#' last stage's, so it is that stage which has to be able to evaluate it, and
#' a chain opening with a derivative-free search can still carry a gradient
#' rule.
#'
#' @param optimizer A `Chain` object.
#'
#' @return A character vector of `state` component names, the last stage's.
#'
#' @examples
#' # A search then a descent: the descent's gradient is what the rule reads.
#' optimizer_provides(chain(sa(), bfgs()))
#'
#' # The other order offers only the simplex-free measure.
#' optimizer_provides(chain(bfgs(), nelder_mead()))
#'
#' @keywords internal
S7::method(optimizer_provides, Chain) <- function(optimizer)
  optimizer_provides(optimizer@stages[[length(optimizer@stages)]])


#' @title Whether a Chain Takes Box Bounds
#' @name optimizer_bounded.Chain
#'
#' @description
#' `TRUE` when every stage takes bounds, `FALSE` otherwise. The bounds are
#' passed to all of them, so a single stage that would ignore them leaves the
#' chain unable to promise the box.
#'
#' @details
#' [prox_grad()] is the one shipped method that answers `FALSE`: it takes its
#' constraint inside the proximal operator, where the constraint composes
#' with the penalty already there. A chain containing it is therefore
#' unbounded whatever else is in it.
#'
#' @param optimizer A `Chain` object.
#'
#' @return A single logical.
#'
#' @examples
#' optimizer_bounded(chain(bfgs(), newton()))
#'
#' pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
#' optimizer_bounded(pg)
#' optimizer_bounded(chain(bfgs(), pg))
#'
#' @keywords internal
S7::method(optimizer_bounded, Chain) <- function(optimizer)
  all(vapply(optimizer@stages, optimizer_bounded, logical(1)))


#' @rdname with_criterion
#' @name with_criterion.Chain
#' @keywords internal
S7::method(with_criterion, Chain) <- function(optimizer, criterion) {
  # The rule that decides the run belongs to the last stage; setting only the
  # outer copy would change the printing and nothing else.
  st <- optimizer@stages
  st[[length(st)]] <- with_criterion(st[[length(st)]], criterion)
  S7::set_props(optimizer, criterion = criterion, stages = st)
}

#' @rdname with_maxit
#' @name with_maxit.Chain
#' @keywords internal
S7::method(with_maxit, Chain) <- function(optimizer, maxit) {
  st <- optimizer@stages
  st[[length(st)]] <- with_maxit(st[[length(st)]], maxit)
  S7::set_props(optimizer, maxit = maxit, stages = st)
}


#' @title Minimize by a Sequence of Optimizers
#' @name minimize.Chain
#'
#' @description
#' Runs the stages in order, handing each the point the previous one reached,
#' and assembles one result from them. The objective, its derivatives and the
#' bounds are passed to every stage unchanged; only the starting point moves.
#'
#' @details
#' Two things happen before the first stage runs. Every stage's stopping rule
#' is put to [check_criterion()], so a chain whose last stage cannot evaluate
#' its own rule fails without spending the earlier ones. And the
#' gradient-consistency check is switched off for the duration, the generic
#' having already made it once at the caller's `par`; without that, the same
#' warning would print once per stage.
#'
#' @param optimizer A `Chain` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. Passed to every
#'   stage as given, with `par` replaced by the previous stage's answer.
#'
#' @return An [optimizer_result()] carrying the last stage's point, value,
#'   gradient, iteration count and verdict, the summed evaluation counts, the
#'   stacked trace and the elapsed time of the whole chain. `seed` is the
#'   **first** stage's, that being the state a repeat of the chain has to
#'   start from.
#'
#' @keywords internal
S7::method(minimize, Chain) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    stages <- optimizer@stages
    # Every stage's rule is checked before any of them runs, so a chain whose
    # third stage cannot evaluate its criterion fails before spending the first
    # two rather than after.
    for (st in stages) check_criterion(st)

    # The generic checked fn against gr once already, at the caller's par;
    # letting each stage repeat it would print the same warning once per stage.
    old_opt <- options(optimizers7.check_gradient = FALSE)
    on.exit(options(old_opt), add = TRUE)

    x <- as.numeric(par)
    t0 <- proc.time()[["elapsed"]]
    results <- vector("list", length(stages))
    for (k in seq_along(stages)) {
      if (optimizer@verbose) {
        cat(sprintf("  stage %d/%d: %s\n", k, length(stages),
                    stages[[k]]@name))
      }
      res <- minimize(stages[[k]], fn, x, gr = gr, he = he,
                      lower = lower, upper = upper, ...)
      results[[k]] <- res
      x <- res@par
    }
    elapsed <- proc.time()[["elapsed"]] - t0

    build_chain_result(results, optimizer, elapsed)
  }


#' Assemble the Result of a Chain
#'
#' @description
#' Combines the per-stage results into one: the last stage's point, value,
#' gradient, iteration count and verdict, the evaluation counts added over
#' every stage, the traces stacked with a `stage` column, and the messages
#' joined with their stage numbers.
#'
#' @details
#' The counts are added with `Reduce("+")` rather than summed coordinate by
#' coordinate, so that the storage type survives: a chain of one stage must
#' report that stage's counts identically, and a `sum()` over a `vapply()`
#' would hand back doubles where the stage had integers.
#'
#' Traces are stacked only when every stage reports the same columns. A
#' derivative-free stage has a `stationarity` column where a descent has
#' `gnorm`, and when they disagree the last stage's trace is kept alone.
#'
#' @param results The per-stage [optimizer_result()] objects, in order.
#' @param optimizer The `Chain` that produced them.
#' @param elapsed Total seconds for the whole chain.
#'
#' @return An [optimizer_result()].
#'
#' @keywords internal
build_chain_result <- function(results, optimizer, elapsed) {
  last <- results[[length(results)]]
  # counts is a named c(f =, g =, h =): the work of the whole chain is what a
  # caller comparing it against a single optimizer wants to see. Added rather
  # than summed coordinate by coordinate so that the STORAGE TYPE survives --
  # a chain of one stage must report that stage's counts identically, and
  # sum() over a vapply would hand back a double where the stage had integers.
  counts <- Reduce(`+`, lapply(results, function(r) r@counts))

  tr <- NULL
  if (optimizer@keep_trace) {
    parts <- list()
    for (k in seq_along(results)) {
      t_k <- results[[k]]@trace
      if (is.data.frame(t_k) && nrow(t_k)) {
        t_k$stage <- k
        parts[[length(parts) + 1L]] <- t_k
      }
    }
    # the stages need not report the same columns: a derivative-free one has a
    # stationarity where a descent has a gradient norm, so they are stacked
    # only when they agree and otherwise the last stage's is kept
    if (length(parts)) {
      same <- length(unique(lapply(parts, names))) == 1L
      tr <- if (same) do.call(rbind, parts) else parts[[length(parts)]]
    }
  }

  msg <- vapply(seq_along(results), function(k) {
    m <- results[[k]]@message
    if (nzchar(m)) sprintf("stage %d: %s", k, m) else ""
  }, character(1))
  msg <- paste(msg[nzchar(msg)], collapse = "; ")

  optimizer_result(
    par = last@par, value = last@value, gradient = last@gradient,
    counts = counts, iterations = last@iterations,
    converged = last@converged, criterion_met = last@criterion_met,
    message = msg, trace = tr, optimizer = optimizer, elapsed = elapsed,
    seed = results[[1L]]@seed
  )
}
