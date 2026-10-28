#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
#' @include check_optimizer.R
NULL

#' @title S7 Class for Multi-Start
#'
#' @description
#' A `MultiStart` object wraps another optimizer and runs it from several
#' starting points. It is built by [multistart()]. Its seven shared
#' properties describe the wrapper: `criterion` and `max_eval` are copied
#' from the inner optimizer, so that the printed stopping rule and evaluation
#' budget are those in force; `verbose`, `refresh` and `keep_trace` are the
#' arguments of [multistart()]; and `maxit` counts starts.
#'
#' @details
#' Beyond the seven shared properties a `MultiStart` carries six of its own:
#' `optimizer` is the wrapped optimizer, `n`, `starts` and `spread` describe
#' where the runs begin, `ncores` how they are distributed over processes,
#' and `distinct_tol` how close two values must be to count as one optimum.
#'
#' Because the rule that is evaluated belongs to the inner optimizer,
#' [with_criterion()] has a method for this class that sets both; setting the
#' outer property without the inner one changes only the printing.
#'
#' @param optimizer The inner [optimizer()], run from each starting point.
#' @param n How many starts.
#' @param starts An optional matrix of starting points, one per row.
#' @param spread How widely the random starts are scattered.
#' @param ncores The number of processes over which the starts are spread.
#' @param distinct_tol Objective values closer than this count as one
#'   optimum.
#'
#' @return An S7 object of class `MultiStart` inheriting from [optimizer()],
#'   with the six properties above beside the seven shared ones.
#'
#' @seealso [multistart()] for the constructor, [chain()] for the other
#'   wrapper, [with_criterion()] for the reason why setting the outer rule is
#'   not enough.
#' @name MultiStart-class
#' @aliases MultiStart
#' @keywords internal
MultiStart <- S7::new_class("MultiStart", parent = optimizer,
  properties = list(
    optimizer    = S7::class_any,
    n            = S7::class_numeric,
    starts       = S7::class_any,
    spread       = S7::class_numeric,
    ncores       = S7::class_any,
    distinct_tol = S7::class_numeric
  ))


#' @title Run an Optimizer From Many Starting Points
#'
#' @description
#' Wraps any optimizer and runs it from several starting points, returning
#' the best result together with the number of distinct optimal values
#' reached. Two runs count as reaching the same optimum when their objective
#' values differ by at most `distinct_tol`. A count above one shows that the
#' surface has more than one local minimum; a count of one does not exclude
#' other minima with the same value or minima that all the starts missed.
#'
#' @param optimizer The [optimizer()] to run, any of them, including another
#'   `multistart()`. Anything else raises an error naming `bfgs()` as an
#'   example.
#' @param n How many starts, the caller's own `par` among them. Defaults to
#'   `10`. It becomes the `maxit` of the wrapper, so `n = 0` is rejected with
#'   an error that names `maxit`.
#' @param starts An optional matrix of starting points, one per row, used
#'   verbatim; `n` and `spread` are then ignored and `n` is taken from the
#'   number of rows. A matrix whose number of columns differs from the length
#'   of `par` is rejected when the run starts, with an error whose message
#'   gives the required number of columns.
#' @param spread How widely the random starts are scattered, in units of the
#'   unconstrained scale. Defaults to `1`.
#' @param ncores The number of processes over which to spread the starts. Defaults to
#'   `NULL`, meaning `min(n, max(1, parallel::detectCores() - 2))`. Pass `1`
#'   for a sequential run. For a deterministic inner optimizer the answer
#'   does not depend on the value; for a stochastic one, such as [sa()], it
#'   does.
#' @param distinct_tol Objective values differing by less than this count as
#'   the same optimum. Defaults to `1e-6`.
#' @param verbose Report each start as it finishes? Defaults to `FALSE`. Only
#'   a sequential run reports; a run over several processes prints
#'   nothing.
#' @param refresh Report every this many starts. Defaults to `1`.
#' @param keep_trace Keep the per-start summary? Defaults to `TRUE`, unlike
#'   every other optimizer, because here the trace is one row per start,
#'   not one per iteration.
#'
#' @details
#' # What the result reports
#'
#' The result has the `par`, `value`, `gradient`, `converged` and
#' `criterion_met` of the best run. Its `counts` are the sums over all
#' successful starts, its `iterations` is the number of starts, and its
#' `message` and `trace` are those of the wrapper. The summary of the starts
#' is in `trace`, one row for each start, with the columns `start`, `value`,
#' `converged` and `iterations`. The `message` counts the starts that
#' succeeded and those that converged, the distinct optimal values found and
#' how often the best was reached, as in
#'
#' ```
#' 12 starts, 12 succeeded, 12 converged, 2 distinct optima;
#' the best was found 6 times.
#' ```
#'
#' The count is computed by sorting the values reached and cutting wherever
#' consecutive ones differ by more than `distinct_tol`.
#'
#' # Starting points
#'
#' Unless `starts` is given, the first starting point is `par`, and the
#' remaining \eqn{n - 1} form a Latin hypercube: the range of each coordinate
#' is divided into equal strata and each stratum is used exactly once, which
#' spreads the starts more evenly than independent draws. They are generated
#' on the unconstrained scale and mapped back through the bounds, so
#' every start is admissible by construction, however tight the box. The
#' rows of a matrix given in `starts` are not checked against the bounds in
#' advance: a row outside the box fails in its own run and is recorded as a
#' failed start.
#'
#' # Parallel execution
#'
#' The starts are independent and are run over `ncores` processes. The
#' function creates the workers, loads the package, assigns the random
#' streams and shuts the workers down internally: on Unix-alikes the
#' workers are forks, on Windows they form a socket cluster. If the workers
#' cannot load \pkg{optimizers7}, the run gives a warning and proceeds
#' sequentially. If \pkg{optimizers7} was loaded from source with
#' \pkg{pkgload}, the starts are run sequentially on every platform, without
#' a warning, because a worker would load an installed copy of the package
#' instead of the copy in this session.
#'
#' Processes are used because the stopping rule is an R object evaluated at
#' every iteration, and R cannot be called from several threads.
#'
#' The starting points are drawn in the calling session before dispatch, and
#' each worker receives a random stream derived from the seed of the session,
#' so [set.seed()] reproduces the run. When the inner optimizer is
#' deterministic, such as [bfgs()] or [nelder_mead()], the result is the same
#' for every value of `ncores`. When the inner optimizer draws random
#' numbers, such as [sa()], the result depends on the assignment of the
#' streams to the workers, and therefore on `ncores`.
#'
#' # Failed starts
#'
#' A start where the objective is undefined is recorded as failed and the
#' others proceed; the `message` then reports how many succeeded and quotes
#' the first failure. An error is raised only when every start fails, and
#' it quotes the first error.
#'
#' @return An S7 object of class [MultiStart], inheriting from [optimizer()],
#'   to be handed to [minimize()].
#'
#' @examples
#' multistart(bfgs())
#'
#' # The surface has two minima, and one of them is lower. A single run from the
#' # origin finds the lower one; the count shows that there is another.
#' f <- function(p) (p[1]^2 - 1)^2 + p[2]^2 + 0.3 * p[1]
#' set.seed(1)
#' r <- minimize(multistart(bfgs(), n = 12), f, c(0, 0))
#' r@message
#' table(round(r@trace$value, 6))
#' r@par
#'
#' # With a deterministic inner optimizer the answer does not depend on the
#' # number of processes.
#' set.seed(11); one  <- minimize(multistart(bfgs(), n = 8, ncores = 1), f, c(0, 0))
#' set.seed(11); many <- minimize(multistart(bfgs(), n = 8, ncores = 2), f, c(0, 0))
#' identical(one@value, many@value)
#'
#' # A start where the objective is undefined is recorded, not fatal.
#' g <- function(p) if (p[1] > 0.5) stop("undefined here") else sum(p^2)
#' set.seed(3)
#' minimize(multistart(bfgs(), n = 8, ncores = 1), g, c(0, 0))@message
#'
#' # If every start fails, the run raises an error that quotes the first
#' # failure.
#' try(minimize(multistart(bfgs(), n = 4, ncores = 1),
#'              function(p) stop("never works"), c(0, 0)))
#'
#' # Starting points may be given outright, and then n is the row count.
#' S <- rbind(c(-2, 0), c(2, 0), c(0, 3))
#' minimize(multistart(bfgs(), starts = S, ncores = 1), f, c(0, 0))@trace
#'
#' @seealso [chain()] for running optimizers in sequence instead of in
#'   parallel, [sa()] for a global search that needs no restarts,
#'   [start_runif()] for another way of generating starting points.
#' @references
#' McKay, M. D., Beckman, R. J. and Conover, W. J. (1979). A
#' comparison of three methods for selecting values of input variables
#' in the analysis of output from a computer code.
#' *Technometrics* **21**, 239--245.
#'
#' @export
multistart <- function(optimizer, n = 10, starts = NULL, spread = 1,
                       ncores = NULL, distinct_tol = 1e-6,
                       verbose = FALSE, refresh = 1, keep_trace = TRUE) {
  if (!S7::S7_inherits(optimizer, optimizer_class())) {
    stop("'optimizer' must be an optimizer object, e.g. bfgs().", call. = FALSE)
  }
  if (!is.null(starts)) {
    if (!is.matrix(starts) || !nrow(starts)) {
      stop("'starts' must be a matrix with one starting point per row.",
           call. = FALSE)
    }
    n <- nrow(starts)
  }
  # Its budgets are the inner optimizer's, except maxit, which counts starts.
  check_optimizer_args(optimizer@criterion, n, optimizer@max_eval, verbose,
                       refresh, keep_trace)
  check_tol(spread, "spread")
  check_tol(distinct_tol, "distinct_tol")
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) || length(ncores) != 1L || is.na(ncores) ||
        ncores < 1 || !is.finite(ncores) || ncores != round(ncores)) {
      stop("'ncores' must be a single positive whole number, or NULL.",
           call. = FALSE)
    }
    ncores <- as.integer(ncores)
  }

  MultiStart(
    name = paste0("multistart (", optimizer@name, ")"),
    criterion = optimizer@criterion,
    maxit = n, max_eval = optimizer@max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    optimizer = optimizer, n = n, starts = starts, spread = spread,
    ncores = ncores, distinct_tol = distinct_tol
  )
}

#' @title What Multi-Start Can Offer a Stopping Rule
#' @name optimizer_provides.MultiStart
#' @description
#' Returns the answer of the inner optimizer, since the inner optimizer
#' evaluates the rule.
#' @param optimizer A `MultiStart` object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, MultiStart) <- function(optimizer)
  optimizer_provides(optimizer@optimizer)


#' Whether a Multistart Accepts Box Bounds
#'
#' @name optimizer_bounded.MultiStart
#' @description
#' Returns the answer of the inner optimizer, since every start is a run of
#' that optimizer with the same bounds. A multistart over [prox_grad()] therefore
#' does not accept box bounds.
#' @param optimizer A `MultiStart` object.
#' @return A single logical.
#' @keywords internal
S7::method(optimizer_bounded, MultiStart) <- function(optimizer)
  optimizer_bounded(optimizer@optimizer)


# Changing a setting has to reach the optimizer that will actually use it. The
# criterion is the one that matters: the outer copy exists so that printing a
# multistart tells the truth, but the rule that is evaluated belongs to the
# optimizer inside, so setting only the outer one would change nothing while
# looking as though it had.

#' @rdname with_criterion
#' @name with_criterion.MultiStart
#' @keywords internal
S7::method(with_criterion, MultiStart) <- function(optimizer, criterion) {
  S7::set_props(optimizer, criterion = criterion,
                optimizer = with_criterion(optimizer@optimizer, criterion))
}

#' @rdname with_maxit
#' @name with_maxit.MultiStart
#' @keywords internal
S7::method(with_maxit, MultiStart) <- function(optimizer, maxit) {
  # Its own maxit counts starts; the budget being varied is the inner one.
  S7::set_props(optimizer, optimizer = with_maxit(optimizer@optimizer, maxit))
}


#' Latin Hypercube Starting Points
#'
#' @description
#' Scatters `n - 1` starting points around `par`, generated on the
#' unconstrained scale so that any bounds are respected by construction.
#'
#' @details
#' The range of each coordinate is cut into `n - 1` strata, and each stratum
#' is used exactly once, so the starts cannot all fall in one corner, as independent
#' draws can. On the unconstrained scale the range is the starting value
#' \eqn{\eta_0} plus or minus
#' \eqn{3\, \texttt{spread}\, \max(1, \lvert\eta_0\rvert)}. The
#' unconstrained scale is the parameter itself for an unbounded coordinate
#' and its shifted log or scaled logit for a bounded one, so a start for a
#' variance is drawn on the log scale and comes back positive, without a
#' rejected draw.
#'
#' @param par The starting value of the caller.
#' @param n Total number of starts, including `par`.
#' @param spread Scale of the half-width of the sampling range in
#'   unconstrained units: the half-width is
#'   \eqn{3\, \texttt{spread}\, \max(1, \lvert\eta_0\rvert)}.
#' @param bounds Box constraints in the shape that [check_bounds()]
#'   returns, possibly empty.
#'
#' @return A matrix with `n` rows, the first of which is `par`.
#'
#' @keywords internal
make_starts <- function(par, n, spread, bounds) {
  p <- length(par)
  out <- matrix(NA_real_, n, p)
  out[1, ] <- par
  if (n == 1L) return(out)

  m <- n - 1L
  for (j in seq_len(p)) {
    b <- if (length(bounds)) bounds[[j]] else c(-Inf, Inf)
    eta0 <- bounded_forward(b, par[j])
    half <- 3 * spread * max(1, abs(eta0))
    # One draw per stratum, in random order.
    u <- (sample.int(m) - 1 + stats::runif(m)) / m
    eta <- eta0 - half + 2 * half * u
    out[-1L, j] <- bounded_transform(b, eta)$h
  }
  out
}


#' @title Minimize From Many Starting Points
#' @name minimize.MultiStart
#'
#' @description
#' Runs [multistart()] on the objective: it generates or takes the starting
#' points, runs the inner optimizer from each, and assembles the result from
#' the best run.
#'
#' @details
#' Before the first start runs, the objective, the bounds and `par` are
#' validated once, so a bad criterion or a `par` outside its box is rejected
#' before the \eqn{n} runs are launched; a row of a user-supplied `starts`
#' matrix is not checked in advance, and a row outside the box fails in its
#' own run. The state of the generator is recorded for the `seed` of the
#' result. The gradient check is switched off for the inner runs, since the
#' generic has already made it at the `par` of the caller; without that the
#' same warning would be printed once per start.
#'
#' @param optimizer A `MultiStart` object.
#' @param fn,par,gr,hess,lower,upper,... As in [minimize()]. Unless `starts`
#'   was given, `par` is the first starting point and is used as given, and
#'   the others are generated around it. When `starts` was given, `par` is
#'   used only to check the objective and the bounds and to fix the number of
#'   parameters.
#'
#' @return An [optimizer_result()] with the `par`, `value`, `gradient`,
#'   `converged` and `criterion_met` of the best run, the `counts` summed
#'   over the successful starts, `iterations` equal to the number of starts,
#'   the summary of the starts in `trace` and the counts in `message`. `seed`
#'   is the state from which the whole run began.
#'
#' @keywords internal
S7::method(minimize, MultiStart) <-
  function(optimizer, fn, par, gr = NULL, hess = NULL,
           lower = -Inf, upper = Inf, ...) {
    inner <- optimizer@optimizer
    # Validate once, here, so that a bad criterion or a start outside its bounds
    # is refused before n runs are launched rather than n times over.
    prepare_objective(inner, fn, par, gr, hess)
    bx <- check_bounds(lower, upper, par)

    seed <- capture_seed()
    S <- optimizer@starts
    if (is.null(S)) {
      S <- make_starts(as.numeric(par), as.integer(optimizer@n),
                       optimizer@spread, bx)
    } else if (ncol(S) != length(par)) {
      stop("'starts' must have ", length(par), " columns, one per parameter.",
           call. = FALSE)
    }

    n <- nrow(S)
    # The generic checked fn against gr once already, at the caller's par;
    # repeating it from every start would print the same warning n times.
    old_opt <- options(optimizers7.check_gradient = FALSE)
    on.exit(options(old_opt), add = TRUE)
    one <- function(i) {
      tryCatch(minimize(inner, fn, S[i, ], gr = gr, hess = hess,
                        lower = lower, upper = upper),
               error = function(e) conditionMessage(e))
    }

    t0 <- proc.time()[["elapsed"]]
    res <- run_starts(one, n, resolve_ncores(optimizer@ncores, n),
                      optimizer@verbose, optimizer@refresh)
    elapsed <- proc.time()[["elapsed"]] - t0

    build_multistart_result(res, S, optimizer, seed, elapsed)
  }


#' Record the Random Number Generator's State
#'
#' @description
#' Returns the value of `.Random.seed` before a run that draws random
#' numbers, so that the run can be repeated exactly.
#'
#' @details
#' Assigning the recorded value back into the global environment reproduces
#' the run. When `.Random.seed` does not exist yet, one uniform is drawn to
#' make \R create it, which costs a single discarded number on the first
#' stochastic call of a session.
#'
#' The stream is not created with `set.seed(NULL)`, which reseeds from
#' the clock and discards the seed that the caller had set.
#'
#' @return An integer vector, the value of `.Random.seed`.
#'
#' @keywords internal
capture_seed <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}


#' How Many Processes to Use
#'
#' @description
#' Turns `ncores = NULL` into a number: as many processes as there are
#' starts, but not more than the number of cores minus two.
#'
#' @details
#' The rule is `min(n, max(1, detectCores() - 2))`. Two cores are left free,
#' one for the calling session and one for the rest of the machine, which
#' otherwise stops responding. Starting more processes than there are starts
#' costs the start-up time of the extra processes, which on Windows is
#' seconds, whereas a fork takes microseconds.
#'
#' It is capped at two under `R CMD check`, which sets
#' `_R_CHECK_LIMIT_CORES_` and fails a package that ignores it, and it
#' falls back to one process where `detectCores()` returns `NA`.
#'
#' @param ncores The value given by the caller, possibly `NULL`.
#' @param n The number of starts.
#'
#' @return A single integer, at least one.
#'
#' @keywords internal
resolve_ncores <- function(ncores, n) {
  if (is.null(ncores)) {
    avail <- suppressWarnings(parallel::detectCores())
    if (!is.finite(avail)) avail <- 1L
    ncores <- min(n, max(1L, avail - 2L))
  }
  ncores <- min(as.integer(ncores), n)
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) ncores <- min(ncores, 2L)
  max(1L, ncores)
}


#' Run the Starts, in Parallel or Not
#'
#' @description
#' Evaluates `one(i)` for `i` in `1:n`, over `ncores`
#' processes, and cleans up after itself.
#'
#' @details
#' The function selects one of three routes, so that the caller does not
#' have to.
#'
#' One process is the sequential loop, and it is the only route that reports
#' progress during the run, since a worker has no output that the caller
#' would see. When \pkg{optimizers7} was loaded with \pkg{pkgload}, this
#' route is used whatever `ncores` is.
#'
#' On a Unix-alike the workers are forks, through `parallel::mclapply()`.
#' A fork starts in microseconds and inherits the whole session, so the
#' function does not have to load the package or export objects.
#'
#' On Windows there is no `fork`, so a socket cluster is started here and
#' stopped on exit. Its workers are new sessions, so \pkg{optimizers7} has to
#' be loaded on them; when that is not possible, because the package is not
#' installed in a library that they can see, the function gives a warning
#' and runs the starts sequentially instead of failing.
#'
#' Both parallel routes seed their workers from the stream of this session,
#' so `set.seed()` reproduces the run. For a deterministic inner optimizer
#' the result is the same for every `ncores`; for one that draws random
#' numbers, such as [sa()], it depends on the assignment of the streams to
#' the workers. The generator is set to `"L'Ecuyer-CMRG"` for the duration
#' and restored afterwards, since it is the only kind that \R can split into
#' independent streams.
#'
#' @param one A function of the index of the start, returning a result or a
#'   message.
#' @param n How many starts.
#' @param ncores How many processes, already resolved.
#' @param verbose Report each start as it finishes?
#' @param refresh Report every this many.
#'
#' @return A list of `n` results.
#'
#' @keywords internal
run_starts <- function(one, n, ncores, verbose, refresh) {
  sequential <- function() {
    res <- vector("list", n)
    for (i in seq_len(n)) {
      res[[i]] <- one(i)
      if (verbose && refresh > 0 && i %% refresh == 0) {
        v <- if (is.character(res[[i]])) "failed"
             else format(res[[i]]@value, digits = 8)
        cat(sprintf("  start %4d of %d : %s\n", i, n, v))
      }
    }
    res
  }
  if (ncores <= 1L) return(sequential())

  # Independent streams for the workers, derived from the caller's seed so that
  # set.seed() still governs everything.
  old_kind <- RNGkind()[1]
  on.exit(RNGkind(old_kind), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  iseed <- sample.int(.Machine$integer.max, 1L)

  # A package loaded from source with pkgload exists only in this session: a
  # worker process would load whatever copy is installed, which may be another
  # version entirely, and S7 objects built here do not dispatch correctly
  # against methods registered there. Sequential is the only correct route.
  if (isNamespaceLoaded("optimizers7") &&
      exists(".__DEVTOOLS__", asNamespace("optimizers7"),
             inherits = FALSE)) {
    return(sequential())
  }

  if (.Platform$OS.type != "windows") {
    set.seed(iseed)
    return(parallel::mclapply(seq_len(n), one, mc.cores = ncores,
                              mc.set.seed = TRUE))
  }

  cl <- try(parallel::makePSOCKcluster(ncores), silent = TRUE)
  if (inherits(cl, "try-error")) {
    warning("Could not start ", ncores, " worker processes; running the ",
            "starts sequentially.", call. = FALSE)
    return(sequential())
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)

  have <- try(parallel::clusterEvalQ(
    cl, requireNamespace("optimizers7", quietly = TRUE)), silent = TRUE)
  if (inherits(have, "try-error") || !all(unlist(have))) {
    warning("The worker processes could not load optimizers7, so the starts ",
            "were run\n  sequentially. They are separate R sessions and do ",
            "not inherit this one's\n  loaded packages, so a package loaded ",
            "with pkgload, or installed in a\n  library not on their ",
            ".libPaths(), is invisible to them.", call. = FALSE)
    return(sequential())
  }
  parallel::clusterSetRNGStream(cl, iseed)
  parallel::parLapply(cl, seq_len(n), one)
}


#' Assemble the Result of a Multi-Start Run
#'
#' @description
#' Picks the run with the lowest value and takes its `par`, `value`,
#' `gradient`, `converged` and `criterion_met`. The evaluation counts are
#' summed over the successful starts, `iterations` is the number of starts,
#' and the table of the starts, the counts of starts and the number of
#' distinct optimal values are added in `trace` and `message`.
#'
#' @details
#' The number of distinct optimal values is computed by sorting the values
#' reached and cutting wherever consecutive values differ by more than
#' `distinct_tol`. It describes the objective and not the optimizer, and a
#' single run cannot provide it.
#'
#' @param res The list of results, in which the entry of a failed start is a
#'   character message.
#' @param S The matrix of starting points.
#' @param optimizer The `MultiStart` object.
#' @param seed The generator state with which the run began.
#' @param elapsed Seconds.
#'
#' @return An [optimizer_result()].
#'
#' @keywords internal
build_multistart_result <- function(res, S, optimizer, seed, elapsed) {
  ok <- !vapply(res, is.character, logical(1))
  if (!any(ok)) {
    stop("Every start failed. The first error was: ", res[[1]], call. = FALSE)
  }

  vals <- rep(NA_real_, length(res))
  conv <- rep(NA, length(res))
  iters <- rep(NA_integer_, length(res))
  vals[ok] <- vapply(res[ok], function(r) r@value, numeric(1))
  conv[ok] <- vapply(res[ok], function(r) r@converged, logical(1))
  iters[ok] <- vapply(res[ok], function(r) as.integer(r@iterations), integer(1))

  best <- which.min(replace(vals, is.na(vals), Inf))
  out <- res[[best]]

  # Distinct optima: sort the values reached and cut wherever consecutive ones
  # separate by more than the tolerance.
  v <- sort(vals[ok])
  n_distinct <- if (length(v) <= 1L) length(v)
                else 1L + sum(diff(v) > optimizer@distinct_tol)
  n_best <- sum(abs(vals[ok] - vals[best]) <= optimizer@distinct_tol)

  msg <- sprintf(
    "%d starts, %d succeeded, %d converged, %d distinct optim%s; the best was found %d time%s",
    length(res), sum(ok), sum(conv, na.rm = TRUE), n_distinct,
    if (n_distinct == 1L) "um" else "a", n_best, if (n_best == 1L) "" else "s")
  if (any(!ok)) {
    msg <- paste0(msg, ". First failure: ", sub("[.]$", "", res[!ok][[1]]))
  }
  if (nzchar(out@message)) msg <- paste0(msg, ". Best run: ", out@message)

  trace <- NULL
  if (optimizer@keep_trace) {
    trace <- data.frame(start = seq_along(res), value = vals,
                        converged = conv, iterations = iters)
  }

  counts <- c(f = 0, g = 0, h = 0)
  for (r in res[ok]) counts <- counts + r@counts

  optimizer_result(
    par = out@par, value = out@value, gradient = out@gradient,
    counts = counts,
    # One iteration of this optimizer is one start; the inner iteration counts
    # are in the trace, where they can be read per start rather than summed into
    # a number that means nothing.
    iterations = length(res),
    converged = out@converged,
    criterion_met = out@criterion_met,
    message = msg, trace = trace, optimizer = optimizer, elapsed = elapsed,
    seed = seed
  )
}
