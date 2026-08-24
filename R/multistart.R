#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
#' @include check_optimizer.R
NULL

#' @title S7 Class for Multi-Start
#'
#' @description
#' An optimizer wrapping another and running it from several starting points.
#' Built by [multistart()]. Its seven shared properties describe the wrapper
#' and not the inner run: `criterion`, `max_eval` and the rest are copied
#' from the optimizer inside so that printing tells the truth, while `maxit`
#' counts **starts**.
#'
#' @details
#' Beyond the seven shared properties a `MultiStart` carries six of its own:
#' `optimizer` is the one being wrapped, `n`, `starts` and `spread` say where
#' the runs begin, `ncores` how they are spread over processes, and
#' `distinct_tol` how close two answers must be to count as one.
#'
#' Because the rule that is evaluated belongs to the inner optimizer,
#' [with_criterion()] has a method for this class that sets both; setting the
#' outer property alone changes the printing and nothing else.
#'
#' @param optimizer The inner [optimizer()], run from each starting point.
#' @param n How many starts.
#' @param starts An optional matrix of starting points, one per row.
#' @param spread How widely the random starts are scattered.
#' @param ncores How many processes the starts are spread over.
#' @param distinct_tol Objective values closer than this count as one
#'   optimum.
#'
#' @return An S7 object of class `MultiStart` inheriting from [optimizer()],
#'   with the six properties above beside the seven shared ones.
#'
#' @seealso [multistart()] for the constructor, [chain()] for the other
#'   wrapper, [with_criterion()] for why setting the outer rule is not
#'   enough.
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
#' the best result together with the number of **distinct optima** found. That
#' count is evidence about the objective that no single run can supply: it
#' says whether the surface has one minimum or several.
#'
#' @param optimizer The [optimizer()] to run, any of them, including another
#'   `multistart()`. Anything else raises an error naming `bfgs()` as an
#'   example.
#' @param n How many starts, the caller's own `par` among them. Defaults to
#'   `10`. It becomes the wrapper's `maxit`, so `n = 0` is refused with a
#'   message naming `maxit`.
#' @param starts An optional matrix of starting points, one per row, used
#'   verbatim; `n` and `spread` are then ignored and `n` is taken from the
#'   number of rows. A matrix whose column count does not match the parameter
#'   vector is refused when the run starts, naming both.
#' @param spread How widely the random starts are scattered, in units of the
#'   unconstrained scale. Defaults to `1`.
#' @param ncores How many processes to spread the starts over. Defaults to
#'   `NULL`, meaning `min(n, max(1, parallel::detectCores() - 2))`. Pass `1`
#'   for a sequential run; the answer does not depend on the value.
#' @param distinct_tol Objective values differing by less than this count as
#'   the same optimum. Defaults to `1e-6`.
#' @param verbose Report each start as it finishes? Defaults to `FALSE`.
#' @param refresh Report every this many starts. Defaults to `1`.
#' @param keep_trace Keep the per-start summary? Defaults to `TRUE`, unlike
#'   every other optimizer, because here the trace is one row per **start**,
#'   not one per iteration.
#'
#' @details
#' # What the result reports
#'
#' The result is the best run, with everything that run carried. The
#' per-start summary is in `trace`, one row each with the columns `start`,
#' `value`, `converged` and `iterations`. The `message` counts the starts
#' that succeeded, the ones that converged, the distinct optima found and how
#' often the best was reached, as in
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
#' The first starting point is always `par`. The remaining \eqn{n - 1} form a
#' Latin hypercube: each coordinate's range is divided into equal strata and
#' each stratum is used exactly once, which spreads the starts more evenly
#' than independent draws would. They are generated on the **unconstrained**
#' scale and mapped back through the bounds, so every start is admissible by
#' construction however tight the box.
#'
#' # Parallel execution
#'
#' The starts are independent and are run over `ncores` processes. Worker
#' creation, package loading, random-stream assignment and shutdown are
#' handled internally: on Unix-alikes the workers are forks, on Windows a
#' socket cluster. If the workers cannot load the package the run warns and
#' proceeds sequentially. That is what happens under `pkgload`, the workers
#' being separate sessions that do not inherit this one's loaded packages.
#'
#' Processes are used because the stopping rule is an R object consulted at
#' every iteration, and R cannot be called from several threads.
#'
#' The starting points are drawn in the calling session before dispatch, and
#' each worker receives a random stream derived from the session's seed, so
#' [set.seed()] reproduces the run identically at any `ncores`: measured on
#' the example below, `ncores` of 1, 2 and 4 give `-0.3054284837` and the
#' same distinct count.
#'
#' # Failed starts
#'
#' A start where the objective is undefined is recorded as failed and the
#' rest proceed; the `message` then reports how many succeeded and quotes the
#' first failure. An error is raised only when **every** start fails, and it
#' quotes the first message.
#'
#' @return An S7 object of class [MultiStart], inheriting from [optimizer()],
#'   to be handed to [minimize()].
#'
#' @examples
#' multistart(bfgs())
#'
#' # A surface with two minima, one of them better. A single run from the
#' # origin happens to find the good one; the count is what says there is
#' # another.
#' f <- function(p) (p[1]^2 - 1)^2 + p[2]^2 + 0.3 * p[1]
#' set.seed(1)
#' r <- minimize(multistart(bfgs(), n = 12), f, c(0, 0))
#' r@message
#' table(round(r@trace$value, 6))
#' r@par
#'
#' # The answer does not depend on how many processes ran it.
#' set.seed(11); one  <- minimize(multistart(bfgs(), n = 8, ncores = 1), f, c(0, 0))
#' set.seed(11); many <- minimize(multistart(bfgs(), n = 8, ncores = 2), f, c(0, 0))
#' identical(one@value, many@value)
#'
#' # A start where the objective is undefined is recorded, not fatal.
#' g <- function(p) if (p[1] > 0.5) stop("undefined here") else sum(p^2)
#' set.seed(3)
#' minimize(multistart(bfgs(), n = 8, ncores = 1), g, c(0, 0))@message
#'
#' # Every start failing is fatal, and the message says why.
#' try(minimize(multistart(bfgs(), n = 4, ncores = 1),
#'              function(p) stop("never works"), c(0, 0)))
#'
#' # Starting points may be given outright, and then n is the row count.
#' S <- rbind(c(-2, 0), c(2, 0), c(0, 3))
#' minimize(multistart(bfgs(), starts = S, ncores = 1), f, c(0, 0))@trace
#'
#' @seealso [chain()] for running optimizers in sequence rather than in
#'   parallel, [sa()] for a global search that needs no restarts,
#'   [start_runif()] for the other way of generating starting points.
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
  check_tol(spread)
  check_tol(distinct_tol)
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) || length(ncores) != 1L || is.na(ncores) ||
        ncores < 1 || ncores != round(ncores)) {
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
#' Whatever the inner optimizer offers, since it is the inner optimizer that
#' evaluates the rule.
#' @param optimizer A `MultiStart` object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, MultiStart) <- function(optimizer)
  optimizer_provides(optimizer@optimizer)


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
#' Each coordinate's range is cut into `n - 1` strata, each used exactly
#' once, so the starts cannot all fall in one corner the way independent draws
#' can. The range is `par` plus or minus `3 * spread` on the
#' unconstrained scale, which for an unbounded parameter is the parameter itself
#' and for a bounded one is its log or logit. A start for a variance is
#' therefore drawn as a log and comes back positive with no rejected draw.
#'
#' @param par The user's starting value.
#' @param n Total number of starts, including `par`.
#' @param spread Half-width of the sampling range, in unconstrained units.
#' @param bounds Box constraints in the shape [check_bounds()]
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
#' Runs [multistart()] on the objective: generate or take the starting
#' points, run the inner optimizer from each, and assemble the best.
#'
#' @details
#' Three things happen before the first start runs. The objective and the
#' bounds are validated **once**, so a bad criterion or a start outside its
#' box is refused before \eqn{n} runs are launched rather than \eqn{n} times
#' over. The generator state is captured for the result's `seed`. And the
#' gradient-consistency check is switched off for the duration, the generic
#' having already made it at the caller's `par`; without that the same
#' warning would print once per start.
#'
#' @param optimizer A `MultiStart` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `par` is the first
#'   starting point and is used as given; the rest are generated around it.
#'
#' @return An [optimizer_result()]: the best run, carrying its own `par`,
#'   `value`, `gradient` and `converged`, with the per-start summary in
#'   `trace` and the counts in `message`. `seed` is the state the whole run
#'   began from.
#'
#' @keywords internal
S7::method(minimize, MultiStart) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    inner <- optimizer@optimizer
    # Validate once, here, so that a bad criterion or a start outside its bounds
    # is refused before n runs are launched rather than n times over.
    prepare_objective(inner, fn, par, gr, he)
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
      tryCatch(minimize(inner, fn, S[i, ], gr = gr, he = he,
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
#' The state `.Random.seed` held before a run that draws random numbers, so
#' that the run can be repeated exactly.
#'
#' @details
#' Assigning the recorded value back into the global environment reproduces the
#' run. When no stream exists yet one uniform is drawn to force \R to create
#' one, which costs a single discarded number on the first stochastic call of a
#' session.
#'
#' It is emphatically **not** done with `set.seed(NULL)`, which
#' reseeds from the clock and throws away whatever the caller had set. That is a
#' defect this toolkit has already met once: it made a check in \pkg{distributions7}
#' silently random, and the resulting test failed about one CI run in several
#' hundred, on whichever platform happened to draw it.
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
#' starts, but never more than the machine can spare.
#'
#' @details
#' The rule is `min(n, max(1, detectCores() - 2))`. Two are held back
#' and not one, because the session doing the asking is itself one of them,
#' and a machine with nothing left over is a machine that stops responding.
#' Asking for more processes than there are starts wastes the cost of starting
#' them, which on Windows is seconds rather than microseconds.
#'
#' It is capped at two under `R CMD check`, which sets
#' `_R_CHECK_LIMIT_CORES_` and fails a package that ignores it, and it
#' falls back to one process where `detectCores()` cannot tell.
#'
#' @param ncores What the caller asked for, possibly `NULL`.
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
#' Three routes, chosen here so that the caller need not.
#'
#' One process is the sequential loop, and it is the only route that can report
#' progress as it goes, a worker having nowhere to print to that the caller
#' would see.
#'
#' On a Unix-alike the workers are **forks**, through
#' `parallel::mclapply()`. A fork starts in microseconds and inherits this
#' session entire, so there is nothing to load and nothing to export. That
#' includes a package loaded with \pkg{pkgload}, so this route works during
#' development where a socket cluster does not.
#'
#' On Windows there is no `fork`, so a **socket cluster** is started
#' here and stopped on exit. Its workers are fresh sessions that know nothing, so
#' \pkg{optimizers7} has to be loaded on them; when it cannot be, because it is
#' not installed anywhere they can see, this warns and runs sequentially rather
#' than failing. A slower answer beats an error about `checkForRemoteErrors`
#' for someone who only asked for several starting points.
#'
#' Both parallel routes seed their workers from this session's stream, so
#' `set.seed()` reproduces the run and reproduces it identically whatever
#' `ncores` was. The generator is set to `"L'Ecuyer-CMRG"` for the
#' duration and put back afterwards, that being the only kind \R can split into
#' independent streams.
#'
#' @param one A function of the start's index, returning a result or a message.
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
#' Picks the run with the lowest value, keeps everything that run carried,
#' and adds what the others found: the per-start table, the counts, and the
#' number of distinct optima.
#'
#' @details
#' The count of distinct optima is the reason to run this at all, so it is
#' computed here: the values reached are sorted and
#' cut wherever consecutive ones differ by more than `distinct_tol`. It is
#' a statement about the objective, not about the optimizer, and it is the one
#' piece of evidence a single run can never supply.
#'
#' @param res The list of results, entries that failed being character messages.
#' @param S The matrix of starting points.
#' @param optimizer The `MultiStart` object.
#' @param seed The generator state the run began with.
#' @param elapsed Seconds.
#'
#' @return An [optimizer_result()].
#'
#' @keywords internal
build_multistart_result <- function(res, S, optimizer, seed, elapsed) {
  ok <- !vapply(res, is.character, logical(1))
  if (!any(ok)) {
    stop("Every start failed. The first said: ", res[[1]], call. = FALSE)
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
    msg <- paste0(msg, ". First failure: ", res[!ok][[1]])
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
