#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
#' @include check_optimizer.R
NULL

#' @title S7 Class for Multi-Start
#' @description The class \code{\link{multistart}} instantiates.
#' @param optimizer The inner optimiser, run from each starting point.
#' @param n How many starts.
#' @param starts An optional matrix of starting points.
#' @param spread How widely the random starts are scattered.
#' @param cluster An optional \pkg{parallel} cluster.
#' @param distinct_tol Objective values closer than this count as one optimum.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{multistart}}
#' @keywords internal
MultiStart <- S7::new_class("MultiStart", parent = optimizer,
  properties = list(
    optimizer    = S7::class_any,
    n            = S7::class_numeric,
    starts       = S7::class_any,
    spread       = S7::class_numeric,
    cluster      = S7::class_any,
    distinct_tol = S7::class_numeric
  ))


#' @title Run an Optimiser From Many Starting Points
#'
#' @description
#' Wraps any optimiser and runs it from several starts, returning the best
#' result — and, as importantly, telling you how many different answers it found.
#'
#' @param optimizer The optimiser to run. Any of them, including another
#'   \code{multistart()}.
#' @param n How many starts, the user's own \code{par} among them. Defaults
#'   to \code{10}.
#' @param starts An optional matrix of starting points, one per row, used
#'   verbatim; \code{n} and \code{spread} are then ignored.
#' @param spread How widely the random starts are scattered, in units of the
#'   unconstrained scale. Defaults to \code{1}.
#' @param cluster An optional cluster from \pkg{parallel}. Defaults to
#'   \code{NULL}, meaning run the starts one after another; see Details.
#' @param distinct_tol Objective values differing by less than this are counted
#'   as the same optimum. Defaults to \code{1e-6}.
#' @param verbose Report each start as it finishes? Defaults to \code{FALSE}.
#' @param refresh Report every this many starts. Defaults to \code{1}.
#' @param keep_trace Keep the per-start summary? Defaults to \code{TRUE} — it is
#'   one row per start, not one per iteration, and it is the point of the method.
#'
#' @details
#' Every method in this package finds a \emph{local} minimum, and none of them
#' can tell you whether it is the global one; that is not a shortcoming of the
#' implementations but of local optimisation. Running from several starts is the
#' only general answer, and the most valuable thing it produces is not the best
#' point but the \strong{count of distinct optima}, which is direct evidence
#' about whether the question you asked has one answer. A run reporting one
#' optimum from twenty starts is worth far more than a single run reporting
#' convergence.
#'
#' The result is the best run, with everything it carried. The per-start summary
#' is in \code{trace}: one row each, with the value reached, whether that start
#' converged, and how many iterations it took. The message counts the starts that
#' succeeded, the ones that converged, and the distinct optima found.
#'
#' \subsection{Where the starts come from}{
#' The first is always \code{par}: your guess is a hypothesis worth testing, and
#' silently discarding it would be rude.
#'
#' The rest are a Latin hypercube — each coordinate's range is cut into \code{n}
#' equal strata and each is used exactly once, so the starts cannot all cluster
#' in one corner the way independent draws can. They are generated on the
#' \emph{unconstrained} scale and mapped back, which is what makes bounds
#' automatic: a start for a variance is drawn as a log and comes back positive,
#' and no draw is ever rejected for being outside the box.
#' }
#'
#' \subsection{Parallelism, and why it is at the level of processes}{
#' Pass a cluster made by \code{parallel::makeCluster()} and the starts are
#' distributed across it. \pkg{optimizers7} is loaded on the workers for you; a
#' cluster you made is a cluster you keep, so stopping it is yours to do.
#'
#' Threading the starts inside C++ would be faster and is not possible, and the
#' reason is worth stating because it is not the obvious one. It is easy to
#' assume that a \code{\link{cpp_objective}} would make the run callback-free and
#' so safe to thread. It would not: \strong{the stopping rule is an R object},
#' by design, and it is consulted at every iteration. That is the feature this
#' package was built around — a criterion the user can write — and it means every
#' run returns to R whatever the objective is made of. R is single-threaded, so
#' calling it from several threads is undefined behaviour that crashes rather
#' than errors.
#'
#' Separate processes have no such problem, so that is the route taken, and it
#' has the advantage of working for \emph{every} objective rather than only for
#' compiled ones.
#'
#' For reproducibility across workers use
#' \code{parallel::clusterSetRNGStream()}; ordinary \code{set.seed()} governs
#' only the sequential path. The \code{seed} recorded in the result is the
#' master's state, and reproduces a \emph{sequential} run only: the workers
#' draw from streams of their own, which the master never sees. Reproducing a
#' cluster run means recording the argument you gave
#' \code{clusterSetRNGStream()}, and there is no way for this function to do
#' that for you.
#' }
#'
#' \subsection{A start that fails is not a run that fails}{
#' A random start can easily land where the objective is undefined. Such a start
#' is recorded as failed and the others carry on; only if \emph{every} start
#' fails is that an error. Otherwise one bad draw would throw away nineteen good
#' answers.
#' }
#'
#' @return An S7 object of class \code{MultiStart}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @examples
#' multistart(bfgs())
#'
#' # a surface with two minima, one of them better
#' f <- function(p) (p[1]^2 - 1)^2 + p[2]^2 + 0.3 * p[1]
#' set.seed(1)
#' r <- minimize(multistart(bfgs(), n = 12), f, c(0, 0))
#' r@message
#' table(round(r@trace$value, 6))
#'
#' @seealso \code{\link{minimize}}, \code{\link{bfgs}}
#' @export
multistart <- function(optimizer, n = 10, starts = NULL, spread = 1,
                       cluster = NULL, distinct_tol = 1e-6,
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
  # Its budgets are the inner optimiser's, except maxit, which counts starts.
  check_optimizer_args(optimizer@criterion, n, optimizer@max_eval, verbose,
                       refresh, keep_trace)
  check_tol(spread)
  check_tol(distinct_tol)

  MultiStart(
    name = paste0("multistart (", optimizer@name, ")"),
    criterion = optimizer@criterion,
    maxit = n, max_eval = optimizer@max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    optimizer = optimizer, n = n, starts = starts, spread = spread,
    cluster = cluster, distinct_tol = distinct_tol
  )
}

#' @title What Multi-Start Can Offer a Stopping Rule
#' @name optimizer_provides.MultiStart
#' @description
#' Whatever the inner optimiser offers, since it is the inner optimiser that
#' evaluates the rule.
#' @param optimizer A \code{MultiStart} object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, MultiStart) <- function(optimizer)
  optimizer_provides(optimizer@optimizer)


# Changing a setting has to reach the optimiser that will actually use it. The
# criterion is the one that matters: the outer copy exists so that printing a
# multistart tells the truth, but the rule that is evaluated belongs to the
# optimiser inside, so setting only the outer one would change nothing while
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
#' Scatters \code{n - 1} starting points around \code{par}, generated on the
#' unconstrained scale so that any bounds are respected by construction.
#'
#' @details
#' Each coordinate's range is cut into \code{n - 1} strata, each used exactly
#' once, so the starts cannot all fall in one corner the way independent draws
#' can. The range is \code{par} plus or minus \code{3 * spread} on the
#' unconstrained scale, which for an unbounded parameter is the parameter itself
#' and for a bounded one is its log or logit — so a start for a variance is drawn
#' as a log and comes back positive without a single rejected draw.
#'
#' @param par The user's starting value.
#' @param n Total number of starts, including \code{par}.
#' @param spread Half-width of the sampling range, in unconstrained units.
#' @param bounds Box constraints in the shape \code{\link{check_bounds}}
#'   returns, possibly empty.
#'
#' @return A matrix with \code{n} rows, the first of which is \code{par}.
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


#' @title Minimise From Many Starting Points
#' @name minimize.MultiStart
#' @description Runs \code{\link{multistart}} on the objective.
#' @param optimizer A \code{MultiStart} object.
#' @param fn,par,gr,he,bounds,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}: the best run, with the per-start
#'   summary in its \code{trace}.
#' @keywords internal
S7::method(minimize, MultiStart) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    inner <- optimizer@optimizer
    # Validate once, here, so that a bad criterion or a start outside its bounds
    # is refused before n runs are launched rather than n times over.
    prepare_objective(inner, fn, par, gr, he)
    bx <- check_bounds(bounds, par)

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
    one <- function(i) {
      tryCatch(minimize(inner, fn, S[i, ], gr = gr, he = he, bounds = bounds),
               error = function(e) conditionMessage(e))
    }

    t0 <- proc.time()[["elapsed"]]
    cl <- optimizer@cluster
    if (is.null(cl)) {
      res <- vector("list", n)
      for (i in seq_len(n)) {
        res[[i]] <- one(i)
        if (optimizer@verbose && optimizer@refresh > 0 &&
            i %% optimizer@refresh == 0) {
          v <- if (is.character(res[[i]])) "failed"
               else format(res[[i]]@value, digits = 8)
          cat(sprintf("  start %4d of %d : %s\n", i, n, v))
        }
      }
    } else {
      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Running on a cluster needs the 'parallel' package.",
             call. = FALSE)
      }
      # A worker is a fresh session and knows nothing, so the package has to be
      # loaded there. Checking rather than assuming: without this the failure
      # arrives as checkForRemoteErrors() reporting "there is no package called
      # 'optimizers7'" from inside parLapply, which says nothing about what the
      # caller should do. It is a real case, not a hypothetical one -- a package
      # loaded with pkgload during development is not installed anywhere the
      # workers can see.
      have <- try(parallel::clusterEvalQ(
        cl, requireNamespace("optimizers7", quietly = TRUE)), silent = TRUE)
      if (inherits(have, "try-error") || !all(unlist(have))) {
        stop("optimizers7 must be installed where the cluster workers can ",
             "load it.\n  They are separate R sessions and do not inherit ",
             "this one's loaded packages;\n  a package loaded with pkgload, ",
             "or installed in a library the workers\n  do not have on their ",
             ".libPaths(), will not be found.", call. = FALSE)
      }
      res <- parallel::parLapply(cl, seq_len(n), one)
    }
    elapsed <- proc.time()[["elapsed"]] - t0

    build_multistart_result(res, S, optimizer, seed, elapsed)
  }


#' Record the Random Number Generator's State
#'
#' @description
#' The state \code{.Random.seed} held before a run that draws random numbers, so
#' that the run can be repeated exactly.
#'
#' @details
#' Assigning the recorded value back into the global environment reproduces the
#' run. When no stream exists yet one uniform is drawn to force \R to create
#' one, which costs a single discarded number on the first stochastic call of a
#' session.
#'
#' It is emphatically \strong{not} done with \code{set.seed(NULL)}, which
#' reseeds from the clock and throws away whatever the caller had set. That is a
#' defect this toolkit has already met once: it made a check in \pkg{distributions7}
#' silently random, and the resulting test failed about one CI run in several
#' hundred, on whichever platform happened to draw it.
#'
#' @return An integer vector, the value of \code{.Random.seed}.
#'
#' @keywords internal
capture_seed <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}


#' Assemble the Result of a Multi-Start Run
#'
#' @description
#' Picks the best run, and summarises what the others found.
#'
#' @details
#' The count of distinct optima is the reason to run this at all, so it is
#' computed rather than left to the caller: the values reached are sorted and
#' cut wherever consecutive ones differ by more than \code{distinct_tol}. It is
#' a statement about the objective, not about the optimiser, and it is the one
#' piece of evidence a single run can never supply.
#'
#' @param res The list of results, entries that failed being character messages.
#' @param S The matrix of starting points.
#' @param optimizer The \code{MultiStart} object.
#' @param seed The generator state the run began with.
#' @param elapsed Seconds.
#'
#' @return An \code{\link{optimizer_result}}.
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
    # One iteration of this optimiser is one start; the inner iteration counts
    # are in the trace, where they can be read per start rather than summed into
    # a number that means nothing.
    iterations = length(res),
    converged = out@converged,
    criterion_met = out@criterion_met,
    message = msg, trace = trace, optimizer = optimizer, elapsed = elapsed,
    seed = seed
  )
}
