#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

# The two methods that use no derivative at all. Both report a stationarity
# measure in place of a gradient norm; see crit_stationary().


#' @title S7 Class for Nelder-Mead
#'
#' @description
#' An optimizer holding the size and shape of the initial simplex, which set
#' of reflection coefficients to use, and the safeguard against a simplex
#' that has collapsed. Built by [nelder_mead()]. The simplex itself moves
#' inside the run; the `simplex` property is a starting one, or `NULL`.
#'
#' @details
#' Beyond the seven properties of every optimizer, a `NelderMead` carries
#' five of its own: `step` and `simplex` describe the starting simplex,
#' `adaptive` selects the coefficients, and `max_restarts` and
#' `degenerate_tol` control the safeguard.
#'
#' @param step Relative size of the initial simplex.
#' @param adaptive Logical; whether the dimension-dependent coefficients of
#'   Gao and Han are used.
#' @param max_restarts How many times a degenerate simplex may be rebuilt.
#' @param degenerate_tol The conditioning below which it is rebuilt.
#' @param simplex An optional starting simplex, a matrix with one vertex per
#'   row, or `NULL`.
#'
#' @return An S7 object of class `NelderMead` inheriting from [optimizer()],
#'   with the five properties above beside the seven shared ones.
#'
#' @seealso [nelder_mead()] for the constructor, [Compass] for the other
#'   derivative-free method.
#' @name NelderMead-class
#' @aliases NelderMead
#' @keywords internal
NelderMead <- S7::new_class("NelderMead", parent = optimizer,
  properties = list(
    step           = S7::class_numeric,
    adaptive       = S7::class_logical,
    max_restarts   = S7::class_numeric,
    degenerate_tol = S7::class_numeric,
    simplex        = S7::class_any
  ))


#' @title The Nelder-Mead Simplex Method
#'
#' @description
#' The Nelder-Mead simplex method: \eqn{p+1} points are maintained and the
#' worst is reflected, expanded, contracted or shrunk through the centroid of
#' the others. Only objective values are used; no derivative is required.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_stationary(1e-8)`, read on the size of the simplex (the largest
#'   distance from the best vertex to another vertex), so the tolerance is on
#'   the parameter scale. A gradient rule is rejected when the run starts.
#' @param step Size of the initial simplex: the vertex that moves coordinate
#'   \eqn{j} is displaced by `step` times the larger of 1 and the absolute
#'   value of that coordinate of the starting value. A single positive
#'   number, default `0.1`.
#' @param adaptive Use Gao and Han's dimension-dependent coefficients?
#'   `TRUE` or `FALSE`, default `TRUE`; see below. At \eqn{p = 2} the two
#'   settings are identical.
#' @param max_restarts How many times a collapsed simplex may be rebuilt. A
#'   single non-negative number, default `3`; `0` turns the safeguard off.
#' @param degenerate_tol Rebuild the simplex when its conditioning falls
#'   below this, a single positive number, default `1e-6`. A value that is
#'   not positive is rejected.
#' @param simplex An optional starting simplex: a matrix with \eqn{p+1} rows,
#'   one vertex per row. Defaults to `NULL`, meaning build one from `par` and
#'   `step`. A value that is not a matrix or `NULL` is rejected.
#' @param maxit Maximum iterations. Defaults to 2000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 50.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' Reflect the worst vertex through the centroid of the rest; if that is the
#' best point yet, try going further; if it is no better than the second
#' worst, pull back; and if even that fails, shrink everything towards the
#' best vertex. No model of the function is built and no derivative is
#' assumed to exist, so the method works at a kink. It is slow because the
#' ordering of \eqn{p+1} values carries little information about the
#' surface; the first example compares its evaluation count with that of
#' [bfgs()].
#'
#' # Degenerate simplices
#'
#' Nelder-Mead can converge to a point that is not a minimizer. McKinnon
#' (1998) exhibited a strictly convex function with continuous derivatives on
#' which it performs inside contractions indefinitely: the simplex flattens
#' onto a line through a point where the gradient is not zero, the vertices
#' approach a common value, and every ordinary stopping rule reports success.
#' The function values behave as they would at convergence, and only the
#' shape of the simplex shows the failure.
#'
#' The method therefore monitors the shape. The conditioning measured is
#' \eqn{\lvert \det E \rvert} divided by the product of the edge lengths,
#' where \eqn{E} holds the edges from the best vertex. It is 1 for a
#' right-angled simplex, 0 for one that has collapsed into a lower dimension,
#' and unchanged by rescaling, so one threshold serves at every size. When it
#' falls below `degenerate_tol` the simplex is rebuilt with right angles at
#' the current best vertex, with edges equal to the size that the simplex had
#' reached. Each rebuild is counted, appears in the trace as `restart` and is
#' reported in the `message` of the result.
#'
#' The safeguard is not free: a restart costs \eqn{p+1} evaluations and can
#' delay a genuine convergence. `max_restarts = 0` turns it off.
#'
#' # Adaptive coefficients
#'
#' The classical reflection, expansion, contraction and shrink factors are
#' \eqn{1, 2, 1/2, 1/2}, chosen when the method was proposed for two or three
#' parameters. In higher dimension a fixed expansion of 2 makes the simplex
#' overshoot along whichever direction it happened to try. Gao and Han (2012)
#' replace them by \eqn{1,\ 1 + 2/p,\ 3/4 - 1/(2p),\ 1 - 1/p}, which at
#' \eqn{p = 2} are *exactly* the classical values, so `TRUE` is the default
#' without changing the method on problems with two parameters. At
#' \eqn{p = 10} they are \eqn{1, 1.2, 0.7, 0.9}, and on a quadratic in ten
#' unknowns the adaptive run needs fewer iterations than the classical one
#' (see the examples).
#'
#' # Where the method applies
#'
#' On a smooth objective the gradient-based methods that use curvature, such
#' as [newton()], [bfgs()] and [lbfgs()], need far fewer evaluations. The
#' method is intended for an objective that is non-smooth or noisy and whose
#' subgradients are unavailable. Where subgradients *are* available,
#' [bundle()] uses them.
#'
#' @return An S7 object of class [NelderMead], inheriting from [optimizer()],
#'   to be handed to [minimize()].
#'
#' @references
#' Nelder, J. A. and Mead, R. (1965). A simplex method for function
#' minimization. *The Computer Journal* **7**, 308--313.
#'
#' McKinnon, K. I. M. (1998). Convergence of the Nelder-Mead simplex method to a
#' nonstationary point. *SIAM Journal on Optimization* **9**, 148--158.
#'
#' Gao, F. and Han, L. (2012). Implementing the Nelder-Mead simplex algorithm
#' with adaptive parameters. *Computational Optimization and Applications*
#' **51**, 259--277.
#'
#' @examples
#' nelder_mead()
#'
#' # A smooth problem, to show the cost: BFGS reaches the same point in far
#' # fewer evaluations.
#' q <- function(p) sum((p - c(1, 2))^2)
#' minimize(nelder_mead(), q, c(0, 0))@counts[["f"]]
#' minimize(bfgs(), q, c(0, 0), gr = function(p) 2 * (p - c(1, 2)))@counts[["f"]]
#'
#' # A sum of absolute deviations, whose minimizer is the median and whose
#' # derivative does not exist there.
#' set.seed(1)
#' y <- rcauchy(101)
#' r <- minimize(nelder_mead(), function(p) sum(abs(y - p)), par = 0)
#' c(nelder_mead = r@par, median = median(y))
#' abs(r@par - median(y))
#'
#' # The adaptive coefficients equal the classical ones at p = 2 and differ
#' # above it; here they save iterations.
#' f10 <- function(x) sum((x - 1:10)^2)
#' c(adaptive = minimize(nelder_mead(maxit = 20000), f10, numeric(10))@iterations,
#'   classical = minimize(nelder_mead(adaptive = FALSE, maxit = 20000),
#'                        f10, numeric(10))@iterations)
#'
#' # The trace names the simplex operation taken at each iteration.
#' unique(minimize(nelder_mead(keep_trace = TRUE), q, c(0, 0))@trace$safeguard)
#'
#' @seealso [compass()] for the derivative-free method with a convergence
#'   theorem, [bundle()] when subgradients are available,
#'   [crit_stationary()] for the default rule of this method.
#' @export
nelder_mead <- function(criterion = crit_stationary(),
                        step = 0.1, adaptive = TRUE,
                        max_restarts = 3, degenerate_tol = 1e-6,
                        simplex = NULL,
                        maxit = 2000, max_eval = Inf,
                        verbose = FALSE, refresh = 50, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  if (length(adaptive) != 1L || !is.logical(adaptive) || is.na(adaptive)) {
    stop("'adaptive' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(max_restarts) != 1L || !is.numeric(max_restarts) ||
      is.na(max_restarts) || max_restarts < 0) {
    stop("'max_restarts' must be a single non-negative number.", call. = FALSE)
  }
  check_tol(degenerate_tol, "degenerate_tol")
  if (!is.null(simplex) && !is.matrix(simplex)) {
    stop("'simplex' must be a matrix with one vertex per row, or NULL.",
         call. = FALSE)
  }
  NelderMead(
    name = "nelder-mead", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    step = step, adaptive = adaptive, max_restarts = max_restarts,
    degenerate_tol = degenerate_tol, simplex = simplex
  )
}


#' @title S7 Class for Pattern Search
#'
#' @description
#' An optimizer holding the poll size and how it changes, which set of poll
#' directions is used, and whether the poll stops at the first improvement.
#' Built by [compass()]. With `directions = "mads"` a run draws from \R's
#' generator and records its seed; with `"coordinate"` it is deterministic
#' and records none.
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Compass` carries five
#' of its own: `step`, `expand` and `shrink` govern the poll radius,
#' `directions` selects the set that is polled, and `opportunistic` decides
#' when the poll stops.
#'
#' @param step Initial poll size, relative to the starting value.
#' @param directions Either `"mads"` or `"coordinate"`.
#' @param opportunistic Logical; whether the poll stops at the first
#'   improvement it finds.
#' @param expand,shrink Factors applied to the poll size after a success and
#'   after a failure. After a success the poll size is capped at \eqn{10^4}
#'   times its initial value.
#'
#' @return An S7 object of class `Compass` inheriting from [optimizer()],
#'   with the five properties above beside the seven shared ones.
#'
#' @seealso [compass()] for the constructor, [NelderMead] for the other
#'   derivative-free method.
#' @name Compass-class
#' @aliases Compass
#' @keywords internal
Compass <- S7::new_class("Compass", parent = optimizer,
  properties = list(
    step          = S7::class_numeric,
    directions    = S7::class_character,
    opportunistic = S7::class_logical,
    expand        = S7::class_numeric,
    shrink        = S7::class_numeric
  ))


#' @title Pattern Search, With Coordinate or Random Polling
#'
#' @description
#' Evaluates the objective around the current point along a set of
#' directions, moves to the first or to the best improvement, and shrinks
#' the radius when it finds none. It uses no derivative, and, unlike
#' [nelder_mead()], it has a convergence theorem.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_stationary(1e-8)`, read on the poll size \eqn{\Delta}, the
#'   quantity in which the convergence theory is stated: for a continuously
#'   differentiable objective, the limit points of a coordinate pattern
#'   search with \eqn{\Delta \to 0} are stationary (Torczon 1997), and for a
#'   locally Lipschitz objective the limit points of MADS are Clarke
#'   stationary (Audet and Dennis 2006).
#' @param step Initial poll size: `step` times the larger of 1 and the
#'   largest absolute coordinate of the starting value. A single positive
#'   number, default `0.1`.
#' @param directions `"mads"` (default) or `"coordinate"`. Partial matching
#'   applies, and any other string is rejected by [match.arg()].
#' @param opportunistic Move to the first improvement found instead of
#'   polling every direction? `TRUE` or `FALSE`, default `TRUE`.
#' @param expand Factor applied to the poll size after a success. A single
#'   number of at least 1, default `2`.
#' @param shrink Factor applied after a failure. A single number strictly
#'   inside \eqn{(0, 1)}, default `0.5`.
#' @param maxit Maximum iterations. Defaults to 2000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 50.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' The directions form a *positive spanning set*: every vector in the space
#' is a non-negative combination of them. The theorem rests on this. If the
#' current point is not stationary, some direction in the set goes downhill,
#' so a poll that fails in every direction at a small radius indicates that
#' the point is close to stationary. A failed poll therefore justifies
#' shrinking the radius, and the limit points of a run whose radius goes to
#' zero are stationary.
#'
#' # Poll directions
#'
#' \describe{
#'   \item{`"coordinate"`}{the \eqn{2p} signed axes, which is compass search.
#'     This poll is cheap, deterministic and reproducible without a seed. Its
#'     convergence theorem assumes that \eqn{f} is continuously
#'     differentiable.}
#'   \item{`"mads"`}{a fresh random orthonormal basis at every poll, taken
#'     plus and minus.}
#' }
#'
#' The choice decides the outcome on the problems for which the method is
#' intended. Where \eqn{f} is only Lipschitz, a *fixed* set of directions can
#' fail: if a kink has a diagonal ridge, no coordinate direction descends
#' it, so the poll fails at a point that is not stationary and the run stops
#' there. On \eqn{\lvert x_1 + x_2 \rvert + 0.1\lVert x \rVert^2}, whose
#' minimum is 0, coordinate polling from \eqn{(1, 0.5)} stops at `0.05`,
#' whereas the random poll reaches values several orders of magnitude lower,
#' which depend on the seed (see the examples).
#'
#' The remedy, which is the idea behind MADS, is to let the directions used
#' over the whole run become dense in the sphere, so that every direction of
#' descent is eventually polled; drawing a new orthonormal basis at each poll
#' achieves this with probability one. The implementation follows this idea
#' and is not LTMADS as published. The cost is reproducibility: a random poll
#' draws from the generator of \R, so `set.seed()` governs the run and the
#' state is recorded in the `seed` of the result.
#'
#' The poll size is capped at \eqn{10^4} times its initial value, and the
#' trace records `poll capped` when the cap applies. A run whose poll size
#' falls below \eqn{10^{-300}} stops with `converged = FALSE` and the
#' message `the poll size underflowed`.
#'
#' # Opportunistic polling
#'
#' Accepting the first improvement instead of the best gives a possibly
#' worse direction and saves up to \eqn{2p - 1} evaluations per poll. The
#' default is `TRUE`; the last example compares the evaluation counts of the
#' two settings.
#'
#' @return An S7 object of class [Compass], inheriting from [optimizer()], to
#'   be handed to [minimize()].
#'
#' @references
#' Torczon, V. (1997). On the convergence of pattern search algorithms.
#' *SIAM Journal on Optimization* **7**, 1--25.
#'
#' Audet, C. and Dennis, J. E. (2006). Mesh adaptive direct search algorithms
#' for constrained optimization. *SIAM Journal on Optimization*
#' **17**, 188--217.
#'
#' @examples
#' compass()
#'
#' # A kink along a diagonal, which no coordinate direction descends. The
#' # true minimum is 0. Coordinate polling stops at the same wrong point every
#' # time; the random poll passes it.
#' f <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
#' minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#' minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#' set.seed(1)
#' minimize(compass(directions = "mads"), f, c(1, 0.5))@value
#'
#' # A coordinate poll needs no seed; a random poll needs one, and records
#' # the state that it used.
#' a <- minimize(compass(directions = "coordinate"), f, c(1, 0.5))
#' b <- minimize(compass(directions = "coordinate"), f, c(1, 0.5))
#' c(identical(a@par, b@par), is.null(a@seed))
#'
#' set.seed(7); m1 <- minimize(compass(), f, c(1, 0.5))
#' set.seed(7); m2 <- minimize(compass(), f, c(1, 0.5))
#' c(identical(m1@par, m2@par), length(m1@seed) > 0)
#'
#' # Polling every direction costs more evaluations for the same answer.
#' q <- function(p) sum((p - c(1, 2))^2)
#' set.seed(5); first <- minimize(compass(), q, c(0, 0))
#' set.seed(5); every <- minimize(compass(opportunistic = FALSE), q, c(0, 0))
#' c(first = first@counts[["f"]], every = every@counts[["f"]])
#'
#' @seealso [nelder_mead()] for the simplex, [bundle()] when subgradients are
#'   available, [crit_stationary()] for the default rule of this method.
#' @export
compass <- function(criterion = crit_stationary(),
                    step = 0.1,
                    directions = c("mads", "coordinate"),
                    opportunistic = TRUE, expand = 2, shrink = 0.5,
                    maxit = 2000, max_eval = Inf,
                    verbose = FALSE, refresh = 50, keep_trace = FALSE) {
  directions <- match.arg(directions)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  if (length(opportunistic) != 1L || !is.logical(opportunistic) ||
      is.na(opportunistic)) {
    stop("'opportunistic' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(expand) != 1L || !is.numeric(expand) || is.na(expand) ||
      expand < 1) {
    stop("'expand' must be a single number of at least 1.", call. = FALSE)
  }
  if (length(shrink) != 1L || !is.numeric(shrink) || is.na(shrink) ||
      shrink <= 0 || shrink >= 1) {
    stop("'shrink' must be a single number in (0, 1).", call. = FALSE)
  }
  Compass(
    name = paste0("pattern search (", directions, ")"), criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    step = step, directions = directions, opportunistic = opportunistic,
    expand = expand, shrink = shrink
  )
}


#' @title What a Derivative-Free Method Can Offer a Stopping Rule
#' @name optimizer_provides.NelderMead
#'
#' @description
#' [nelder_mead()], [compass()] and [brent()] report `"stationarity"` and no
#' gradient, since they compute none. A rule that reads a gradient is
#' rejected when the run starts; without this check it would test `NULL` at
#' every iteration and never fire.
#'
#' @details
#' The measure differs by method, and [crit_stationary()] reads the measure
#' of the method: for [nelder_mead()] the **size of the simplex**, the
#' largest distance from the best vertex to another vertex, so the tolerance
#' is on the parameter scale; for [compass()] the **poll size**
#' \eqn{\Delta}, the quantity in which Torczon's theorem is stated; and for
#' [brent()] the half-width of the bracket.
#'
#' @param optimizer A `NelderMead`, `Compass` or `Brent` object.
#'
#' @return The character vector `"stationarity"`.
#'
#' @examples
#' optimizer_provides(nelder_mead())
#' optimizer_provides(compass())
#'
#' # A gradient rule is therefore rejected, and the error names the method.
#' try(minimize(nelder_mead(criterion = crit_grad()),
#'              function(p) sum(p^2), c(1, 1)))
#'
#' @keywords internal
S7::method(optimizer_provides, NelderMead) <- function(optimizer)
  "stationarity"

#' @rdname optimizer_provides.NelderMead
#' @name optimizer_provides.Compass
#' @keywords internal
S7::method(optimizer_provides, Compass) <- function(optimizer)
  "stationarity"


#' @title Minimize by Nelder-Mead
#' @name minimize.NelderMead
#'
#' @description
#' Runs [nelder_mead()] on the objective: build a simplex of \eqn{p+1}
#' vertices from `par` and `step`, then reflect, expand, contract or shrink
#' it until the stopping rule fires.
#'
#' @details
#' A `simplex` supplied on the optimizer is checked here against the
#' starting value, since only now is the number of parameters known: it must
#' have \eqn{p+1} rows and \eqn{p} columns, and a mismatch is rejected with
#' an error that states the required numbers of rows and columns.
#'
#' @param optimizer A `NelderMead` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` and `he` are
#'   accepted and ignored, since the method uses no derivative; rejecting
#'   them would force calling code to branch on the algorithm. Bounds are
#'   taken and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is `NULL` and whose trace
#'   has a `stationarity` column holding the size of the simplex. The
#'   `safeguard` column names the operation taken at each iteration:
#'   `reflect`, `expand`, `contract in`, `contract out`, `shrink`, or
#'   `restart` when a degenerate simplex was rebuilt.
#'
#' @keywords internal
S7::method(minimize, NelderMead) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
    bounds <- check_bounds(lower, upper, par)

    sx <- optimizer@simplex
    if (is.null(sx)) {
      sx <- matrix(numeric(0), 0, 0)
    } else if (nrow(sx) != length(par) + 1L || ncol(sx) != length(par)) {
      stop("'simplex' must have ", length(par) + 1L, " rows and ",
           length(par), " columns for this starting value.", call. = FALSE)
    }

    t0 <- proc.time()[["elapsed"]]
    out <- nelder_mead_run(
      spec = spec, par = as.numeric(par),
      criterion = optimizer@criterion, crit_fn = crit_met,
      step = optimizer@step, adaptive = optimizer@adaptive,
      max_restarts = as.integer(optimizer@max_restarts),
      degenerate_tol = optimizer@degenerate_tol,
      start_simplex = sx,
      maxit = as.integer(optimizer@maxit),
      max_eval = budget_int(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0)
  }


#' @title Minimize by Pattern Search
#' @name minimize.Compass
#'
#' @description
#' Runs [compass()] on the objective: poll the directions around the current
#' point at the current radius, move to an improvement if one is found and
#' expand the radius, and shrink the radius when the poll fails in every
#' direction.
#'
#' @param optimizer A `Compass` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` and `he` are
#'   accepted and ignored, the method using no derivative. Bounds are taken
#'   and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is `NULL` and whose trace
#'   has a `stationarity` column holding the poll size. With
#'   `directions = "mads"` the `seed` is recorded, since the poll draws from
#'   the generator of \R; with `"coordinate"` it is `NULL` and the run
#'   repeats without a seed.
#'
#' @keywords internal
S7::method(minimize, Compass) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
    bounds <- check_bounds(lower, upper, par)

    # A mads poll draws a fresh basis at every iteration, so the run is
    # reproducible only if the state it started from is known. Recording it is
    # cheaper than discovering afterwards that the interesting run cannot be
    # repeated.
    mads <- identical(optimizer@directions, "mads")
    seed <- if (mads) capture_seed() else NULL

    t0 <- proc.time()[["elapsed"]]
    out <- compass_run(
      spec = spec, par = as.numeric(par),
      criterion = optimizer@criterion, crit_fn = crit_met,
      step = optimizer@step,
      random_directions = mads,
      opportunistic = optimizer@opportunistic,
      expand = optimizer@expand, shrink = optimizer@shrink,
      maxit = as.integer(optimizer@maxit),
      max_eval = budget_int(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0, seed)
  }
