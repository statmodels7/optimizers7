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
#' Beyond the seven properties every optimizer has, a `NelderMead` carries
#' five of its own: `step` and `simplex` describe where the run begins,
#' `adaptive` which coefficients it uses, and `max_restarts` with
#' `degenerate_tol` the safeguard.
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
#'   `crit_stationary(1e-8)`, read on the diameter of the simplex, so the
#'   tolerance is on the parameter scale. A gradient rule is refused when the
#'   run starts.
#' @param step Size of the initial simplex, relative to each coordinate of
#'   the starting value. A single positive number, default `0.1`.
#' @param adaptive Use Gao and Han's dimension-dependent coefficients?
#'   `TRUE` or `FALSE`, default `TRUE`; see below. At \eqn{p = 2} the two
#'   settings are identical.
#' @param max_restarts How many times a collapsed simplex may be rebuilt. A
#'   single non-negative number, default `3`; `0` turns the safeguard off.
#' @param degenerate_tol Rebuild the simplex when its conditioning falls
#'   below this, a single positive number, default `1e-6`. A non-positive
#'   value is refused, with a message naming `tol` rather than the argument.
#' @param simplex An optional starting simplex: a matrix with \eqn{p+1} rows,
#'   one vertex per row. Defaults to `NULL`, meaning build one from `par` and
#'   `step`. Anything that is not a matrix or `NULL` is refused.
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
#' assumed to exist, which is why the method survives a kink. It is also why
#' it is slow: an ordering of \eqn{p+1} values is very little information
#' about a surface, and on the quadratic below it costs 133 evaluations where
#' [bfgs()] costs 3.
#'
#' # Degenerate simplices
#'
#' Nelder-Mead can converge to a point that is not a minimizer. McKinnon
#' (1998) exhibited a strictly convex function with continuous derivatives on
#' which it performs inside contractions for ever: the simplex flattens onto
#' a line through a point where the gradient is not zero, every vertex
#' agrees, and every ordinary stopping rule reports success. The *values*
#' offer no defense, behaving exactly as convergence would; what has gone
#' wrong is the *shape*.
#'
#' So the shape is watched. The conditioning measured is
#' \eqn{\lvert \det E \rvert} divided by the product of the edge lengths,
#' where \eqn{E} holds the edges from the best vertex. It is 1 for a
#' right-angled simplex, 0 for one that has collapsed into a lower dimension,
#' and unchanged by rescaling, so one threshold serves at every size. When it
#' falls below `degenerate_tol` the simplex is rebuilt right-angled at the
#' current best vertex, at the diameter it had reached, so the scale the run
#' has earned is kept. Each rebuild is counted, appears in the trace as
#' `restart` and is reported in the result's `message`.
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
#' at no cost to the small problems anyone would recognize. At \eqn{p = 10}
#' they are \eqn{1, 1.2, 0.7, 0.9}, and on a quadratic in ten unknowns the
#' adaptive run takes 959 iterations against 1230.
#'
#' # When to reach for it
#'
#' Rarely, and knowingly. On a smooth objective every gradient-based method
#' here beats it by orders of magnitude. Its place is an objective that is
#' genuinely non-smooth or noisy and whose subgradients are unavailable.
#' Where subgradients *are* available, [bundle()] uses them and this method
#' uses nothing.
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
#' # A smooth problem, to show that it works and what it costs. BFGS reaches
#' # the same point in 3 evaluations.
#' q <- function(p) sum((p - c(1, 2))^2)
#' minimize(nelder_mead(), q, c(0, 0))@counts[["f"]]
#' minimize(bfgs(), q, c(0, 0), gr = function(p) 2 * (p - c(1, 2)))@counts[["f"]]
#'
#' # What it is for: a sum of absolute deviations, whose minimizer is the
#' # median and whose derivative does not exist there.
#' set.seed(1)
#' y <- rcauchy(101)
#' r <- minimize(nelder_mead(), function(p) sum(abs(y - p)), par = 0)
#' c(nelder_mead = r@par, median = median(y))
#' abs(r@par - median(y))
#'
#' # The adaptive coefficients are the classical ones at p = 2 and differ
#' # above it, where they save iterations.
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
#'   [crit_stationary()] for the rule this method reads.
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
  check_tol(degenerate_tol)
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
#' Beyond the seven properties every optimizer has, a `Compass` carries five
#' of its own: `step`, `expand` and `shrink` govern the poll radius,
#' `directions` says which set is polled, and `opportunistic` when the poll
#' stops.
#'
#' @param step Initial poll size, relative to the starting value.
#' @param directions Either `"mads"` or `"coordinate"`.
#' @param opportunistic Logical; whether the poll stops at the first
#'   improvement it finds.
#' @param expand,shrink Factors applied to the poll size after a success and
#'   after a failure.
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
#' Looks around the current point along a set of directions; moves to the first
#' or best improvement it finds, and shrinks the radius when it finds none.
#' Uses no derivative, and unlike [nelder_mead()] it comes with a
#' convergence theorem.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_stationary(1e-8)`, read on the poll size \eqn{\Delta}. This is the
#'   rule with a theorem behind it: the limit points of a pattern search with
#'   \eqn{\Delta \to 0} are Clarke stationary.
#' @param step Initial poll size, scaled by the largest coordinate of the
#'   starting value. A single positive number, default `0.1`.
#' @param directions `"mads"` (default) or `"coordinate"`. Partial matching
#'   applies and any other string is refused, naming both.
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
#' is a non-negative combination of them. That is what buys the theorem. If
#' the current point is not stationary then some direction in the set goes
#' downhill, so a poll that fails everywhere is evidence about the point and
#' not about the directions. A failed poll therefore licenses shrinking the
#' radius, and the limit points of a run whose radius goes to zero are
#' stationary.
#'
#' # Poll directions
#'
#' \describe{
#'   \item{`"coordinate"`}{the \eqn{2p} signed axes, which is compass search.
#'     Cheap, deterministic and reproducible without a seed. The theorem it
#'     enjoys assumes \eqn{f} is continuously differentiable.}
#'   \item{`"mads"`}{a fresh random orthonormal basis at every poll, taken
#'     plus and minus.}
#' }
#'
#' The difference decides the outcome on the problems this method exists for.
#' Where \eqn{f} is merely Lipschitz, a *fixed* set of directions can fail: a
#' kink whose ridge runs diagonally is descended by no coordinate direction,
#' so the poll fails at a point that is not stationary and the run stops
#' there. Measured on \eqn{\lvert x_1 + x_2 \rvert + 0.1\lVert x \rVert^2},
#' whose minimum is 0, coordinate polling from \eqn{(1, 0.5)} stops at
#' `0.05` on every one of ten runs, the same value each time; the random
#' poll has a median of `0.014` and reaches `9.8e-05` at best.
#'
#' The repair, which is the idea behind MADS, is to let the directions used
#' over the whole run become dense in the sphere, so no direction of descent
#' is missed for ever, and drawing a new orthonormal basis at each poll
#' achieves that with probability one. What is implemented is that idea
#' rather than LTMADS as published, and it is the property the convergence
#' proof rests on. The cost is reproducibility: a random poll draws from \R's
#' generator, so `set.seed()` governs the run and the state is recorded in
#' the result's `seed`.
#'
#' # Opportunistic polling
#'
#' Accepting the first improvement instead of the best costs a worse
#' direction and saves up to \eqn{2p - 1} evaluations per poll. Measured on
#' the quadratic below: 332 evaluations against 393, in 96 iterations against
#' 98. The saving in evaluations is real and the cost in iterations was not
#' visible here, so `TRUE` is the default.
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
#' # A kink running diagonally, which no coordinate direction descends. The
#' # true minimum is 0. Coordinate polling stops at the same wrong point every
#' # time; the random poll gets past it.
#' f <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
#' minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#' minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#' set.seed(1)
#' minimize(compass(directions = "mads"), f, c(1, 0.5))@value
#'
#' # A coordinate poll needs no seed, a random one does, and it records the
#' # state it used.
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
#'   available, [crit_stationary()] for the rule this method reads.
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
#' Both derivative-free methods report `"stationarity"` and no gradient,
#' having none to report. A rule reading a gradient is refused when the run
#' starts, so it cannot sit there testing `NULL` at every iteration and never
#' firing.
#'
#' @details
#' The measure differs by method and [crit_stationary()] reads whichever is
#' offered: for [nelder_mead()] the **diameter of the simplex**, so the
#' tolerance is on the parameter scale, and for [compass()] the **poll size**
#' \eqn{\Delta}, the quantity Torczon's theorem is stated in.
#'
#' @param optimizer A `NelderMead` or `Compass` object.
#'
#' @return The character vector `"stationarity"`.
#'
#' @examples
#' optimizer_provides(nelder_mead())
#' optimizer_provides(compass())
#'
#' # So a gradient rule is refused, naming the method.
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
#' it until its diameter falls below the tolerance.
#'
#' @details
#' A `simplex` supplied on the optimizer is checked here against the
#' starting value, since only now is the number of parameters known: it must
#' have \eqn{p+1} rows and \eqn{p} columns, and a mismatch reports both
#' shapes.
#'
#' @param optimizer A `NelderMead` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` and `he` are
#'   accepted and ignored, the method using no derivative; refusing them
#'   would force calling code to branch on the algorithm. Bounds are taken
#'   and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is `NULL` and whose trace
#'   carries a `stationarity` column holding the simplex diameter. The
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
#' expand, and shrink the radius when the poll fails everywhere.
#'
#' @param optimizer A `Compass` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` and `he` are
#'   accepted and ignored, the method using no derivative. Bounds are taken
#'   and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is `NULL` and whose trace
#'   carries a `stationarity` column holding the poll size. With
#'   `directions = "mads"` the `seed` is recorded, the poll drawing from \R's
#'   generator; with `"coordinate"` it is `NULL` and the run repeats without
#'   one.
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
