#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

# The two methods that use no derivative at all. Both report a stationarity
# measure in place of a gradient norm; see crit_stationary().


#' @title S7 Class for Nelder-Mead
#' @description The class \code{\link{nelder_mead}} instantiates.
#' @param step Relative size of the initial simplex.
#' @param adaptive Whether to use dimension-dependent coefficients.
#' @param max_restarts How many times a degenerate simplex may be rebuilt.
#' @param degenerate_tol The conditioning below which it is rebuilt.
#' @param simplex An optional starting simplex.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{nelder_mead}}
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
#' Keeps \eqn{p+1} points, and at each step reflects the worst of them through
#' the face of the others. It uses no derivative and asks nothing of the
#' objective but that it return a number.
#'
#' @param criterion The stopping rule. Defaults to
#'   \code{crit_stationary(1e-8)}, on the diameter of the simplex.
#' @param step Size of the initial simplex, relative to each coordinate of the
#'   starting value. Defaults to \code{0.1}.
#' @param adaptive Use Gao and Han's dimension-dependent coefficients? Defaults
#'   to \code{TRUE}; see Details.
#' @param max_restarts How many times a collapsed simplex may be rebuilt.
#'   Defaults to \code{3}. Zero disables the safeguard.
#' @param degenerate_tol Rebuild the simplex when its conditioning falls below
#'   this. Defaults to \code{1e-6}.
#' @param simplex An optional starting simplex: a matrix with \eqn{p+1} rows,
#'   one vertex per row. Defaults to \code{NULL}, meaning build one from
#'   \code{par} and \code{step}.
#' @param maxit Maximum iterations. Defaults to 2000.
#' @param max_eval Maximum objective evaluations. Defaults to 20000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 50.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' Reflect the worst vertex through the centroid of the rest; if that is the
#' best point yet, try going further; if it is no better than the second worst,
#' pull back; and if even that fails, shrink everything towards the best vertex.
#' No model of the function is built and no derivative is assumed to exist,
#' which is why it survives a kink — and also why it is slow, since an ordering
#' of \eqn{p+1} values is very little information about a surface.
#'
#' \subsection{Degenerate simplices}{
#' Nelder-Mead can converge to a point that is not a minimiser. McKinnon (1998)
#' exhibited a strictly convex function with continuous derivatives on which it
#' performs inside contractions for ever: the simplex flattens onto a line
#' through a point where the gradient is not zero, every vertex agrees, and
#' every ordinary stopping rule reports success. There is no defence in the
#' \emph{values} — they behave exactly as convergence would — because what has
#' gone wrong is the \emph{shape} of the simplex.
#'
#' So the shape is what is watched. The conditioning measured is
#' \eqn{\lvert \det E \rvert} divided by the product of the edge lengths, where
#' \eqn{E} holds the edges from the best vertex: it is 1 for a right-angled
#' simplex, 0 for one that has collapsed into a lower dimension, and unchanged
#' by rescaling, so one threshold serves at every size. When it falls below
#' \code{degenerate_tol} the simplex is rebuilt right-angled at the current best
#' vertex, at the diameter it had reached — keeping the scale the run has earned
#' rather than starting over. Each rebuild is counted, reported in the trace as
#' \code{"restart"} and in the result's message.
#'
#' The safeguard is not free: a restart costs \eqn{p+1} evaluations and can
#' delay a genuine convergence. \code{max_restarts = 0} turns it off.
#' }
#'
#' \subsection{Adaptive coefficients}{
#' The classical reflection, expansion, contraction and shrink factors are
#' \eqn{1, 2, 1/2, 1/2}, chosen when the method was proposed for two or three
#' parameters. In higher dimension a fixed expansion of 2 makes the simplex
#' overshoot along whichever direction it happened to try. Gao and Han (2012)
#' replace them by \eqn{1,\ 1 + 2/p,\ 3/4 - 1/(2p),\ 1 - 1/p}, which at
#' \eqn{p = 2} reduce \emph{exactly} to the classical values — so the default is
#' \code{TRUE} at no cost to the small problems anyone would recognise.
#' }
#'
#' \subsection{Scope}{
#' Rarely, and knowingly. If the objective is smooth, every gradient-based
#' method here will beat it by orders of magnitude. Its place is an objective
#' that is genuinely non-smooth or noisy and whose subgradients are not
#' available; if they \emph{are} available, \code{\link{bundle}} is the better
#' tool, since it uses them and this uses nothing.
#' }
#'
#' @return An S7 object of class \code{NelderMead}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Nelder, J. A. and Mead, R. (1965). A simplex method for function
#' minimization. \emph{The Computer Journal} \strong{7}, 308--313.
#'
#' McKinnon, K. I. M. (1998). Convergence of the Nelder-Mead simplex method to a
#' nonstationary point. \emph{SIAM Journal on Optimization} \strong{9}, 148--158.
#'
#' Gao, F. and Han, L. (2012). Implementing the Nelder-Mead simplex algorithm
#' with adaptive parameters. \emph{Computational Optimization and Applications}
#' \strong{51}, 259--277.
#'
#' @examples
#' nelder_mead()
#'
#' # a smooth problem, to show it works and to show what it costs
#' minimize(nelder_mead(), function(p) sum((p - c(1, 2))^2), c(0, 0))
#'
#' # what it is actually for: a sum of absolute deviations, whose minimiser is
#' # the median and whose derivative does not exist there
#' set.seed(1)
#' y <- rcauchy(101)
#' minimize(nelder_mead(), function(p) sum(abs(y - p)), par = 0)@par
#' median(y)
#'
#' @seealso \code{\link{compass}}, \code{\link{bundle}},
#'   \code{\link{crit_stationary}}
#' @export
nelder_mead <- function(criterion = crit_stationary(1e-8),
                        step = 0.1, adaptive = TRUE,
                        max_restarts = 3, degenerate_tol = 1e-6,
                        simplex = NULL,
                        maxit = 2000, max_eval = 20000,
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
#' @description The class \code{\link{compass}} instantiates.
#' @param step Initial poll size, relative to the starting value.
#' @param directions Either \code{"mads"} or \code{"coordinate"}.
#' @param opportunistic Whether to accept the first improvement found.
#' @param expand,shrink Factors applied to the poll size.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{compass}}
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
#' Uses no derivative, and unlike \code{\link{nelder_mead}} it comes with a
#' convergence theorem.
#'
#' @param criterion The stopping rule. Defaults to
#'   \code{crit_stationary(1e-8)}, on the poll size.
#' @param step Initial poll size, scaled by the largest coordinate of the
#'   starting value. Defaults to \code{0.1}.
#' @param directions \code{"mads"} (default) or \code{"coordinate"}; see
#'   Details.
#' @param opportunistic Move to the first improvement found rather than polling
#'   every direction? Defaults to \code{TRUE}.
#' @param expand Factor applied to the poll size after a success. Defaults to
#'   \code{2}.
#' @param shrink Factor applied after a failure. Defaults to \code{0.5}.
#' @param maxit Maximum iterations. Defaults to 2000.
#' @param max_eval Maximum objective evaluations. Defaults to 20000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 50.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' The directions form a \emph{positive spanning set}: any vector in the space
#' is a non-negative combination of them. That is the whole idea, and it is what
#' buys the theorem — if the current point is not stationary then some direction
#' in the set goes downhill, so a poll that fails everywhere is evidence about
#' the point rather than about the directions. A failed poll therefore licenses
#' shrinking the radius, and the limit points of a run in which the radius goes
#' to zero are stationary.
#'
#' \subsection{Poll directions}{
#' \describe{
#'   \item{\code{"coordinate"}}{the \eqn{2p} signed axes: this is compass
#'     search. Cheap, deterministic, reproducible. The theorem it enjoys assumes
#'     \eqn{f} is continuously differentiable.}
#'   \item{\code{"mads"}}{a fresh random orthonormal basis at every poll,
#'     taken plus and minus.}
#' }
#' The difference is not cosmetic on the problems this method exists for. When
#' \eqn{f} is merely Lipschitz — has kinks — a \emph{fixed} set of directions can
#' fail: a kink whose ridge runs diagonally is descended by no coordinate
#' direction, the poll fails at a point that is not stationary, and the run
#' stops there. The repair, which is the idea behind MADS, is to let the set of
#' directions used over the whole run become dense in the sphere, so no
#' direction of descent can be missed for ever; drawing a new orthonormal basis
#' at each poll achieves that with probability one. What is implemented here is
#' that idea rather than LTMADS as published, but it is the property the
#' convergence proof rests on.
#'
#' The cost is reproducibility: a random poll draws from \R's generator, so
#' \code{set.seed()} governs the run.
#' }
#'
#' \subsection{Opportunistic polling}{
#' Accepting the first improvement rather than the best costs a worse direction
#' and saves up to \eqn{2p - 1} evaluations. On an expensive objective that
#' trade is usually favourable, so it is the default; on a cheap one,
#' \code{opportunistic = FALSE} tends to need fewer iterations.
#' }
#'
#' @return An S7 object of class \code{Compass}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Torczon, V. (1997). On the convergence of pattern search algorithms.
#' \emph{SIAM Journal on Optimization} \strong{7}, 1--25.
#'
#' Audet, C. and Dennis, J. E. (2006). Mesh adaptive direct search algorithms
#' for constrained optimization. \emph{SIAM Journal on Optimization}
#' \strong{17}, 188--217.
#'
#' @examples
#' compass()
#'
#' # a kink running diagonally, which no coordinate direction descends
#' f <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
#' minimize(compass(), f, c(1, 0.5))@value
#' set.seed(1)
#' minimize(compass(directions = "mads"), f, c(1, 0.5))@value
#'
#' @seealso \code{\link{nelder_mead}}, \code{\link{bundle}},
#'   \code{\link{crit_stationary}}
#' @export
compass <- function(criterion = crit_stationary(1e-8),
                    step = 0.1,
                    directions = c("mads", "coordinate"),
                    opportunistic = TRUE, expand = 2, shrink = 0.5,
                    maxit = 2000, max_eval = 20000,
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
#' @description
#' The objective and a stationarity measure, but no gradient.
#' @details
#' A rule reading a gradient is refused rather than left testing \code{NULL} at
#' every iteration and never firing. \code{\link{crit_stationary}} is what
#' takes its place.
#' @param optimizer A \code{NelderMead} or \code{Compass} object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, NelderMead) <- function(optimizer)
  "stationarity"

#' @rdname optimizer_provides.NelderMead
#' @name optimizer_provides.Compass
#' @keywords internal
S7::method(optimizer_provides, Compass) <- function(optimizer)
  "stationarity"


#' @title Minimise by Nelder-Mead
#' @name minimize.NelderMead
#' @description Runs \code{\link{nelder_mead}} on the objective.
#' @param optimizer A \code{NelderMead} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}. \code{gr} and
#'   \code{he} are accepted and ignored: the method uses no derivative, and
#'   refusing them would force calling code to branch on the algorithm.
#' @return An \code{\link{optimizer_result}}.
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
      max_eval = as.integer(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0)
  }


#' @title Minimise by Pattern Search
#' @name minimize.Compass
#' @description Runs \code{\link{compass}} on the objective.
#' @param optimizer A \code{Compass} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}. \code{gr} and
#'   \code{he} are accepted and ignored.
#' @return An \code{\link{optimizer_result}}.
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
      max_eval = as.integer(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0, seed)
  }
