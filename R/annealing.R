#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for Simulated Annealing
#' @description The class \code{\link{sa}} instantiates.
#' @param visiting Which proposal the walk uses.
#' @param t0 The initial temperature, or \code{NULL} to calibrate it.
#' @param cooling The geometric cooling factor.
#' @param cycles,steps The work done at each temperature.
#' @param step The initial step, relative to the starting value.
#' @param target_accept The acceptance rate the step adaptation aims at.
#' @param adjust How hard the step is adjusted towards that rate.
#' @param n_eps How many temperature levels the stopping rule looks back over.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{sa}}
#' @name Sa-class
#' @aliases Sa
#' @keywords internal
Sa <- S7::new_class("Sa", parent = optimizer,
  properties = list(
    visiting       = S7::class_character,
    t0             = S7::class_any,
    cooling        = S7::class_numeric,
    cycles         = S7::class_numeric,
    steps          = S7::class_numeric,
    step           = S7::class_numeric,
    target_accept  = S7::class_numeric,
    adjust         = S7::class_numeric,
    n_eps          = S7::class_numeric
  ))


#' @title Simulated Annealing With an Adaptive Step
#'
#' @description
#' A global search that accepts uphill moves with a probability falling as the
#' run cools, with one step length per coordinate adjusted in flight to hold
#' that coordinate's acceptance rate near a target.
#'
#' @details
#' The method exists here for the problems the local ones cannot start on: a
#' multimodal objective where the answer depends on which basin the run began
#' in. It is not a competitor to \code{\link{bfgs}} or \code{\link{newton}} on a
#' smooth problem, where it will be beaten by orders of magnitude, and the
#' intended use is to hand its result to one of them --- see
#' \code{\link{chain}}.
#'
#' \strong{The adaptive step.} The parameters are moved one coordinate at a
#' time, and after every \code{steps} sweeps each coordinate's step is
#' multiplied or divided according to how often its moves were accepted, so
#' that the rate is held inside a band around \code{target_accept}
#' (Corana et al. 1987). This is what makes the method usable on a statistical
#' objective: the coordinates of an unconstrained parameter vector sit on
#' scales orders of magnitude apart, and one step length is wrong for all of
#' them. A coordinate accepting almost everything is being proposed too
#' timidly to explore; one accepting almost nothing is being thrown too far to
#' land.
#'
#' \strong{The proposal.} \code{"uniform"} draws the move uniformly on the
#' coordinate's current step, which with the adaptation above is Corana's
#' algorithm. \code{"cauchy"} draws it from a Cauchy, whose heavy tail lets a
#' run leave a basin in one move rather than walking out of it; that is fast
#' simulated annealing (Szu and Hartley 1987), and it is the \eqn{q = 2} member
#' of the Tsallis family. The general Tsallis visiting distribution at
#' arbitrary \eqn{q} is deliberately not offered: its generator would have to
#' be transcribed and could not be checked against anything already here, and
#' an unverified generator is worse than a case that can be verified.
#'
#' \strong{The temperature.} With \code{t0 = NULL} the initial temperature is
#' calibrated from the objective's own variation, by sampling proposals around
#' the starting value and setting \eqn{T_0} so that an average uphill move is
#' accepted about four times in five. A fixed number cannot serve: on an
#' objective of size \eqn{10^6} every move is accepted and the run is a random
#' walk, on one of size \eqn{10^{-6}} none is and it is a poor local search.
#'
#' \strong{What the run returns is the best point SEEN}, not the last one. An
#' annealing run wanders by construction, so its final iterate is a draw and
#' not an answer.
#'
#' \strong{Whether it converged is a separate question} and is never answered
#' by the schedule having finished. The stationarity measure reported is
#' Corana's own termination rule --- by how much the best value has moved over
#' the last \code{n_eps} temperature levels --- so \code{\link{crit_stationary}}
#' IS that rule rather than a second convention invented beside it. A run that
#' merely exhausts \code{maxit} reports \code{converged = FALSE}, which for a
#' global search is the ordinary outcome and not a failure.
#'
#' \strong{The run is stochastic}, so \code{set.seed()} governs it and the
#' state it started from is recorded in the result.
#'
#' @param criterion The stopping rule. Defaults to
#'   \code{crit_stationary(1e-8)}, read on the measure described above.
#' @param visiting \code{"uniform"} (default) or \code{"cauchy"}.
#' @param t0 The initial temperature. \code{NULL}, the default, calibrates it
#'   from the objective.
#' @param cooling The factor the temperature is multiplied by at each level.
#'   Defaults to \code{0.85}.
#' @param cycles How many step-adjustment cycles per temperature level.
#'   Defaults to \code{3}.
#' @param steps How many sweeps of every coordinate per cycle. Defaults to
#'   \code{10}. One temperature level therefore costs
#'   \code{cycles * steps * length(par)} evaluations.
#' @param step The initial step, relative to the starting value. Defaults to
#'   \code{1}.
#' @param target_accept The acceptance rate the adaptation aims at, held inside
#'   a band of 0.1 either side. Defaults to \code{0.5}.
#' @param adjust How hard the step is moved towards that rate. Defaults to
#'   \code{2}.
#' @param n_eps How many temperature levels the stationarity measure looks back
#'   over. Defaults to \code{4}.
#' @param maxit Maximum temperature levels. Defaults to 100. It is the budget
#'   that decides how tightly the run finishes: the step adaptation shrinks the
#'   proposal as the acceptance rate falls with the temperature, so on a
#'   quadratic the answer improves from \eqn{2\times10^{-1}} at 15 levels to
#'   \eqn{5\times10^{-5}} at 100.
#' @param max_eval Maximum objective evaluations. Defaults to \code{Inf}.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many levels. Defaults to 10.
#' @param keep_trace Store the path? Defaults to \code{FALSE}.
#'
#' @return An S7 object of class \code{Sa}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Corana, A., Marchesi, M., Martini, C. and Ridella, S. (1987). Minimizing
#' multimodal functions of continuous variables with the simulated annealing
#' algorithm. \emph{ACM Transactions on Mathematical Software} \strong{13},
#' 262--280.
#'
#' Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by
#' simulated annealing. \emph{Science} \strong{220}, 671--680.
#'
#' Szu, H. and Hartley, R. (1987). Fast simulated annealing. \emph{Physics
#' Letters A} \strong{122}, 157--162.
#'
#' @examples
#' sa()
#'
#' # an objective with many local minima, where a local method stops in the
#' # basin it started in and this one does not
#' rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
#' set.seed(1)
#' minimize(sa(), rastrigin, c(4.4, -3.6))@value
#' minimize(bfgs(), rastrigin, c(4.4, -3.6))@value
#'
#' @seealso \code{\link{chain}}, \code{\link{multistart}},
#'   \code{\link{crit_stationary}}
#' @export
sa <- function(criterion = crit_stationary(),
               visiting = c("uniform", "cauchy"),
               t0 = NULL, cooling = 0.85, cycles = 3, steps = 10,
               step = 1, target_accept = 0.5, adjust = 2, n_eps = 4,
               maxit = 100, max_eval = Inf,
               verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  visiting <- match.arg(visiting)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  if (!is.null(t0) && (length(t0) != 1L || !is.numeric(t0) || is.na(t0) ||
                       t0 <= 0)) {
    stop("'t0' must be a single positive number, or NULL to calibrate it.",
         call. = FALSE)
  }
  if (length(cooling) != 1L || !is.numeric(cooling) || is.na(cooling) ||
      cooling <= 0 || cooling >= 1) {
    stop("'cooling' must be a single number in (0, 1).", call. = FALSE)
  }
  whole <- function(v, nm) {
    if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 1 ||
        v != round(v)) {
      stop("'", nm, "' must be a single positive whole number.", call. = FALSE)
    }
  }
  whole(cycles, "cycles")
  whole(steps, "steps")
  whole(n_eps, "n_eps")
  if (length(target_accept) != 1L || !is.numeric(target_accept) ||
      is.na(target_accept) || target_accept <= 0.1 || target_accept >= 0.9) {
    stop("'target_accept' must be a single number in (0.1, 0.9).",
         call. = FALSE)
  }
  if (length(adjust) != 1L || !is.numeric(adjust) || is.na(adjust) ||
      adjust <= 0) {
    stop("'adjust' must be a single positive number.", call. = FALSE)
  }
  Sa(name = paste0("simulated annealing (", visiting, ")"),
     criterion = criterion, maxit = maxit, max_eval = max_eval,
     verbose = verbose, refresh = refresh, keep_trace = keep_trace,
     visiting = visiting, t0 = t0, cooling = cooling, cycles = cycles,
     steps = steps, step = step, target_accept = target_accept,
     adjust = adjust, n_eps = n_eps)
}


#' @title What Simulated Annealing Can Offer a Stopping Rule
#' @name optimizer_provides.Sa
#' @description
#' The objective and a stationarity measure, but no gradient.
#' @details
#' The measure is Corana's termination rule, by how much the best value has
#' moved over the last \code{n_eps} temperature levels. A rule reading a
#' gradient is rejected at construction rather than left testing \code{NULL}
#' and never firing.
#' @param optimizer An \code{Sa} object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, Sa) <- function(optimizer) "stationarity"


#' @title Minimize by Simulated Annealing
#' @name minimize.Sa
#' @description Runs \code{\link{sa}} on the objective.
#' @param optimizer An \code{Sa} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}. \code{gr}
#'   and \code{he} are accepted and ignored: the method uses no derivative, and
#'   rejecting them would force calling code to branch on the algorithm.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Sa) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
    bounds <- check_bounds(lower, upper, par)

    # Every proposal is drawn from R's generator, so the run is reproducible
    # only if the state it started from is known. Recorded for the same reason
    # a mads poll records it.
    seed <- capture_seed()

    t0 <- proc.time()[["elapsed"]]
    out <- sa_run(
      spec = spec, par = as.numeric(par),
      criterion = optimizer@criterion, crit_fn = crit_met,
      cauchy = identical(optimizer@visiting, "cauchy"),
      t0 = if (is.null(optimizer@t0)) -1 else as.numeric(optimizer@t0),
      cooling = optimizer@cooling,
      cycles = as.integer(optimizer@cycles),
      steps = as.integer(optimizer@steps),
      step = optimizer@step,
      target_accept = optimizer@target_accept,
      adjust = optimizer@adjust,
      n_eps = as.integer(optimizer@n_eps),
      maxit = as.integer(optimizer@maxit),
      max_eval = budget_int(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0, seed)
  }


#' The Annealing Loop in R
#'
#' @description
#' The loop \code{sa_run()} replaces, kept so the compiled route has something
#' to be compared against that shares none of its code.
#'
#' @details
#' It draws from \R's generator in the same order as the kernel, so from one
#' seed the two runs are the same run and the comparison needs no tolerance.
#' The order is what has to match: one uniform per proposal, and a second one
#' for the Metropolis test ONLY when the proposal is uphill, which is where a
#' transcription most easily drifts.
#'
#' @param fn The objective.
#' @param par The starting value.
#' @param cauchy Whether the proposal is Cauchy rather than uniform.
#' @param t0 The initial temperature, or a non-positive value to calibrate it.
#' @param cooling,cycles,steps,step,target_accept,adjust,n_eps,maxit As in
#'   \code{\link{sa}}.
#'
#' @return A list with \code{par}, \code{value} and \code{n_value}.
#'
#' @keywords internal
sa_run_r <- function(fn, par, cauchy = FALSE, t0 = -1, cooling = 0.85,
                     cycles = 3, steps = 10, step = 1, target_accept = 0.5,
                     adjust = 2, n_eps = 4, maxit = 100) {
  n_value <- 0L
  safe <- function(x) {
    n_value <<- n_value + 1L
    v <- fn(x)
    if (is.finite(v)) v else Inf
  }
  x <- as.numeric(par)
  p <- length(x)
  f <- safe(x)
  best_x <- x
  best_f <- f
  v <- step * pmax(1, abs(x))

  if (!(t0 > 0)) {
    total <- 0
    m <- 0L
    n_calib <- max(10L, 10L * p)
    for (k in seq_len(n_calib)) {
      j <- min(p, floor(stats::runif(1) * p) + 1)
      trial <- x
      trial[j] <- trial[j] + v[j] * (2 * stats::runif(1) - 1)
      ft <- safe(trial)
      if (is.finite(ft)) { total <- total + abs(ft - f); m <- m + 1L }
    }
    mean_jump <- if (m > 0L) total / m else 0
    t0 <- if (mean_jump > 0) mean_jump / (-log(0.8)) else 1
  }
  temp <- t0

  recent <- numeric(0)
  stat <- Inf
  for (it in seq_len(maxit)) {
    for (cyc in seq_len(cycles)) {
      n_acc <- numeric(p)
      for (m in seq_len(steps)) {
        for (j in seq_len(p)) {
          dev <- if (cauchy) stats::rcauchy(1) else 2 * stats::runif(1) - 1
          trial <- x
          trial[j] <- trial[j] + v[j] * dev
          ft <- safe(trial)
          take <- if (ft <= f) TRUE else if (!is.finite(ft) || temp <= 0) FALSE
            else stats::runif(1) < exp(-(ft - f) / temp)
          if (take) {
            x <- trial
            f <- ft
            n_acc[j] <- n_acc[j] + 1
            if (f < best_f) { best_f <- f; best_x <- x }
          }
        }
      }
      lo <- target_accept - 0.1
      hi <- target_accept + 0.1
      for (j in seq_len(p)) {
        ratio <- n_acc[j] / steps
        if (ratio > hi) {
          v[j] <- v[j] * (1 + adjust * (ratio - hi) / (1 - hi))
        } else if (ratio < lo) {
          v[j] <- v[j] / (1 + adjust * (lo - ratio) / lo)
        }
        if (!is.finite(v[j]) || v[j] <= 0) v[j] <- step
      }
    }
    recent <- c(recent, f)
    if (length(recent) > n_eps) {
      recent <- recent[-1L]
      stat <- max(abs(f - best_f), max(abs(f - recent)))
    }
    temp <- temp * cooling
  }
  list(par = best_x, value = best_f, n_value = n_value)
}
