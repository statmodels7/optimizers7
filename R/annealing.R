#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for Simulated Annealing
#'
#' @description
#' An optimizer holding the annealing schedule and the step adaptation: which
#' proposal the walk draws from, the temperature and how it cools, how much
#' work is done at each level, and the acceptance rate the per-coordinate
#' steps are tuned to. Built by [sa()]. It is the package's only global
#' search and the only method whose result depends on the state of the random
#' number generator.
#'
#' @details
#' Beyond the seven properties every optimizer has, an `Sa` carries nine of
#' its own. They fall into three groups: the proposal (`visiting`, `step`),
#' the schedule (`t0`, `cooling`, `cycles`, `steps`), and the adaptation and
#' its stopping rule (`target_accept`, `adjust`, `n_eps`).
#'
#' @param visiting `"uniform"` or `"cauchy"`, which proposal the walk draws
#'   its moves from.
#' @param t0 The initial temperature, or `NULL` to calibrate it from the
#'   objective.
#' @param cooling The geometric cooling factor, in \eqn{(0, 1)}.
#' @param cycles,steps The work done at each temperature level:
#'   `cycles * steps * p` evaluations.
#' @param step The initial step, relative to the starting value.
#' @param target_accept The acceptance rate the step adaptation aims at.
#' @param adjust How hard the step is adjusted towards that rate.
#' @param n_eps How many temperature levels the stopping rule looks back over.
#'
#' @return An S7 object of class `Sa` inheriting from [optimizer()], with the
#'   nine properties above beside the seven shared ones.
#'
#' @seealso [sa()] for the constructor, [chain()] for handing its answer to a
#'   local method.
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
#' in. It is not a competitor to [bfgs()] or [newton()] on a
#' smooth problem, where it will be beaten by orders of magnitude, and the
#' intended use is to hand its result to one of them; see [chain()].
#'
#' # The adaptive step
#'
#' The parameters are moved one coordinate at a time, and after every `steps`
#' sweeps each coordinate's step is multiplied or divided according to how
#' often its moves were accepted, so that the rate is held inside a band
#' around `target_accept` (Corana et al. 1987). The adaptation makes the
#' method usable on a statistical objective, whose unconstrained coordinates
#' sit on scales orders of magnitude apart: one step length is wrong for all
#' of them. A coordinate accepting almost everything is being
#' proposed too timidly to explore, and one accepting almost nothing is being
#' thrown too far to land.
#'
#' # The proposal
#'
#' `"uniform"` draws the move uniformly on the coordinate's current step,
#' which together with the adaptation above is Corana's algorithm.
#' `"cauchy"` draws it from a Cauchy, whose heavy tail lets a run leave a
#' basin in one move instead of walking out of it. That is fast simulated
#' annealing (Szu and Hartley 1987), the \eqn{q = 2} member of the Tsallis
#' family.
#'
#' The difference is large on a surface with many basins. Over thirty seeds
#' on Rastrigin in two dimensions from \eqn{(4.4, -3.6)}, with 40 temperature
#' levels: the Cauchy proposal has a median value of `0.0051` and gets below
#' `0.01` on 19 of 30 runs, the uniform one a median of `0.94` and 1 of 30.
#' The default is nevertheless `"uniform"`, that being Corana's algorithm as
#' published and the better-behaved proposal where the objective is defined
#' on a bounded region.
#'
#' The general Tsallis visiting distribution at arbitrary \eqn{q} is not
#' offered. Its generator would have to be transcribed and could not be
#' checked against anything already here.
#'
#' # The temperature
#'
#' With `t0 = NULL` the initial temperature is calibrated from the
#' objective's own variation, by sampling proposals around the starting value
#' and setting \eqn{T_0} so that an average uphill move is accepted about
#' four times in five. That costs 20 evaluations, once, before the run
#' begins.
#'
#' A fixed number cannot serve, because the Metropolis probability
#' \eqn{\exp(-\Delta f / T)} compares the temperature against the objective's
#' own scale. At \eqn{T = 1} a step costing \eqn{10^{-6}} is accepted with
#' probability 1 and the run is a random walk; one costing \eqn{10^{6}} is
#' accepted with probability 0 and the method has become a local search.
#' Measured over eight seeds on \eqn{10^{-6}\sum(x - (1,2))^2}, the
#' calibrated run ends a median `0.0044` from the answer and `t0 = 1` a
#' median `0.375`.
#'
#' The other side of that is worth knowing too: on a **convex** objective too
#' cold a temperature costs nothing, there being one basin to find. On
#' \eqn{10^{6}\sum(x - (1,2))^2} the fixed `t0 = 1` reaches `4.6e-07` where
#' the calibration reaches `0.0013`. The calibration is insurance against a
#' multimodal surface, and on a problem that has no second basin it is a
#' small tax.
#'
#' # What the run returns, and what convergence means here
#'
#' The result is the best point **seen**, not the last one. An annealing run
#' wanders by construction, so its final iterate is a draw. The `value`
#' column of the trace is correspondingly the best so far, and is monotone
#' non-increasing.
#'
#' Whether the run converged is a separate question, never answered by the
#' schedule having finished. The stationarity measure reported is Corana's
#' own termination rule, by how much the best value has moved over the last
#' `n_eps` temperature levels, so [crit_stationary()] is that rule and no
#' second convention invented beside it. A run that merely exhausts `maxit`
#' reports `converged = FALSE`, which for a global search is the ordinary
#' outcome.
#'
#' The run is stochastic, so `set.seed()` governs it, and the state it began
#' from is recorded in the result's `seed` for repeating it.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_stationary(1e-8)`, read on Corana's measure above. A rule reading
#'   a gradient is refused when the run starts, this method computing none.
#' @param visiting `"uniform"` (default) or `"cauchy"`. Partial matching
#'   applies, and any other string is refused by name.
#' @param t0 The initial temperature, a single positive number, or `NULL`
#'   (the default) to calibrate it from the objective at a cost of 20
#'   evaluations.
#' @param cooling The factor the temperature is multiplied by at each level, a
#'   single number **strictly** inside \eqn{(0, 1)}. Defaults to `0.85`. `1`
#'   would never cool and `0` would freeze at once, so both endpoints are
#'   refused.
#' @param cycles How many step-adjustment cycles per temperature level, a
#'   positive whole number. Defaults to `3`.
#' @param steps How many sweeps of every coordinate per cycle, a positive
#'   whole number. Defaults to `10`. One temperature level therefore costs
#'   `cycles * steps * length(par)` evaluations, and a whole run costs that
#'   times `maxit`, plus one for the starting value and 20 more if the
#'   temperature is calibrated.
#' @param step The initial step, relative to the starting value. Defaults to
#'   `1`. It is only a starting point: the adaptation moves each coordinate's
#'   step from here within the first few cycles.
#' @param target_accept The acceptance rate the adaptation aims at, held
#'   inside a band of 0.1 either side. A single number strictly inside
#'   \eqn{(0.1, 0.9)}; anything outside is refused, an extreme target being
#'   unreachable by any step. Defaults to `0.5`.
#' @param adjust How hard the step is moved towards that rate, a single
#'   positive number. Defaults to `2`, which is Corana's own constant.
#' @param n_eps How many temperature levels the stationarity measure looks
#'   back over, a positive whole number. Defaults to `4`.
#' @param maxit Maximum temperature levels. Defaults to 100. This is the
#'   budget that decides how tightly the run finishes, the step adaptation
#'   shrinking the proposal as the acceptance rate falls with the
#'   temperature. Measured on \eqn{\sum(x - (1, 2))^2} from the origin at one
#'   seed, the distance to the answer goes `2.7e-02`, `1.8e-02`, `3.0e-03`,
#'   `2.3e-05` at 15, 30, 60 and 100 levels, for 921, 1821, 3621 and 6021
#'   evaluations.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many levels. Defaults to 10.
#' @param keep_trace Store the path? Defaults to `FALSE`. The `value` column
#'   is the best so far.
#'
#' @return An S7 object of class [Sa], inheriting from [optimizer()], to be
#'   handed to [minimize()].
#'
#' @references
#' Corana, A., Marchesi, M., Martini, C. and Ridella, S. (1987). Minimizing
#' multimodal functions of continuous variables with the simulated annealing
#' algorithm. *ACM Transactions on Mathematical Software* **13**,
#' 262--280.
#'
#' Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by
#' simulated annealing. *Science* **220**, 671--680.
#'
#' Szu, H. and Hartley, R. (1987). Fast simulated annealing. *Physics
#' Letters A* **122**, 157--162.
#'
#' @examples
#' sa()
#'
#' # Rastrigin has 121 local minima. A local method stops in the basin it
#' # started in; this one leaves it.
#' rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
#' set.seed(1)
#' minimize(sa(), rastrigin, c(4.4, -3.6))@value    # 0.914
#' minimize(bfgs(), rastrigin, c(4.4, -3.6))@value  # 7.96
#'
#' # The heavy-tailed proposal leaves a basin in one move, and on this surface
#' # that is worth two orders of magnitude.
#' set.seed(2)
#' minimize(sa(visiting = "cauchy", maxit = 40), rastrigin, c(4.4, -3.6))@value
#' set.seed(2)
#' minimize(sa(visiting = "uniform", maxit = 40), rastrigin, c(4.4, -3.6))@value
#'
#' # The reported value is the best seen, and the trace records the same:
#' # it never rises, though the walk itself does.
#' set.seed(5)
#' r <- minimize(sa(maxit = 30, keep_trace = TRUE), rastrigin, c(4.4, -3.6))
#' all(diff(r@trace$value) <= 0)
#' all.equal(r@value, min(r@trace$value))
#'
#' # The run repeats exactly from the state it recorded.
#' set.seed(11)
#' a <- minimize(sa(maxit = 20), rastrigin, c(4.4, -3.6))
#' assign(".Random.seed", a@seed, globalenv())
#' b <- minimize(sa(maxit = 20), rastrigin, c(4.4, -3.6))
#' identical(a@par, b@par)
#'
#' # A gradient rule is refused when the run starts, this method computing none.
#' try(minimize(sa(criterion = crit_grad()), rastrigin, c(4.4, -3.6)))
#'
#' @seealso [chain()] for handing the answer to a local method,
#'   [multistart()] for the other way of covering several basins,
#'   [crit_stationary()] for the rule this method's measure feeds.
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
#'
#' @description
#' Reports `"stationarity"` and nothing else. The method draws proposals and
#' compares values; it computes no derivative, so a rule reading a gradient
#' is refused when the run starts instead of sitting there testing `NULL` and
#' never firing.
#'
#' @details
#' The measure offered is Corana's termination rule: by how much the best
#' value has moved over the last `n_eps` temperature levels. It is read by
#' [crit_stationary()], which is this method's default rule.
#'
#' @param optimizer An `Sa` object.
#'
#' @return The character vector `"stationarity"`.
#'
#' @examples
#' optimizer_provides(sa())
#'
#' # So a gradient rule is refused, with both names in the message.
#' rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
#' try(minimize(sa(criterion = crit_grad()), rastrigin, c(1, 1)))
#'
#' @keywords internal
S7::method(optimizer_provides, Sa) <- function(optimizer) "stationarity"


#' @title Minimize by Simulated Annealing
#' @name minimize.Sa
#'
#' @description
#' Runs [sa()] on the objective: a Metropolis walk over one coordinate at a
#' time, with the step of each coordinate adapted to hold its acceptance rate
#' near the target, and the temperature cooled geometrically. Records the
#' state of the random number generator before the first draw, so the run can
#' be repeated.
#'
#' @param optimizer An `Sa` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` and `he` are
#'   accepted and ignored, the method using no derivative; refusing them
#'   would force calling code to branch on the algorithm. Bounds are taken
#'   and removed by reparametrization like any other method's.
#'
#' @return An [optimizer_result()] holding the best point seen, with `seed`
#'   filled in and `converged` `TRUE` only when Corana's rule fired.
#'
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
#' The annealing loop written in R, the twin of the compiled `sa_run()`. It
#' exists so that the kernel has something to be compared against that shares
#' none of its code, and it is used by the tests alone.
#'
#' @details
#' It draws from \R's generator in the same order as the kernel, so from one
#' seed the two are the same run and the comparison needs no tolerance at
#' all. The order is the thing that has to match: one uniform per proposal,
#' and a second one for the Metropolis test only when the proposal is uphill.
#' A transcription that consumes a uniform unconditionally passes every test
#' written on the answer and fails this one immediately.
#'
#' @param fn The objective, a function of the parameter vector.
#' @param par The starting value, a numeric vector.
#' @param cauchy `TRUE` for the Cauchy proposal, `FALSE` for the uniform one.
#' @param t0 The initial temperature, or a non-positive value to calibrate
#'   it. `-1` is what [minimize.Sa()] passes for `t0 = NULL`.
#' @param cooling,cycles,steps,step,target_accept,adjust,n_eps,maxit As in
#'   [sa()].
#'
#' @return A list with `par` (the best point seen), `value` (the objective
#'   there) and `n_value` (the evaluation count). It reports no trace and no
#'   stationarity measure, the comparison being of the walk alone.
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
