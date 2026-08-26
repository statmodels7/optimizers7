#' @include optimizer_class.R
#' @include generics.R
#' @include line_search.R
NULL

# Newton, BFGS and L-BFGS. Each is a constructor and a method that names its
# direction; the loop, the line search, the stopping rule and the reporting are
# all shared, which is what the frame was built for.


#' @title S7 Class for Newton's Method
#'
#' @description
#' An optimizer holding the initial step length, the line search, and the two
#' settings that decide how an indefinite Hessian is repaired. Built by
#' [newton()]. It is the only shipped method that reads the `he` argument of
#' [minimize()].
#'
#' @details
#' Beyond the seven properties every optimizer has, a `Newton` carries four:
#' `step` and `line_search`, shared with the other line-search methods, and
#' `hessian_mod` and `floor`, which are its own.
#'
#' @param hessian_mod How an indefinite Hessian is repaired, `"eigen"` or
#'   `"ridge"`.
#' @param floor The smallest eigenvalue the repaired Hessian may have.
#' @param step,line_search The initial step length and the [line_search()]
#'   object, as in [newton()].
#'
#' @return An S7 object of class `Newton` inheriting from [optimizer()], with
#'   the four properties above beside the seven shared ones.
#'
#' @seealso [newton()] for the constructor, [bfgs()] for the method that
#'   needs no Hessian.
#' @name Newton-class
#' @aliases Newton
#' @keywords internal
Newton <- S7::new_class("Newton", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    hessian_mod = S7::class_character,
                    floor = S7::class_numeric))


#' @title Newton's Method with a Modified Hessian
#'
#' @description
#' Solves \eqn{H d = -g} for the direction, repairing \eqn{H} when it is not
#' positive definite, and then line searches along it. Converges
#' quadratically near a minimum: on Rosenbrock from the customary start it
#' reaches a value of `3.7e-21` in 21 iterations, 29 objective evaluations
#' and 21 Hessians.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`; see
#'   [crit_any()] for what that disjunction buys and costs.
#' @param hessian_mod How to repair an indefinite Hessian: `"eigen"`
#'   (default) or `"ridge"`. Partial matching applies and any other string is
#'   refused, naming both.
#' @param floor The smallest eigenvalue the repaired Hessian is allowed, a
#'   single positive number. Defaults to `1e-8`, relative to the largest
#'   eigenvalue.
#' @param step Initial step length offered to the line search, a single
#'   positive number. Defaults to `1`, which is the natural Newton step and
#'   is accepted unchanged near the solution.
#' @param line_search A [line_search()] object. Defaults to [armijo()], the
#'   curvature condition [bfgs()] needs being unnecessary here because the
#'   curvature is read from the Hessian.
#' @param maxit,max_eval,verbose,refresh,keep_trace As in [optimizer()].
#'   `maxit` defaults to 200 here where [optimizer()] uses 500, a Newton run
#'   that has not arrived in 200 iterations being in trouble of another kind.
#'
#' @details
#' # Why the Hessian has to be repaired
#'
#' \eqn{H^{-1}g} is a descent direction only when \eqn{H} is positive
#' definite. Where the objective curves downwards the unmodified step points
#' towards a saddle or a maximum, and no line search can rescue it: every
#' step along an ascent direction increases the objective. That is the
#' ordinary situation far from the solution, not an edge case.
#'
#' Both repairs begin with an attempted Cholesky factorization, which when it
#' succeeds is at once the test for positive definiteness and the solve.
#' When it fails:
#'
#' \describe{
#'   \item{`"eigen"`}{decompose \eqn{H} and raise every eigenvalue below
#'     `floor` to it. The direction is then the Newton one in the subspace
#'     where the curvature is trustworthy, and gradient-like in the rest. Costs a
#'     symmetric eigendecomposition and gives the best-conditioned repair.}
#'   \item{`"ridge"`}{add \eqn{\tau I} with \eqn{\tau} doubling until the
#'     factorization succeeds. This is Levenberg's idea: it interpolates between
#'     the Newton step at \eqn{\tau = 0} and a scaled steepest-descent step for
#'     large \eqn{\tau}. Cheaper, and blunter.}
#' }
#'
#' Which repair fired is recorded in the trace, under the names
#' `hessian modified` and `hessian modified (capped)`, so a run that spent
#' its time repairing can be told from one that spent it converging.
#'
#' # A repaired step is capped and an unrepaired one is not
#'
#' Where the factorization succeeds, \eqn{H^{-1}g} is the Newton step and
#' `step = 1` is its own unit, so it is taken as it stands. Where it fails,
#' the component of the direction in the floored subspace is
#' \eqn{g_i/\lambda_{\mathrm{floor}}}, whose size is set by `floor` and by no
#' curvature of the objective. The direction is therefore scaled to
#' \eqn{\min(1, 1/\lVert d\rVert_\infty)}, a displacement of order one in the
#' parameters. That is the rule [bb()] applies when its secant pair reports
#' no curvature and the one [bfgs()] and [lbfgs()] apply to their first
#' direction, and it can only ever shorten a step, so a run whose repaired
#' steps were already of order one is untouched. The same scaling reaches the
#' gradient the method falls back on when a solve fails.
#'
#' Without the cap the length of a repaired step is unbounded and the line
#' search is the only thing bounding it, at one objective evaluation per
#' backtrack. Measured on a marginal criterion whose outer Hessian is
#' indefinite at ordinary points, a gradient of 29 along a direction floored
#' at \eqn{10^{-8}\lambda_{\max}} gave a step of \eqn{4\times 10^{4}} on a
#' log scale, and each of the 23 backtracks that followed was a whole
#' penalized refit that could not be evaluated.
#'
#' # On the Hessian itself
#'
#' When `he` is not supplied to [minimize()] the Hessian is obtained by
#' differencing the gradient. With an analytic gradient that is one numerical
#' differentiation and is acceptable. When the gradient is *also* differenced
#' it is two composed, the one place in the package where that happens, and
#' the cost shows: one Newton iteration on a quadratic in 20 unknowns takes
#' 1682 objective evaluations without a gradient and 42 gradient evaluations
#' with one. An objective with neither derivative is better served by
#' [bfgs()], which needs no Hessian at all.
#'
#' @return An S7 object of class [Newton], inheriting from [optimizer()], to
#'   be handed to [minimize()].
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                           200 * (p[2] - p[1]^2))
#' rosen_he <- function(p) matrix(
#'   c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
#'     -400 * p[1], 200), 2, 2)
#'
#' minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
#'
#' # A saddle, where the unrepaired step points the wrong way. The Hessian at
#' # the start has eigenvalues 2 and -1.99, and both repairs reach the minimum
#' # in five iterations, reporting the repair in the trace.
#' sad <- function(p) p[1]^2 - p[2]^2 + 0.1 * p[2]^4
#' sg  <- function(p) c(2 * p[1], -2 * p[2] + 0.4 * p[2]^3)
#' sh  <- function(p) matrix(c(2, 0, 0, -2 + 1.2 * p[2]^2), 2, 2)
#' eigen(sh(c(0.5, 0.1)), only.values = TRUE)$values
#'
#' r <- minimize(newton(keep_trace = TRUE), sad, c(0.5, 0.1), gr = sg, he = sh)
#' c(value = r@value, iterations = r@iterations)
#' unique(r@trace$safeguard)
#'
#' # summary() counts the repairs, which is how a struggling run is read.
#' summary(minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1),
#'                  gr = rosen_gr, he = rosen_he))
#'
#' @seealso [bfgs()] and [lbfgs()] when no Hessian is available,
#'   [armijo()] for the line search, [summary.optimizer_result()] for the
#'   safeguard counts.
#' @references
#' Gill, P. E., Murray, W. and Wright, M. H. (1981).
#' *Practical Optimization*. Academic Press, London.
#'
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
#' 2nd edition. Springer, New York.
#'
#' @export
newton <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                   hessian_mod = c("eigen", "ridge"), floor = 1e-8,
                   step = 1, line_search = armijo(),
                   maxit = 200, max_eval = Inf,
                   verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  hessian_mod <- match.arg(hessian_mod)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  if (length(floor) != 1L || !is.numeric(floor) || is.na(floor) || floor <= 0) {
    stop("'floor' must be a single positive number.", call. = FALSE)
  }
  Newton(name = "Newton", criterion = criterion, maxit = maxit,
         max_eval = max_eval, verbose = verbose, refresh = refresh,
         keep_trace = keep_trace, step = step, line_search = line_search,
         hessian_mod = hessian_mod, floor = floor)
}


#' @title S7 Class for BFGS
#'
#' @description
#' An optimizer holding the initial step length, the line search, and the two
#' settings that govern what happens when a secant pair carries no usable
#' curvature. Built by [bfgs()]. The inverse-Hessian approximation itself is
#' not a property: it is built inside the run and discarded with it.
#'
#' @details
#' Beyond the seven properties every optimizer has, a `Bfgs` carries four:
#' `step` and `line_search`, shared with the other line-search methods, and
#' `curv_tol` and `max_skip`, which are its own.
#'
#' @param curv_tol The curvature threshold below which the update is skipped.
#' @param max_skip Consecutive skips before the approximation is reset to the
#'   identity.
#' @param step,line_search The initial step length and the [line_search()]
#'   object, as in [bfgs()].
#'
#' @return An S7 object of class `Bfgs` inheriting from [optimizer()], with
#'   the four properties above beside the seven shared ones.
#'
#' @seealso [bfgs()] for the constructor, [Lbfgs] for the limited-memory
#'   version.
#' @name Bfgs-class
#' @aliases Bfgs
#' @keywords internal
Bfgs <- S7::new_class("Bfgs", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    curv_tol = S7::class_numeric,
                    max_skip = S7::class_numeric))


#' @title BFGS
#'
#' @description
#' Builds an approximation to the inverse Hessian from successive gradients,
#' so the direction is a matrix-vector product and no second derivatives are
#' ever required. The default method for a smooth problem of moderate size:
#' on Rosenbrock from the customary start it converges in 35 iterations and
#' 49 objective evaluations with the gradient supplied.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`.
#' @param curv_tol The update is skipped when
#'   \eqn{s^\top y \le \texttt{curv\_tol}\,\lVert s\rVert\lVert y\rVert}.
#'   Defaults to `1e-10`. A larger value skips more often. Zero is admitted
#'   and gives the textbook condition \eqn{s^\top y > 0}; a negative value is
#'   refused, because the comparison would then hold for every pair and the
#'   skip protection would be off without saying so.
#' @param max_skip Consecutive skipped updates before the approximation is
#'   reset to the identity. Defaults to 5, and has to be a whole number.
#'   Zero is admitted and resets on the first skip; a negative value is
#'   refused, the comparison being `>=`, which would make it behave as zero.
#' @param step,line_search,maxit,max_eval,verbose,refresh,keep_trace As in
#'   [newton()], except that `line_search` defaults to [wolfe()] and `maxit`
#'   to 500. See below for why the line search matters here.
#'
#' @details
#' # Notation
#'
#' \eqn{s = x_{new} - x_{old}} is the **secant vector** and
#' \eqn{y = g_{new} - g_{old}} the change in the gradient. The scalar step
#' length a line search chooses is \eqn{\alpha}, a different quantity; the
#' two are related by \eqn{s = \alpha d} for the direction \eqn{d}.
#'
#' # Why the line search is the strong Wolfe one
#'
#' The update is meaningful only when \eqn{s^\top y > 0}, and the Wolfe
#' curvature condition is exactly the guarantee that it holds: Armijo
#' backtracking can accept a step so short that the gradient has barely
#' moved, leaving a pair with no curvature in it.
#'
#' In practice the guarantee is worth less than it sounds, and the honest
#' measurement is worth having. Over the eight problems of
#' [test_problems()], BFGS under Armijo skips an update on three of them,
#' once each, and matches Wolfe elsewhere: Rosenbrock 39 iterations against
#' 35, himmelblau 11 against 10, and on the non-smooth `abs_sum` Armijo is
#' the **better** of the two, reaching `6.8e-08` in 37 iterations where Wolfe
#' stops at `1.3e-02` reporting failure. Wolfe is the default because it
#' makes the update sound by construction; Armijo is a reasonable choice
#' where evaluations are dear.
#'
#' # Skipping and resetting
#'
#' When the curvature condition fails anyway the update is **skipped**. A
#' small \eqn{s^\top y} makes \eqn{\rho = 1/s^\top y} enormous and one bad
#' step destroys the accumulated approximation; a stale but sound matrix
#' beats a fresh but corrupted one. After `max_skip` consecutive skips there
#' is no curvature information left worth keeping and the matrix is reset to
#' the identity, so the method restarts as steepest descent and rebuilds.
#' Both events appear in the trace, as `bfgs update skipped` and
#' `bfgs reset`.
#'
#' The first accepted pair rescales the identity by \eqn{s^\top y / y^\top y}.
#' Without it the first quasi-Newton step is taken with a unit Hessian, which
#' on a badly scaled problem has the wrong magnitude entirely and wastes a
#' line search discovering so. The first direction is also scaled to
#' \eqn{\min(1, 1/\lVert g\rVert_\infty)}, a displacement of order one in the
#' parameters, which can only ever shorten it.
#'
#' @return An S7 object of class [Bfgs], inheriting from [optimizer()], to be
#'   handed to [minimize()].
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rg <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                     200 * (p[2] - p[1]^2))
#' minimize(bfgs(), rosen, c(-1.2, 1), gr = rg)
#'
#' # Armijo instead of Wolfe: the update is skipped once here, which the trace
#' # reports, and the run still arrives.
#' a <- minimize(bfgs(line_search = armijo(), keep_trace = TRUE), rosen,
#'               c(-1.2, 1), gr = rg)
#' c(iterations = a@iterations, skips = sum(a@trace$safeguard == "bfgs update skipped"))
#'
#' # Forcing the skip: a curvature threshold nothing can meet exhausts the
#' # budget, and the trace names both safeguards.
#' b <- minimize(bfgs(curv_tol = 1e10, keep_trace = TRUE), rosen, c(-1.2, 1),
#'               gr = rg)
#' c(converged = b@converged, value = b@value)
#' unique(b@trace$safeguard)
#'
#' @seealso [lbfgs()] for many parameters, [newton()] when a Hessian is
#'   available, [wolfe()] for the line search this method wants.
#' @references
#' Broyden, C. G. (1970). The convergence of a class of double-rank
#' minimization algorithms. *IMA Journal of Applied Mathematics*
#' **6**, 76--90. The update was obtained independently the same
#' year by Fletcher, Goldfarb and Shanno, whence the name.
#'
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
#' 2nd edition. Springer, New York.
#'
#' @export
bfgs <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                 curv_tol = 1e-10, max_skip = 5,
                 step = 1, line_search = wolfe(),
                 maxit = 500, max_eval = Inf,
                 verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_nonneg(curv_tol, "curv_tol")
  check_whole(max_skip, "max_skip")
  Bfgs(name = "BFGS", criterion = criterion, maxit = maxit,
       max_eval = max_eval, verbose = verbose, refresh = refresh,
       keep_trace = keep_trace, step = step, line_search = line_search,
       curv_tol = curv_tol, max_skip = max_skip)
}


#' @title S7 Class for Limited-Memory BFGS
#'
#' @description
#' An optimizer holding how many secant pairs to keep, the curvature
#' threshold a pair must meet to be stored, and the usual step length and
#' line search. Built by [lbfgs()]. The pairs themselves live inside the run.
#'
#' @details
#' Beyond the seven properties every optimizer has, an `Lbfgs` carries four:
#' `step` and `line_search`, and `memory` and `curv_tol`, which are its own.
#' Where [Bfgs] has `max_skip`, this class has nothing corresponding: a pair
#' failing the curvature test is discarded rather than skipped, so there is
#' no accumulated matrix to protect.
#'
#' @param memory How many secant pairs to keep.
#' @param curv_tol The curvature threshold below which a pair is not stored.
#' @param step,line_search The initial step length and the [line_search()]
#'   object, as in [lbfgs()].
#'
#' @return An S7 object of class `Lbfgs` inheriting from [optimizer()], with
#'   the four properties above beside the seven shared ones.
#'
#' @seealso [lbfgs()] for the constructor, [Bfgs] for the full-matrix
#'   version.
#' @name Lbfgs-class
#' @aliases Lbfgs
#' @keywords internal
Lbfgs <- S7::new_class("Lbfgs", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    memory = S7::class_numeric,
                    curv_tol = S7::class_numeric))


#' @title Limited-Memory BFGS
#'
#' @description
#' BFGS without ever forming the matrix: only the last `memory` secant pairs
#' are kept, and the direction comes from the two-loop recursion. Costs
#' \eqn{O(mp)} in time and memory where [bfgs()] costs \eqn{O(p^2)}, and
#' returns the same direction the full matrix would when the pairs are the
#' same.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`.
#' @param memory How many secant pairs to keep, a single positive whole
#'   number. Defaults to 10. A fractional value is refused rather than
#'   truncated, the kernel reading the count through `as.integer()`.
#' @param curv_tol A pair is discarded when
#'   \eqn{s^\top y \le \texttt{curv\_tol}\,\lVert s\rVert\lVert y\rVert}.
#'   Defaults to `1e-10`. Zero is admitted and gives \eqn{s^\top y > 0}; a
#'   negative value is refused, for the reason [bfgs()]'s page gives.
#' @param step,line_search,maxit,max_eval,verbose,refresh,keep_trace As in
#'   [bfgs()].
#'
#' @details
#' # Notation
#'
#' \eqn{s = x_{new} - x_{old}} is the secant vector and
#' \eqn{y = g_{new} - g_{old}} the change in the gradient. The scalar step
#' length is \eqn{\alpha}, a different quantity.
#'
#' # Where the saving is
#'
#' The two-loop recursion returns exactly the product the explicitly
#' assembled inverse would: built from the same four pairs and the same
#' scaling, the two agree to `1.7e-16`. What differs is the cost.
#'
#' Measured on a dense quadratic, both methods taking the gradient and
#' converging to the same point:
#'
#' \tabular{lll}{
#'   **p** \tab **bfgs** \tab **lbfgs** \cr
#'   5 \tab 9 iterations \tab 8 iterations \cr
#'   50 \tab 15 \tab 15 \cr
#'   200 \tab 16, 0.06 s \tab 17, 0.00 s \cr
#'   800 \tab 16, **4.06 s** \tab 17, **0.02 s** \cr
#' }
#'
#' The iteration counts barely differ; the cost per iteration is what
#' diverges, and it does so between 200 and 800 parameters. Beyond a few
#' thousand the full matrix is not storable at all.
#'
#' The recursion is scaled at each iteration by the most recent pair's
#' \eqn{s^\top y / y^\top y}. That single number does the work the full
#' matrix would otherwise do, and is the reason the method converges without
#' one.
#'
#' # Choosing memory
#'
#' More is not better. Old pairs describe curvature at points the iterate has
#' left, and they cost \eqn{O(p)} each per iteration. Measured on Rosenbrock,
#' `memory` of 3, 5, 10, 30 and 100 gives 38, 36, 37, 37 and 37 iterations:
#' past a handful there is nothing left to gain and the arithmetic keeps
#' rising. Ten is the conventional choice for that reason.
#'
#' @return An S7 object of class [Lbfgs], inheriting from [optimizer()], to
#'   be handed to [minimize()].
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rg <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                     200 * (p[2] - p[1]^2))
#' minimize(lbfgs(memory = 5), rosen, c(-1.2, 1), gr = rg)
#'
#' # The two-loop recursion is not an approximation to the full update: on the
#' # same pairs it returns the same product, to machine precision.
#' set.seed(1)
#' p <- 6
#' S <- matrix(rnorm(p * 4), p, 4)
#' Y <- matrix(rnorm(p * 4), p, 4) + 3 * S      # positive curvature
#' g <- rnorm(p)
#' gamma <- sum(S[, 4] * Y[, 4]) / sum(Y[, 4]^2)
#'
#' H <- diag(gamma, p)                          # the explicit inverse
#' for (j in 1:4) {
#'   s <- S[, j]; y <- Y[, j]; rho <- 1 / sum(s * y)
#'   V <- diag(p) - rho * outer(s, y)
#'   H <- V %*% H %*% t(V) + rho * outer(s, s)
#' }
#'
#' two_loop <- function(g, S, Y, gamma) {       # the recursion
#'   m <- ncol(S); a <- numeric(m); q <- g
#'   for (j in m:1) {
#'     rho <- 1 / sum(S[, j] * Y[, j])
#'     a[j] <- rho * sum(S[, j] * q); q <- q - a[j] * Y[, j]
#'   }
#'   r <- gamma * q
#'   for (j in 1:m) {
#'     rho <- 1 / sum(S[, j] * Y[, j])
#'     b <- rho * sum(Y[, j] * r); r <- r + S[, j] * (a[j] - b)
#'   }
#'   r
#' }
#' max(abs(H %*% g - two_loop(g, S, Y, gamma)))
#'
#' @seealso [bfgs()] for a few parameters, [newton()] when a Hessian is
#'   available.
#' @references
#' Nocedal, J. (1980). Updating quasi-Newton matrices with limited
#' storage. *Mathematics of Computation* **35**, 773--782.
#'
#' Liu, D. C. and Nocedal, J. (1989). On the limited memory BFGS method
#' for large scale optimization. *Mathematical Programming*
#' **45**, 503--528.
#'
#' @export
lbfgs <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                  memory = 10, curv_tol = 1e-10,
                  step = 1, line_search = wolfe(),
                  maxit = 500, max_eval = Inf,
                  verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_count(memory, "memory")
  check_nonneg(curv_tol, "curv_tol")
  Lbfgs(name = "L-BFGS", criterion = criterion, maxit = maxit,
        max_eval = max_eval, verbose = verbose, refresh = refresh,
        keep_trace = keep_trace, step = step, line_search = line_search,
        memory = memory, curv_tol = curv_tol)
}


# --- the methods ------------------------------------------------------------
#
# Each names its direction and hands everything else to run_descent(); that they
# are three lines apiece is the point of the frame.

#' @title Minimize by Newton's Method
#' @name minimize.Newton
#' @param optimizer A `Newton` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @return An [optimizer_result()].
#' @keywords internal
S7::method(minimize, Newton) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    # A numerical Hessian differences the gradient once per coordinate and
    # direction; without an analytic gradient each of those is itself 2p
    # objective values, so one iteration costs about 4p^2 evaluations. When the
    # evaluation budget admits fewer than two such iterations the run can only
    # end on the budget, so that is said here rather than discovered from a
    # one-iteration result.
    if (is.null(he) && is.function(fn)) {
      p <- length(par)
      per_iter <- if (is.null(gr)) 4 * p^2 + 2 * p else 0
      if (per_iter > 0 && 2 * per_iter > optimizer@max_eval) {
        warning("newton() without 'gr' and 'he' costs about 4*p^2 = ",
                format(4 * p^2), " objective evaluations per iteration at p = ",
                p, ",
  and max_eval = ", format(optimizer@max_eval),
                " admits fewer than two iterations. Supply 'gr' (or 'he'),
",
                "  raise 'max_eval', or use bfgs(), which needs no Hessian.",
                call. = FALSE)
      }
    }
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "newton", hessian_mod = optimizer@hessian_mod,
                     floor = optimizer@floor))
  }

#' @title Minimize by BFGS
#' @name minimize.Bfgs
#' @param optimizer A `Bfgs` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @return An [optimizer_result()].
#' @keywords internal
S7::method(minimize, Bfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "bfgs", curv_tol = optimizer@curv_tol,
                     max_skip = as.integer(optimizer@max_skip)))
  }

#' @title Minimize by Limited-Memory BFGS
#' @name minimize.Lbfgs
#' @param optimizer An `Lbfgs` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @return An [optimizer_result()].
#' @keywords internal
S7::method(minimize, Lbfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "lbfgs", memory = as.integer(optimizer@memory),
                     curv_tol = optimizer@curv_tol))
  }
