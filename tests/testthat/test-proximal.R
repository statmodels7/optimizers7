# The proximal gradient method. The solution is checked against the KKT
# conditions of the problem, which share no code with the iteration, and the
# smooth special case is checked against gradient descent on the same
# objective.

lasso_problem <- function(n = 200, p = 8, lambda = 0.4, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  b0 <- c(2, -1.5, 0, 0, 0.8, rep(0, max(0, p - 5)))[seq_len(p)]
  y <- as.numeric(X %*% b0 + rnorm(n))
  list(X = X, y = y, lambda = lambda, n = n, p = p,
       fn = function(b) sum((y - X %*% b)^2) / (2 * n),
       gr = function(b) -as.numeric(crossprod(X, y - X %*% b)) / n,
       prox = function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0),
       g = function(b) lambda * sum(abs(b)))
}

# max violation of the lasso stationarity conditions: the smooth gradient
# must equal -lambda*sign(b) where b is non-zero and lie inside the interval
# elsewhere.
kkt_violation <- function(pr, b, tol = 1e-8) {
  gsm <- pr$gr(b)
  on <- abs(b) > tol
  v_on <- if (any(on)) max(abs(gsm[on] + pr$lambda * sign(b[on]))) else 0
  v_off <- if (any(!on)) max(pmax(abs(gsm[!on]) - pr$lambda, 0)) else 0
  max(v_on, v_off)
}

test_that("the lasso solution satisfies its own stationarity conditions", {
  # The tolerance is one the run can reach: with momentum the attainable
  # mapping is bounded by the rounding of the objective, measured at ~4e-9
  # on this problem, so a tighter ask would report failure at the answer.
  pr <- lasso_problem()
  fit <- minimize(prox_grad(prox = pr$prox, g = pr$g, criterion = crit_grad(1e-8)),
                  fn = pr$fn, gr = pr$gr, par = rep(0, pr$p))
  expect_true(fit@converged)
  expect_lt(kkt_violation(pr, fit@par), 1e-7)
  # the reported mapping IS the stationarity measure: it agrees with the KKT
  # violation, computed here from the conditions rather than from the iteration
  expect_equal(max(abs(fit@gradient)), kkt_violation(pr, fit@par),
               tolerance = 1e-6)
  # the value reported is the total objective, not the smooth part
  expect_equal(fit@value, pr$fn(fit@par) + pr$g(fit@par))
  # and the zeros are exact, which is what the operator buys over a smooth
  # approximation of the absolute value
  expect_true(sum(fit@par == 0) >= 3)
})

test_that("acceleration pays where the problem is ill conditioned", {
  # Measured across conditioning: at a condition number of 3 the plain
  # method wins narrowly (39 iterations against 24 for the accelerated one,
  # which pays an extra gradient per iteration for its mapping), by 55 it is
  # 4153 against 126, and by 480 the plain method does not converge at all.
  # The claim being tested is therefore about the ill-conditioned case, and a
  # well-conditioned problem would confirm the opposite.
  set.seed(3)
  n <- 400; p <- 30; rho <- 0.5
  Z <- matrix(rnorm(n * p), n, p)
  X <- sqrt(1 - rho) * Z + sqrt(rho) * matrix(rnorm(n), n, p)
  b0 <- c(2, -1.5, rep(0, p - 3), 0.8)
  y <- as.numeric(X %*% b0 + rnorm(n))
  lambda <- 0.05
  fn <- function(b) sum((y - X %*% b)^2) / (2 * n)
  gr <- function(b) -as.numeric(crossprod(X, y - X %*% b)) / n
  prox <- function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0)
  gv <- function(b) lambda * sum(abs(b))

  slow <- minimize(prox_grad(prox, gv, accelerate = FALSE,
                             criterion = crit_grad(1e-8), maxit = 20000),
                   fn = fn, gr = gr, par = rep(0, p))
  fast <- minimize(prox_grad(prox, gv, accelerate = TRUE,
                             criterion = crit_grad(1e-8), maxit = 20000),
                   fn = fn, gr = gr, par = rep(0, p))
  expect_true(slow@converged && fast@converged)
  expect_equal(fast@par, slow@par, tolerance = 1e-4)

  # An iteration count is a DISCONTINUOUS function of the data, so one draw
  # cannot carry this claim. The restart is a discrete decision taken on a
  # continuous quantity, and the plain method's stopping index is the crossing
  # time of a nearly flat sequence: perturbing the response by a few ulps moves
  # the plain count between 753 and 4153 at this seed alone, and across seeds
  # it runs from 598 to 11964. Asserted on one draw, fast < slow/2 asserts the
  # draw. Measured over twelve it fails on two of them here, and macOS is in
  # effect a thirteenth: it reached 596 against a plain 652 and reddened this
  # line while the other four platforms passed.
  #
  # What survives every draw is the claim in the median, which is what
  # "acceleration pays" means. Over nine the ratios are 1.0, 7.8, 32.7, 6.6,
  # 16.4, 1.4, 4.8, 6.6 and 7.7, so the median is 6.6 against a bound of 2 and
  # five of the nine would have to collapse before this failed. They are
  # printed so a failure elsewhere can be read rather than guessed at.
  ratio <- function(seed) {
    set.seed(seed)
    Z <- matrix(rnorm(n * p), n, p)
    X <- sqrt(1 - rho) * Z + sqrt(rho) * matrix(rnorm(n), n, p)
    y <- as.numeric(X %*% b0 + rnorm(n))
    fn <- function(b) sum((y - X %*% b)^2) / (2 * n)
    gr <- function(b) -as.numeric(crossprod(X, y - X %*% b)) / n
    count <- function(accel) {
      minimize(prox_grad(prox, gv, accelerate = accel,
                         criterion = crit_grad(1e-8), maxit = 20000),
               fn = fn, gr = gr, par = rep(0, p))@iterations
    }
    count(FALSE) / count(TRUE)
  }
  ratios <- vapply(1:9, ratio, numeric(1))
  report <- sprintf("plain/accelerated over nine draws: %s",
                    paste(sprintf("%.1f", ratios), collapse = ", "))
  expect_gt(median(ratios), 2, label = report)
})

test_that("the restart is not fired by the objective's own rounding", {
  # Near the solution the objective moves at the rounding level, so a restart
  # test reading a bare > discards the momentum over and over and the run
  # creeps. The allowance is eight units in the last place of the current
  # objective.
  set.seed(1)
  n <- 200L
  p <- 8L
  X <- matrix(rnorm(n * p), n, p)
  b0 <- c(2, -1.5, 0, 0, 0.8, 0, 0, 0)
  y <- as.numeric(X %*% b0 + rnorm(n))
  lambda <- 0.4
  fn <- function(b) sum((y - X %*% b)^2) / (2 * n)
  gr <- function(b) -as.numeric(crossprod(X, y - X %*% b)) / n
  prox <- function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0)
  gv <- function(b) lambda * sum(abs(b))

  run <- function(...) {
    minimize(prox_grad(prox, gv, criterion = crit_grad(1e-9),
                       maxit = 100000L, ...),
             fn = fn, gr = gr, par = rep(0, p))
  }
  kept <- run()
  off <- run(restart = FALSE)
  plain <- run(accelerate = FALSE)

  # all three reach the same answer
  expect_true(kept@converged && off@converged && plain@converged)
  expect_equal(kept@par, plain@par, tolerance = 1e-7)
  expect_identical(sum(kept@par != 0), sum(plain@par != 0))

  # and the guarded restart is nowhere near the count the unguarded one took,
  # which was 21646 on this problem. The bound is loose because an iteration
  # count on a tight tolerance is platform arithmetic, and the counts are
  # printed so a failure elsewhere can be read.
  report <- sprintf("restart %d, restart=FALSE %d, plain %d",
                    kept@iterations, off@iterations, plain@iterations)
  expect_lt(kept@iterations, 5000L, label = report)

  # the default tolerance never saw the defect and must not move
  quick <- minimize(prox_grad(prox, gv, criterion = crit_grad(1e-6)),
                    fn = fn, gr = gr, par = rep(0, p))
  expect_lt(quick@iterations, 30L)
})

test_that("a real increase still restarts", {
  # The negative control. On an ill-conditioned quadratic the increases the
  # restart exists for are far above the objective's rounding, so the
  # allowance must not suppress them: the accelerated run with the restart
  # has to stay far ahead of the one without it.
  p <- 8L
  d <- exp(seq(0, log(2400), length.out = p))
  set.seed(7)
  b <- rnorm(p)
  fn <- function(x) sum(d * (x - b)^2) / 2
  gr <- function(x) d * (x - b)
  run <- function(res) {
    minimize(prox_grad(function(v, t) v, function(z) 0, restart = res,
                       criterion = crit_grad(1e-8), maxit = 200000L),
             fn = fn, gr = gr, par = rep(0, p))
  }
  on_ <- run(TRUE)
  off <- run(FALSE)
  expect_true(on_@converged && off@converged)
  expect_equal(on_@par, b, tolerance = 1e-5)
  report <- sprintf("restart %d, restart=FALSE %d",
                    on_@iterations, off@iterations)
  expect_lt(on_@iterations, off@iterations / 4, label = report)
})

test_that("with no non-smooth part the method is gradient descent", {
  # prox of zero is the identity, so the iteration reduces to the smooth one
  # and must land where a smooth method lands
  fn <- function(p) sum((p - c(1, -2, 0.5))^2) + 0.3 * sum(p^4)
  gr <- function(p) 2 * (p - c(1, -2, 0.5)) + 1.2 * p^3
  pg <- minimize(prox_grad(function(v, t) v, function(b) 0,
                           criterion = crit_grad(1e-8), maxit = 5000),
                 fn = fn, gr = gr, par = c(0, 0, 0))
  ref <- minimize(bfgs(criterion = crit_grad(1e-10)), fn = fn, gr = gr,
                  par = c(0, 0, 0))
  expect_equal(pg@par, ref@par, tolerance = 1e-6)
  # the reported mapping is then the ordinary gradient
  expect_equal(pg@gradient, gr(pg@par), tolerance = 1e-6)
})

test_that("the method passes the optimizer contract", {
  res <- check_optimizer(
    prox_grad(function(v, t) v, function(b) 0),
    verbose = FALSE)
  expect_true(all(res$checks), info = paste(names(res$checks)[!res$checks],
                                            collapse = ", "))
})

test_that("a run without a supplied gradient differences the smooth part", {
  pr <- lasso_problem(n = 100, p = 4)
  fit <- minimize(prox_grad(pr$prox, pr$g, criterion = crit_grad(1e-7)),
                  fn = pr$fn, par = rep(0, 4))
  expect_lt(kkt_violation(pr, fit@par), 1e-5)
  expect_match(fit@message, "finite differences")
})

test_that("the trace records the iteration path", {
  pr <- lasso_problem(n = 80, p = 5)
  fit <- minimize(prox_grad(pr$prox, pr$g, keep_trace = TRUE, maxit = 25),
                  fn = pr$fn, gr = pr$gr, par = rep(0, 5))
  expect_s3_class(fit@trace, "data.frame")
  expect_named(fit@trace, c("iteration", "value", "step", "gradient"))
  expect_true(all(diff(fit@trace$iteration) == 1))
})

test_that("the constructor and the method reject what they cannot honour", {
  expect_error(prox_grad(prox = "not a function", g = function(b) 0),
               "'prox' must be a function")
  expect_error(prox_grad(prox = function(v, t) v), "'g' must be a function")
  expect_error(prox_grad(function(v, t) v, function(b) 0, shrink = 1),
               "strictly between 0 and 1")
  expect_error(prox_grad(function(v, t) v, function(b) 0, accelerate = NA),
               "TRUE or FALSE")
  # box bounds belong inside the operator, not beside it
  expect_error(
    minimize(prox_grad(function(v, t) v, function(b) 0),
             fn = function(p) sum(p^2), gr = function(p) 2 * p,
             par = c(1, 1), lower = 0),
    "through the proximal operator")
})
