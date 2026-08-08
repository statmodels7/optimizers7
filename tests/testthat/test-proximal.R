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
  expect_lt(fast@iterations, slow@iterations / 10)
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
