# Newton, BFGS and L-BFGS. What is tested here is not that they find minima --
# any of them would on a quadratic -- but the safeguards, which is where a
# second-order method is actually difficult.

rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                          200 * (p[2] - p[1]^2))
rosen_he <- function(p) matrix(c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
                                 -400 * p[1], 200), 2, 2)

# x^4 - 2x^2 + y^2: the Hessian is diag(12x^2 - 4, 2), indefinite near x = 0,
# and the function is bounded below with minima at (+-1, 0) and value -1.
saddle    <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
saddle_gr <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
saddle_he <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)


test_that("all three solve Rosenbrock, and far faster than first order", {
  x0 <- c(-1.2, 1)
  n <- minimize(newton(), rosen, x0, gr = rosen_gr, he = rosen_he)
  b <- minimize(bfgs(),   rosen, x0, gr = rosen_gr)
  l <- minimize(lbfgs(),  rosen, x0, gr = rosen_gr)

  for (r in list(n, b, l)) {
    expect_true(r@converged)
    expect_equal(r@par, c(1, 1), tolerance = 1e-6)
  }

  # the whole reason for using curvature: gradient descent needs thousands
  gd <- minimize(gd(maxit = 50000, max_eval = 1e7), rosen, x0,
                 gr = rosen_gr)
  expect_true(gd@converged)
  expect_lt(n@iterations, 100)
  expect_lt(b@iterations, 100)
  expect_lt(l@iterations, 100)
  expect_gt(gd@iterations, 1000)
})


test_that("Newton repairs an indefinite Hessian and says that it did", {
  # At the starting point the Hessian has a negative eigenvalue, so the
  # unmodified Newton direction points away from the minimum. Both repairs must
  # reach a true minimizer, and both must record what they did.
  expect_lt(min(eigen(saddle_he(c(0.05, 0.6)))$values), 0)

  for (mod in c("eigen", "ridge")) {
    r <- minimize(newton(hessian_mod = mod, keep_trace = TRUE),
                  saddle, c(0.05, 0.6), gr = saddle_gr, he = saddle_he)
    expect_true(r@converged, label = mod)
    expect_equal(r@value, -1, tolerance = 1e-8, label = mod)
    expect_equal(abs(r@par[1]), 1, tolerance = 1e-6, label = mod)
    expect_true("hessian modified" %in% r@trace$safeguard, label = mod)
  }
})


test_that("Newton works with the Hessian differenced from the gradient", {
  # One numerical differentiation of an analytic gradient, which is the case
  # newton() is willing to stand behind.
  with_he <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
  no_he   <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr)

  expect_true(no_he@converged)
  expect_equal(no_he@par, c(1, 1), tolerance = 1e-6)
  expect_equal(no_he@counts[["h"]], 0)      # none supplied, so none called
  expect_gt(with_he@counts[["h"]], 0)
})


test_that("BFGS and L-BFGS never ask for a Hessian", {
  b <- minimize(bfgs(),  rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
  l <- minimize(lbfgs(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
  # he is accepted and ignored, so calling code need not branch on the method
  expect_equal(b@counts[["h"]], 0)
  expect_equal(l@counts[["h"]], 0)
  expect_true(b@converged)
  expect_true(l@converged)
})


test_that("the objective never increases, for any of the three", {
  for (o in list(newton(keep_trace = TRUE), bfgs(keep_trace = TRUE),
                 lbfgs(keep_trace = TRUE))) {
    r <- minimize(o, rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
    expect_true(all(diff(r@trace$value) <= 1e-12), label = r@optimizer@name)
  }
})


test_that("L-BFGS with more memory is not worse, and matches BFGS", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(30 * 40), 40, 30)) + diag(30)
  b <- rnorm(30)
  f  <- function(p) as.numeric(0.5 * t(p) %*% A %*% p - sum(b * p))
  g  <- function(p) as.numeric(A %*% p - b)
  truth <- solve(A, b)

  small <- minimize(lbfgs(memory = 3, maxit = 5000),  f, rep(0, 30), gr = g)
  large <- minimize(lbfgs(memory = 30, maxit = 5000), f, rep(0, 30), gr = g)
  full  <- minimize(bfgs(maxit = 5000), f, rep(0, 30), gr = g)

  for (r in list(small, large, full)) {
    expect_true(r@converged)
    expect_equal(r@par, truth, tolerance = 1e-5)
  }
  # more curvature remembered, no more iterations needed
  expect_lte(large@iterations, small@iterations)
})


test_that("a non-finite region is stepped around by every method", {
  f  <- function(p) if (sum(p^2) > 4) Inf else sum((p - c(1, 1))^2)
  g  <- function(p) if (sum(p^2) > 4) c(NaN, NaN) else 2 * (p - c(1, 1))
  for (o in list(newton(step = 50), bfgs(step = 50), lbfgs(step = 50))) {
    r <- minimize(o, f, c(0, 0), gr = g)
    expect_true(is.finite(r@value), label = r@optimizer@name)
    expect_equal(r@par, c(1, 1), tolerance = 1e-4, label = r@optimizer@name)
  }
})


test_that("the three constructors validate their own settings", {
  expect_error(newton(hessian_mod = "nonsense"), "should be one of")
  expect_error(newton(floor = 0), "positive")
  expect_error(newton(step = -1), "positive")
  expect_error(bfgs(line_search = "wolfe"), "line_search object")
  expect_error(lbfgs(memory = 0), "positive")
  expect_error(lbfgs(criterion = "gradient"), "criterion object")
})


test_that("each carries its own settings and prints them", {
  expect_match(capture.output(print(newton()))[1], "Newton")
  expect_match(paste(capture.output(print(newton())), collapse = " "),
               "hessian_mod = eigen")
  expect_match(paste(capture.output(print(bfgs())), collapse = " "),
               "Wolfe")
  expect_match(paste(capture.output(print(lbfgs(memory = 7))), collapse = " "),
               "memory = 7")
})


test_that("every method reaches the same minimizer of a shared problem", {
  # A stronger statement than each being right on its own: they must agree.
  f  <- function(p) sum((p - c(3, -2, 0.5))^2) + 0.1 * sum(p^4)
  g  <- function(p) 2 * (p - c(3, -2, 0.5)) + 0.4 * p^3
  x0 <- c(0, 0, 0)

  rs <- list(
    newton = minimize(newton(), f, x0, gr = g),
    bfgs   = minimize(bfgs(),   f, x0, gr = g),
    lbfgs  = minimize(lbfgs(),  f, x0, gr = g),
    gd     = minimize(gd(maxit = 20000), f, x0, gr = g)
  )
  pars <- do.call(rbind, lapply(rs, function(r) r@par))
  for (j in seq_len(ncol(pars))) {
    expect_lt(diff(range(pars[, j])), 1e-5)
  }
})
