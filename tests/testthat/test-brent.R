test_that("brent() reaches the minimum of a smooth function of one variable", {
  f <- function(x) (x - 2)^2 + exp(-x)
  ref <- stats::optimize(f, c(-5, 10), tol = 1e-12)$minimum
  r <- minimize(brent(), f, par = 0)
  expect_true(r@converged)
  expect_equal(r@par, ref, tolerance = 1e-5)
  expect_null(r@gradient)
  # fewer evaluations than the simplex on the same function
  expect_lt(r@counts[["f"]], minimize(nelder_mead(), f, par = 0)@counts[["f"]])
})

test_that("brent() brackets a minimum far from the start, in either direction", {
  f <- function(x) (x - 40)^2
  expect_equal(minimize(brent(), f, par = 0)@par, 40, tolerance = 1e-5)
  g <- function(x) (x + 40)^2
  expect_equal(minimize(brent(), g, par = 0)@par, -40, tolerance = 1e-5)
})

test_that("brent() respects bounds and rejects more than one parameter", {
  f <- function(x) (x - 3)^2
  r <- minimize(brent(), f, par = 0.5, lower = 0, upper = 1)
  expect_true(r@par > 0 && r@par < 1)
  expect_gt(r@par, 0.99)
  expect_error(minimize(brent(), function(x) sum(x^2), par = c(1, 1)),
               "one variable")
  expect_error(minimize(brent(criterion = crit_grad()), f, par = 0), "gradient")
  expect_error(brent(step = -1), "'step'")
})

test_that("brent() records its steps in the trace", {
  r <- minimize(brent(keep_trace = TRUE), function(x) cos(x), par = 2)
  expect_equal(r@par, pi, tolerance = 1e-5)
  expect_true(all(r@trace$safeguard %in% c("parabolic", "golden")))
  expect_true(all(diff(r@trace$stationarity) <= 1e-12))
})
