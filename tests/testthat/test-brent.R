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

test_that("brent() reports its own stop under its own label", {
  f <- function(x) (x - 1)^2
  # a rule that never fires: the run ends on Brent's bracket test, which is
  # reported as such and not as the rule
  r <- minimize(brent(criterion = crit_never()), f, par = 0)
  expect_true(r@converged)
  expect_identical(r@criterion_met, "bracket at machine resolution")
  r <- minimize(brent(criterion = crit_stationary(1e-14)), f, par = 0)
  expect_identical(r@criterion_met, "bracket at machine resolution")
  r <- minimize(brent(), f, par = 0)
  expect_identical(r@criterion_met, brent()@criterion@label)
})

test_that("brent() passes the previous iterate and the iteration to the rule", {
  seen <- new.env()
  Rec <- S7::new_class("Rec", parent = criterion)
  S7::method(crit_met, Rec) <- function(criterion, state) {
    seen$names <- union(seen$names, names(state)[!vapply(state, is.null, TRUE)])
    FALSE
  }
  minimize(brent(criterion = Rec(label = "rec")), function(x) (x - 1)^2, 0)
  expect_true(all(c("iter", "f_new", "f_old", "x_new", "x_old",
                    "stationarity") %in% seen$names))
  r <- minimize(brent(criterion = crit_abs_par(1e-3)), function(x) (x - 1)^2, 0)
  expect_identical(r@criterion_met, crit_abs_par(1e-3)@label)
})

test_that("brent() does not reduce a bracket it never found", {
  r <- minimize(brent(max_expand = 3), function(x) (x - 1000)^2, par = 0)
  expect_false(r@converged)
  expect_match(r@message, "found no minimum")
  expect_identical(r@iterations, 0L)
  # the point returned is the lowest one evaluated
  expect_equal(r@value, (r@par - 1000)^2)
  r <- minimize(brent(max_eval = 4), function(x) -x, par = 0)
  expect_false(r@converged)
  expect_identical(r@criterion_met, "evaluation budget reached")
})
