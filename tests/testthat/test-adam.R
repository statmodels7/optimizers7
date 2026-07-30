# Adam: the full-sample method, the stochastic one, and the refusals.
#
# The tests below deliberately do NOT check that Adam beats the second-order
# methods, because it does not and should not: on a small smooth problem it is
# the wrong tool. What they check is that it converges to the right point when
# given enough iterations, that the stochastic path finds the answer through the
# noise, and that everything it cannot honestly report is refused rather than
# quietly returned.

quad <- function(p) sum((p - c(1, 2))^2)
quad_g <- function(p) 2 * (p - c(1, 2))


test_that("Adam reaches the minimum of a quadratic", {
  r <- minimize(adam(alpha = 0.1, maxit = 3000), quad, c(0, 0), gr = quad_g)
  expect_equal(r@par, c(1, 2), tolerance = 1e-5)
  expect_equal(r@value, 0, tolerance = 1e-9)
})


test_that("a budget run reports converged = FALSE, because nothing checked", {
  r <- minimize(adam(alpha = 0.1, maxit = 3000), quad, c(0, 0), gr = quad_g)
  expect_false(r@converged)
  expect_match(r@criterion_met, "iteration budget")
  # ...and it is nonetheless at the optimum, which is the point: `converged`
  # records what a rule confirmed, not how good the answer happens to be.
  expect_equal(r@par, c(1, 2), tolerance = 1e-5)
})


test_that("a real criterion is honoured on the full sample", {
  r <- minimize(adam(alpha = 0.1, criterion = crit_grad(1e-7), maxit = 5000),
                quad, c(0, 0), gr = quad_g)
  expect_true(r@converged)
  expect_lt(r@iterations, 5000)
  expect_lt(max(abs(r@gradient)), 1e-6)
})


test_that("the gradient shown to the stopping rule is the one at x_new", {
  # The rule is asked at the top of the following iteration precisely so that
  # crit_grad() tests the point the run has reached, not the one it left. If it
  # were shown the stale gradient, stopping at 1e-7 would leave a gradient at
  # the reported point that need not be small at all.
  r <- minimize(adam(alpha = 0.05, criterion = crit_grad(1e-7), maxit = 20000),
                quad, c(0, 0), gr = quad_g)
  expect_true(r@converged)
  expect_lt(max(abs(quad_g(r@par))), 1e-6)
})


test_that("the first step is alpha, which pins bias correction and the floor", {
  # After one iteration mhat = g and vhat = g^2 exactly, so the step is
  # alpha |g| / (|g| + eps): almost alpha, and short of it by precisely the
  # denominator floor. Two things are pinned at once. Without bias correction
  # the ratio would be (1-b1)/sqrt(1-b2) = 0.1/0.0316, more than three times
  # too long; and asserting the eps term rather than tolerating it means the
  # floor cannot quietly disappear.
  al  <- 0.03
  eps <- 1e-8
  gmax <- max(abs(quad_g(c(0, 0))))
  r <- minimize(adam(alpha = al, eps = eps, maxit = 1, keep_trace = TRUE),
                quad, c(0, 0), gr = quad_g)
  expect_equal(r@trace$step[1], al * gmax / (gmax + eps), tolerance = 1e-14)
  expect_lt(r@trace$step[1], al)
})


test_that("amsgrad changes the path and still finds the optimum", {
  plain <- minimize(adam(alpha = 0.1, maxit = 3000, keep_trace = TRUE),
                    quad, c(0, 0), gr = quad_g)
  ams   <- minimize(adam(alpha = 0.1, maxit = 3000, amsgrad = TRUE,
                         keep_trace = TRUE), quad, c(0, 0), gr = quad_g)
  expect_equal(ams@par, c(1, 2), tolerance = 1e-4)
  expect_false(isTRUE(all.equal(plain@trace$step, ams@trace$step)))
  # the floor is reported where it bites
  expect_true(any(ams@trace$safeguard == "amsgrad floor"))
})


test_that("decay shortens the steps monotonically", {
  r <- minimize(adam(alpha = 0.1, decay = 1, maxit = 200, keep_trace = TRUE),
                quad, c(0, 0), gr = quad_g)
  late <- r@trace$step[150:200]
  expect_lt(max(late), r@trace$step[1])
})


# --- subsampling ------------------------------------------------------------

make_sum <- function(y) {
  finite_sum(fn = function(par, idx) sum((y[idx] - par)^2) / 2,
             gr = function(par, idx) -sum(y[idx] - par),
             n  = length(y))
}

test_that("minibatch Adam finds the mean through the sampling noise", {
  set.seed(11)
  y <- rnorm(1000, mean = 3)
  r <- minimize(adam(alpha = 0.05, resample = 0.05, decay = 0.005,
                     maxit = 4000),
                make_sum(y), par = 0)
  expect_equal(r@par, mean(y), tolerance = 1e-2)
})


test_that("the reported value and gradient are full-sample, not minibatch", {
  set.seed(12)
  y <- rnorm(400, mean = 1.5)
  obj <- make_sum(y)
  r <- minimize(adam(alpha = 0.05, resample = 0.1, maxit = 1500), obj, par = 0)

  full_f <- sum((y - r@par)^2) / 2
  full_g <- -sum(y - r@par)
  expect_equal(r@value, full_f, tolerance = 1e-10)
  expect_equal(r@gradient, full_g, tolerance = 1e-8)
})


test_that("the same seed gives the same stochastic run", {
  set.seed(13)
  y <- rnorm(300, mean = -1)
  obj <- make_sum(y)
  set.seed(99); a <- minimize(adam(resample = 0.1, maxit = 300), obj, par = 0)
  set.seed(99); b <- minimize(adam(resample = 0.1, maxit = 300), obj, par = 0)
  expect_identical(a@par, b@par)
})


test_that("resample below 1 refuses a plain function and a gradient-free sum", {
  expect_error(minimize(adam(resample = 0.5), quad, c(0, 0), gr = quad_g),
               "no terms to draw from")

  no_gr <- finite_sum(fn = function(par, idx) sum((seq_len(10)[idx] - par)^2),
                      n = 10)
  expect_error(minimize(adam(resample = 0.5), no_gr, par = 0),
               "analytic gradient")
})


test_that("a criterion reading a noisy quantity is refused, not left to hang", {
  set.seed(14)
  y <- rnorm(100)
  obj <- make_sum(y)

  expect_error(minimize(adam(resample = 0.2, criterion = crit_grad(1e-6)),
                        obj, par = 0),
               "gradient")
  expect_error(minimize(adam(resample = 0.2, criterion = crit_rel_obj()),
                        obj, par = 0),
               "objective")
  # the explanation, not just the refusal
  expect_error(minimize(adam(resample = 0.2, criterion = crit_grad()),
                        obj, par = 0),
               "sampling noise")

  # a rule on the parameters reads nothing noisy and is accepted
  expect_silent(minimize(adam(resample = 0.2, criterion = crit_abs_par(1e-12),
                             maxit = 20), obj, par = 0))
})


test_that("a full-sample Adam accepts every criterion", {
  expect_silent(minimize(adam(criterion = crit_rel_obj(), maxit = 20),
                         quad, c(0, 0), gr = quad_g))
  expect_silent(minimize(adam(criterion = crit_grad(), maxit = 20),
                         quad, c(0, 0), gr = quad_g))
})


# --- bounds and safeguards --------------------------------------------------

test_that("Adam respects a box, and reports on the user's scale", {
  r <- minimize(adam(alpha = 0.1, maxit = 4000), quad, c(0.5, 0.5),
                gr = quad_g, bounds = list(c(0, 5), c(0, 1)))
  expect_gt(r@par[2], 0.99)
  expect_lt(r@par[2], 1)
  expect_equal(r@par[1], 1, tolerance = 1e-3)
  expect_equal(r@gradient, quad_g(r@par), tolerance = 1e-4)
})


test_that("a non-finite gradient stops the run rather than poisoning it", {
  bad <- function(p) sum(p^2)
  bad_g <- function(p) if (p[1] < 0.5) c(NaN) else 2 * p
  r <- minimize(adam(alpha = 0.5, maxit = 500), bad, par = 1, gr = bad_g)
  expect_false(r@converged)
  expect_match(r@message, "not finite")
  expect_true(is.finite(r@par))
})


test_that("a non-finite objective at the start is an error", {
  expect_error(minimize(adam(), function(p) NaN, par = 1),
               "not finite at the starting value")
})


test_that("the constructor refuses nonsense", {
  expect_error(adam(alpha = -1), "'alpha'")
  expect_error(adam(beta1 = 1), "'beta1'")
  expect_error(adam(beta2 = -0.1), "'beta2'")
  expect_error(adam(decay = -1), "'decay'")
  expect_error(adam(amsgrad = 1), "'amsgrad'")
  expect_error(adam(resample = 0), "'resample'")
  expect_error(adam(resample = 1.5), "'resample'")
})


test_that("crit_never never fires", {
  st <- list(iter = 1, f_new = 0, f_old = 0, x_new = 0, x_old = 0,
             gradient = 0)
  expect_false(crit_met(crit_never(), st))
  expect_identical(crit_needs(crit_never()), character())
})
