# Adam, on an exact objective and on a noisy one.
#
# The tests below deliberately do NOT check that Adam beats the second-order
# methods, because it does not and should not: on a small smooth problem it is
# the wrong tool. What they check is that it converges to the right point when
# given enough iterations, and that it finds the answer through the noise of an
# objective that resamples itself.

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


test_that("a real criterion is honored on the full sample", {
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


# --- a noisy objective, which is now the caller's business -------------------
#
# Adam no longer draws minibatches. It does not need to: an objective that
# resamples is a closure, and Adam optimizing a noisy objective is what
# subsampling ever was. These tests are the same claims as before, made through
# the interface that replaced the machinery.

test_that("Adam finds the mean through the noise of a resampling objective", {
  set.seed(11)
  y <- rnorm(1000, mean = 3)
  m <- 50
  bf <- function(p) { i <- sample.int(1000, m); sum((y[i] - p)^2) / 2 }
  bg <- function(p) { i <- sample.int(1000, m); -sum(y[i] - p) }

  r <- minimize(adam(alpha = 0.05, decay = 0.005, maxit = 4000),
                bf, par = 0, gr = bg)
  expect_equal(r@par, mean(y), tolerance = 1e-2)
})


test_that("set.seed governs it, because the draws happen in the objective", {
  set.seed(13)
  y <- rnorm(300, mean = -1)
  bf <- function(p) { i <- sample.int(300, 30); sum((y[i] - p)^2) / 2 }
  bg <- function(p) { i <- sample.int(300, 30); -sum(y[i] - p) }

  set.seed(99); a <- minimize(adam(maxit = 300), bf, par = 0, gr = bg)
  set.seed(99); b <- minimize(adam(maxit = 300), bf, par = 0, gr = bg)
  expect_identical(a@par, b@par)
})


test_that("restarting Adam each batch throws its state away", {
  # The reason ?adam tells the caller to resample INSIDE the objective. Every
  # restart puts m and v back to zero and the bias correction back to t = 1, so
  # after one iteration mhat/sqrt(vhat) is the SIGN of the gradient exactly, and
  # the step is alpha whatever the surface looks like.
  #
  # The claim is about the step lengths, not about which run ends up closer.
  # A first version of this test compared the two answers and failed, correctly:
  # in one well-scaled dimension there is no adaptivity to destroy, and both
  # runs finish within alpha of the target, so which is nearer is luck.
  al <- 0.05
  f  <- function(p) sum((p - c(1, 200))^2)     # coordinates on wildly
  g  <- function(p) 2 * (p - c(1, 200))        # different scales

  one <- minimize(adam(alpha = al, maxit = 1, keep_trace = TRUE), f,
                  c(0, 0), gr = g)
  # both coordinates move by alpha, though their gradients differ by 200-fold
  expect_equal(one@trace$step[1], al, tolerance = 1e-6)

  # a continuous run adapts: by the end its steps are no longer alpha
  many <- minimize(adam(alpha = al, maxit = 500, keep_trace = TRUE), f,
                   c(0, 0), gr = g)
  expect_false(isTRUE(all.equal(tail(many@trace$step, 1), al,
                                tolerance = 1e-6)))

  # and the restarted loop never gets past alpha per step, so 500 of them can
  # travel at most 500 * alpha -- nowhere near the second coordinate's optimum
  expect_lt(500 * al, 200)
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
                gr = quad_g, lower = c(0, 0), upper = c(5, 1))
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
  # resample is gone: an optimizer does not know what an observation is
  expect_error(adam(resample = 0.5), "unused argument")
})


test_that("crit_never never fires", {
  st <- list(iter = 1, f_new = 0, f_old = 0, x_new = 0, x_old = 0,
             gradient = 0)
  expect_false(crit_met(crit_never(), st))
  expect_identical(crit_needs(crit_never()), character())
})
