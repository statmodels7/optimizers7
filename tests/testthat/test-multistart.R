# Multi-start, and the one thing it can tell you that a single run cannot.

# Two minima: at x = -1 and x = +1, tilted by 0.3x so the left one is lower.
# f(-1) = -0.3, f(+1) = 0.3.
twowell <- function(p) (p[1]^2 - 1)^2 + p[2]^2 + 0.3 * p[1]
twowell_g <- function(p) c(4 * p[1] * (p[1]^2 - 1) + 0.3, 2 * p[2])

quad <- function(p) sum((p - c(1, 2))^2)


test_that("it finds the lower of two minima that a single run can miss", {
  # From (0.5, 0.5) one run settles in the shallower well at x = +0.96 and has
  # no way of knowing there is a better one; the tilt makes the left well lower.
  single <- minimize(bfgs(), twowell, c(0.5, 0.5), gr = twowell_g)
  expect_gt(single@par[1], 0)

  set.seed(1)
  many <- minimize(multistart(bfgs(), n = 15), twowell, c(0.5, 0.5),
                   gr = twowell_g)
  expect_lt(many@par[1], 0)
  expect_lt(many@value, single@value)
})


test_that("it counts the distinct optima, which is the point of running it", {
  set.seed(2)
  r <- minimize(multistart(bfgs(), n = 20), twowell, c(1.5, 0.5),
                gr = twowell_g)
  expect_match(r@message, "2 distinct optima")
  expect_match(r@message, "20 starts")

  # and on a problem with one minimum it says one
  set.seed(2)
  r1 <- minimize(multistart(bfgs(), n = 8), quad, c(0, 0))
  expect_match(r1@message, "1 distinct optimum")
})


test_that("the user's own start is always among them, and is start 1", {
  set.seed(3)
  r <- minimize(multistart(bfgs(), n = 5, keep_trace = TRUE), quad, c(0, 0))
  expect_equal(nrow(r@trace), 5L)
  # started at the answer, start 1 converges immediately
  r0 <- minimize(multistart(bfgs(), n = 5, keep_trace = TRUE), quad, c(1, 2))
  expect_equal(r0@trace$iterations[1], 1L)
})


test_that("the per-start summary is the trace", {
  set.seed(4)
  r <- minimize(multistart(bfgs(), n = 6), twowell, c(0.5, 0.5),
                gr = twowell_g)
  expect_named(r@trace, c("start", "value", "converged", "iterations"))
  expect_equal(nrow(r@trace), 6L)
  expect_equal(min(r@trace$value, na.rm = TRUE), r@value)
})


test_that("random starts respect bounds without a single rejected draw", {
  # Generated on the unconstrained scale and mapped back, so a start for a
  # positive parameter cannot come out negative.
  set.seed(5)
  y <- rnorm(200, mean = 0, sd = 2)
  nll <- function(p) length(y) * log(p[1]) + sum(y^2) / (2 * p[1]^2)
  r <- minimize(multistart(bfgs(), n = 10, spread = 2), nll, par = 1,
                bounds = list(c(0, Inf)))
  expect_true(r@converged)
  expect_equal(r@par, sqrt(mean(y^2)), tolerance = 1e-6)
  expect_gt(r@par, 0)
})


test_that("a start that fails does not take the run down with it", {
  # Undefined for a negative first coordinate, so several random starts land
  # where the objective is not finite.
  f <- function(p) if (p[1] <= 0) stop("undefined here") else (p[1] - 2)^2
  set.seed(6)
  r <- minimize(multistart(nelder_mead(), n = 12, spread = 3), f, par = 1)
  expect_equal(r@par, 2, tolerance = 1e-5)
  expect_true(any(is.na(r@trace$value)))
  expect_match(r@message, "succeeded")
})


test_that("but every start failing is an error", {
  f <- function(p) stop("always broken")
  set.seed(7)
  expect_error(minimize(multistart(bfgs(), n = 3), f, par = 1),
               "Every start failed")
})


test_that("an explicit matrix of starts is used verbatim", {
  S <- rbind(c(-1.5, 0), c(1.5, 0), c(0.2, 0.2))
  r <- minimize(multistart(bfgs(), starts = S, keep_trace = TRUE),
                twowell, c(0, 0), gr = twowell_g)
  expect_equal(nrow(r@trace), 3L)
  expect_lt(r@par[1], 0)
  expect_error(minimize(multistart(bfgs(), starts = matrix(0, 3, 5)),
                        twowell, c(0, 0)),
               "2 columns")
})


test_that("a Latin hypercube uses each stratum once", {
  # The property worth having over independent draws: the starts cannot all
  # cluster in one part of the range. With m = n - 1 random starts, one falls in
  # each of the m equal strata of the sampling interval.
  set.seed(8)
  S <- make_starts(par = 0, n = 21, spread = 1, bounds = list())
  x <- S[-1, 1]
  half <- 3                                   # 3 * spread * max(1, |0|)
  strata <- floor((x + half) / (2 * half) * 20) + 1
  expect_setequal(strata, 1:20)
})


test_that("the whole thing composes: multistart of a derivative-free method", {
  set.seed(9)
  r <- minimize(multistart(nelder_mead(), n = 8), twowell, c(1.5, 0.5))
  expect_lt(r@par[1], 0)
  expect_null(r@gradient)
})


test_that("a criterion the inner method cannot evaluate is refused once", {
  expect_error(minimize(multistart(nelder_mead(criterion = crit_grad()), n = 3),
                        quad, c(0, 0)),
               "needs gradient")
  expect_identical(optimizer_provides(multistart(nelder_mead())),
                   optimizer_provides(nelder_mead()))
})


test_that("the run records the seed it began with, and it reproduces it", {
  set.seed(10)
  a <- minimize(multistart(bfgs(), n = 6), twowell, c(1.5, 0.5),
                gr = twowell_g)
  expect_false(is.null(a@seed))

  # putting the recorded state back repeats the run exactly
  assign(".Random.seed", a@seed, envir = globalenv())
  b <- minimize(multistart(bfgs(), n = 6), twowell, c(1.5, 0.5),
                gr = twowell_g)
  expect_identical(a@trace$value, b@trace$value)
})


test_that("only the methods that draw record a seed", {
  expect_null(minimize(bfgs(), quad, c(0, 0))@seed)
  expect_null(minimize(compass(directions = "coordinate"), quad, c(0, 0))@seed)
  expect_false(is.null(minimize(compass(directions = "mads"), quad, c(0, 0))@seed))

  # Adam draws nothing now, so it records nothing. A noisy objective draws in
  # the caller's own code, where the caller's own set.seed() governs it.
  expect_null(minimize(adam(maxit = 5), quad, c(0, 0))@seed)
})


test_that("a mads run repeats from its own recorded seed", {
  set.seed(11)
  a <- minimize(compass(directions = "mads"), twowell, c(1.5, 0.5))
  assign(".Random.seed", a@seed, envir = globalenv())
  b <- minimize(compass(directions = "mads"), twowell, c(1.5, 0.5))
  expect_identical(a@par, b@par)
})


test_that("the constructor refuses nonsense", {
  expect_error(multistart("bfgs"), "must be an optimizer")
  expect_error(multistart(bfgs(), n = 0), "'maxit'")
  expect_error(multistart(bfgs(), spread = -1), "'tol'")
  expect_error(multistart(bfgs(), starts = 1:3), "one starting point per row")
})


test_that("it prints the inner optimiser by name", {
  out <- capture.output(print(multistart(bfgs(), n = 4)))
  expect_true(any(grepl("multistart", out)))
  expect_true(any(grepl("optimizer = BFGS", out)))
})
