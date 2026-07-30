# Box bounds by reparametrisation.
#
# The transforms are linkfunctions7's bounded_link(), written out in C++ because
# they are evaluated on every objective evaluation and a callback into R there
# would destroy the reason for compiling the loop at all. That duplication is
# the risk, so the first test below pins the C++ implementation to
# linkfunctions7 as the reference: if the two ever drift, this fails. It is not
# hypothetical -- linkfunctions7's exponential floor changed on 2026-07-30, and
# a copy taken before then would now disagree.

skip_if_no_lf7 <- function() {
  testthat::skip_if_not_installed("linkfunctions7")
}

# The transforms as the package applies them, recovered through a bounded run of
# zero iterations is not possible, so they are exercised through the round trip
# that minimize() performs on the starting value and on the reported optimum.

test_that("the transforms agree with linkfunctions7 to machine precision", {
  skip_if_no_lf7()

  cases <- list(
    lower  = list(b = c(0, Inf),  lk = linkfunctions7::bounded_link(lwr = 0)),
    upper  = list(b = c(-Inf, 5), lk = linkfunctions7::bounded_link(upr = 5)),
    both   = list(b = c(-2, 3),   lk = linkfunctions7::bounded_link(-2, 3)),
    shift  = list(b = c(1.5, Inf), lk = linkfunctions7::bounded_link(lwr = 1.5))
  )

  for (nm in names(cases)) {
    b  <- cases[[nm]]$b
    lk <- cases[[nm]]$lk
    eta <- seq(-6, 6, length.out = 41)

    # h, h' and h'' as linkfunctions7 computes them
    ref_h  <- linkfunctions7::linkinv(lk, eta)
    ref_d1 <- linkfunctions7::dlinkinv(lk, eta)
    ref_d2 <- linkfunctions7::d2linkinv(lk, eta)

    got <- bounded_transform(b, eta)

    expect_equal(got$h,  ref_h,  tolerance = 1e-12, label = paste(nm, "h"))
    expect_equal(got$d1, ref_d1, tolerance = 1e-12, label = paste(nm, "h'"))
    expect_equal(got$d2, ref_d2, tolerance = 1e-12, label = paste(nm, "h''"))
  }
})


test_that("the forward map inverts the backward one", {
  for (b in list(c(0, Inf), c(-Inf, 5), c(-2, 3))) {
    theta <- switch(paste(b, collapse = ","),
                    "0,Inf"  = c(0.1, 1, 10),
                    "-Inf,5" = c(-3, 0, 4.9),
                    c(-1.5, 0, 2.5))
    rt <- bounded_transform(b, bounded_forward(b, theta))$h
    expect_equal(rt, theta, tolerance = 1e-10)
  }
})


test_that("a bounded run stays inside its box and finds the right point", {
  # unconstrained minimum at (1, 2); the box pushes the second coordinate
  # against its ceiling of 1 and leaves the first alone
  f <- function(p) sum((p - c(1, 2))^2)
  g <- function(p) 2 * (p - c(1, 2))
  bx <- list(c(0, 5), c(0, 1))

  for (o in list(gradient_descent(maxit = 5000), bfgs(), lbfgs(), newton())) {
    r <- minimize(o, f, c(0.5, 0.5), gr = g, bounds = bx)
    expect_gt(r@par[1], 0); expect_lt(r@par[1], 5)
    expect_gt(r@par[2], 0); expect_lt(r@par[2], 1)
    # the free coordinate reaches its unconstrained optimum
    expect_equal(r@par[1], 1, tolerance = 1e-4, label = r@optimizer@name)
    # the constrained one is pressed against the bound
    expect_gt(r@par[2], 0.99, label = r@optimizer@name)
  }
})


test_that("an interior optimum is found exactly, bounds or not", {
  # When the box does not bind, the answer must be the unconstrained one.
  f <- function(p) sum((p - c(1, 2))^2)
  g <- function(p) 2 * (p - c(1, 2))

  free  <- minimize(bfgs(), f, c(0.5, 0.5), gr = g)
  boxed <- minimize(bfgs(), f, c(0.5, 0.5), gr = g,
                    bounds = list(c(-10, 10), c(-10, 10)))

  expect_equal(boxed@par, c(1, 2), tolerance = 1e-6)
  expect_equal(boxed@par, free@par, tolerance = 1e-5)
  expect_equal(boxed@value, 0, tolerance = 1e-10)
})


test_that("par and gradient are reported on the user's scale together", {
  # The optimiser works in eta; reporting par in theta and the gradient in eta
  # would give a result whose two halves refer to different things.
  f <- function(p) sum((p - c(1, 2))^2)
  g <- function(p) 2 * (p - c(1, 2))
  r <- minimize(bfgs(), f, c(0.5, 0.5), gr = g,
                bounds = list(c(-10, 10), c(-10, 10)))

  expect_equal(r@par, c(1, 2), tolerance = 1e-6)
  # the reported gradient is the gradient of the ORIGINAL objective at par
  expect_equal(r@gradient, g(r@par), tolerance = 1e-5)
})


test_that("a one-sided box keeps a positive parameter positive", {
  # The case the toolkit actually exists for: a scale parameter.
  set.seed(3)
  y <- rnorm(200, mean = 0, sd = 2)
  nll <- function(p) length(y) * log(p[1]) + sum(y^2) / (2 * p[1]^2)
  r <- minimize(bfgs(), nll, par = 1, bounds = list(c(0, Inf)))

  expect_true(r@converged)
  expect_gt(r@par, 0)
  expect_equal(r@par, sqrt(mean(y^2)), tolerance = 1e-5)
})


test_that("a start on the boundary is refused, naming the coordinate", {
  # The transform sends a bound to an infinite eta, so a run started exactly on
  # one begins at infinity and fails far from the cause.
  f <- function(p) sum(p^2)
  expect_error(minimize(bfgs(), f, c(0, 1), bounds = list(c(0, 5), c(0, 5))),
               "parameter 1")
  expect_error(minimize(bfgs(), f, c(1, 5), bounds = list(c(0, 5), c(0, 5))),
               "parameter 2")
})


test_that("malformed bounds are refused", {
  f <- function(p) sum(p^2)
  expect_error(minimize(bfgs(), f, c(1, 1), bounds = list(c(0, 5))),
               "one element per parameter")
  expect_error(minimize(bfgs(), f, c(1, 1), bounds = list(c(0, 5), c(5, 0))),
               "strictly below")
  expect_error(minimize(bfgs(), f, c(1, 1), bounds = list(c(0, 5), 3)),
               "c\\(lower, upper\\)")
  expect_error(minimize(bfgs(), f, c(1, 1), bounds = "positive"),
               "must be a list")
})


test_that("maximize() carries bounds through", {
  f <- function(p) -sum((p - c(1, 2))^2)
  r <- maximize(bfgs(), f, c(0.5, 0.5), bounds = list(c(0, 5), c(0, 1)))
  expect_gt(r@par[2], 0.99)
  expect_lt(r@par[2], 1)
  expect_equal(r@par[1], 1, tolerance = 1e-4)
})
