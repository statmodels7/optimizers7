# Starting values the caller did not have to write out.

quad <- function(p) sum((p - c(1, 2, 3))^2)
quad_g <- function(p) 2 * (p - c(1, 2, 3))


test_that("a starter with npar is resolved to a vector of that length", {
  expect_equal(starting_values(start_zeros(4), 4), rep(0, 4))

  set.seed(1)
  v <- starting_values(start_runif(-2, 2, npar = 5), 5)
  expect_length(v, 5)
  expect_true(all(v > -2 & v < 2))

  # min and max may be given per parameter
  set.seed(1)
  w <- starting_values(start_runif(c(0, 10), c(1, 11), npar = 2), 2)
  expect_true(w[1] > 0 && w[1] < 1)
  expect_true(w[2] > 10 && w[2] < 11)
})


test_that("zeros are zeros on the UNCONSTRAINED scale, which is the point", {
  # A vector of zeros on the parameter scale is not a starting value at all for
  # a positive parameter: it sits on the boundary, where the objective is
  # usually infinite. Resolved through the bound transform it is one.
  r <- minimize(bfgs(), quad, start_zeros(3), gr = quad_g, lower = 0)
  expect_equal(r@par, c(1, 2, 3), tolerance = 1e-6)

  # and a doubly bounded parameter starts in the middle of its interval. The
  # tolerance is on the parameter, so it follows the stopping rule: a gradient
  # of 1e-6 on a unit curvature locates the minimiser to about that, and this
  # run differentiates numerically as well.
  f <- function(p) (p - 0.3)^2
  r2 <- minimize(bfgs(), f, start_zeros(1), lower = 0, upper = 1)
  expect_equal(r2@par, 0.3, tolerance = 1e-5)

  # the value actually handed to the objective, checked directly
  seen <- NULL
  probe <- function(p) { if (is.null(seen)) seen <<- p; sum((p - 1)^2) }
  minimize(gd(maxit = 1), probe, start_zeros(2), lower = c(0, 0),
           upper = c(Inf, 1))
  expect_equal(seen, c(1, 0.5))
})


test_that("a uniform start never lands outside its bounds", {
  set.seed(7)
  for (i in 1:20) {
    r <- minimize(bfgs(maxit = 1), function(p) sum(p^2),
                  start_runif(-8, 8, npar = 3),
                  lower = c(0, 0, -1), upper = c(Inf, 1, 1))
    expect_true(all(r@par > c(0, 0, -1)), label = paste("draw", i))
    expect_true(all(r@par < c(Inf, 1, 1)), label = paste("draw", i))
  }
})


test_that("the number of parameters comes from the bounds when it can", {
  # Bounds are one per parameter, so a vector of them answers the question and
  # nothing has to be probed.
  r <- minimize(bfgs(), quad, start_zeros(), gr = quad_g,
                lower = c(-10, -10, -10), upper = c(10, 10, 10))
  expect_equal(r@par, c(1, 2, 3), tolerance = 1e-6)

  # and lower alone is enough
  r2 <- minimize(bfgs(), quad, start_zeros(), gr = quad_g,
                 lower = c(0, 0, 0))
  expect_equal(r2@par, c(1, 2, 3), tolerance = 1e-6)

  expect_error(
    minimize(bfgs(), quad, start_zeros(), lower = c(0, 0), upper = c(1, 1, 1)),
    "disagree")
})


test_that("a gradient of fixed length pins the count exactly", {
  # A hand-written gradient that spells its components out returns two numbers
  # whatever it is handed, so requiring length(gr(x)) == length(x) leaves one
  # length standing where the objective alone leaves infinitely many.
  f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
  gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                      200 * (p[2] - p[1]^2))
  expect_error(infer_npar(f, NULL, function(k) numeric(k)), "more than one")
  expect_equal(infer_npar(f, gr, function(k) numeric(k)), 2L)

  r <- minimize(bfgs(), f, start_zeros(), gr = gr)
  expect_equal(r@par, c(1, 1), tolerance = 1e-6)
})


test_that("a model with a design matrix decides it outright, which is the real case", {
  # The objective a modelling package actually hands over is not a toy: it
  # multiplies the parameter by something of a fixed width, and X %*% beta with
  # the wrong length is an error rather than a number. This is the case the
  # inference exists for, and it needs no gradient and no help.
  set.seed(2)
  X <- cbind(1, matrix(rnorm(40), 20))
  y <- drop(X %*% c(0.5, -1, 2)) + rnorm(20, sd = 0.2)
  rss <- function(b) sum((y - X %*% b)^2)

  expect_equal(infer_npar(rss, NULL, function(k) numeric(k)), 3L)
  r <- minimize(bfgs(), rss, start_zeros())
  expect_equal(r@par, drop(qr.solve(X, y)), tolerance = 1e-5)
})


test_that("a vectorised objective cannot be decided, and neither can its gradient", {
  # Worth pinning, because both plausible guesses about R are wrong. Recycling
  # warns only when the shorter length is not a DIVISOR of the longer, so
  # sum((p - c(1, 2, 3))^2) accepts a length-one vector in silence and returns a
  # perfectly finite 14; and its gradient 2 * (p - c(1, 2, 3)) returns length
  # six for a length-six argument, so the gradient rule passes it too.
  #
  # Neither is a defect and neither can be guessed at. The refusal names the two
  # lengths it found, which is the information the caller needs.
  expect_warning(quad(numeric(4)), "multiple")     # 4 is not a multiple of 3
  expect_silent(quad(numeric(1)))                  # but 1 divides 3
  expect_length(quad_g(numeric(6)), 6)             # and so does 3 divide 6

  expect_error(infer_npar(quad, NULL, function(k) numeric(k)),
               "more than one length")
  expect_error(infer_npar(quad, quad_g, function(k) numeric(k)),
               "more than one length")

  # Saying so is one word, and the three routes that avoid the probe entirely
  # are npar, vector bounds, and an objective that rejects the wrong length.
  expect_equal(minimize(bfgs(), quad, start_zeros(3), gr = quad_g)@par,
               c(1, 2, 3), tolerance = 1e-6)
})


test_that("an objective that really takes any length is refused, by name", {
  # sum((p - 1)^2) is a perfectly good function of any number of parameters.
  # Guessing would mean optimising a different problem from the one asked.
  expect_error(infer_npar(function(p) sum((p - 1)^2), NULL,
                          function(k) numeric(k)),
               "more than one length")
  expect_error(minimize(bfgs(), function(p) sum((p - 1)^2), start_zeros()),
               "npar")
  # and a gradient does not save this one, since it too takes any length
  expect_error(minimize(bfgs(), function(p) sum((p - 1)^2), start_zeros(),
                        gr = function(p) 2 * (p - 1)),
               "npar")

  # and saying so fixes it
  r <- minimize(bfgs(), function(p) sum((p - 1)^2), start_zeros(4))
  expect_equal(r@par, rep(1, 4), tolerance = 1e-6)
})


test_that("an objective no length satisfies is refused too", {
  expect_error(infer_npar(function(p) stop("no"), NULL, function(k) numeric(k),
                          npar_max = 4),
               "no length")
})


test_that("the probe stops as soon as the answer is known to be ambiguous", {
  # Two accepted lengths are enough to know, so the cost of the bad case is two
  # evaluations rather than npar_max.
  n <- 0
  f <- function(p) { n <<- n + 1; sum(p^2) }
  expect_error(infer_npar(f, NULL, function(k) numeric(k), npar_max = 50),
               "more than one length")
  expect_equal(n, 2)
})


test_that("starters work with every optimiser, including multistart", {
  set.seed(4)
  for (o in list(gd(maxit = 3), cg(maxit = 3), bb(maxit = 3), bfgs(maxit = 3),
                 lbfgs(maxit = 3), newton(maxit = 3), adam(maxit = 3),
                 nelder_mead(maxit = 3), compass(maxit = 3),
                 multistart(bfgs(maxit = 3), n = 2, ncores = 1))) {
    r <- minimize(o, quad, start_zeros(3), gr = quad_g)
    expect_length(r@par, 3)
    expect_true(is.finite(r@value))
  }
  r <- minimize(bfgs(), quad, start_runif(npar = 3), gr = quad_g)
  expect_equal(r@par, c(1, 2, 3), tolerance = 1e-6)
})


test_that("maximize takes a starter too", {
  r <- maximize(bfgs(), function(p) -sum((p - c(1, 2))^2), start_zeros(2))
  expect_equal(r@par, c(1, 2), tolerance = 1e-6)
  expect_equal(r@value, 0, tolerance = 1e-8)
})


test_that("a numeric par is untouched, so the ordinary call pays nothing", {
  a <- minimize(bfgs(), quad, c(0, 0, 0), gr = quad_g)
  b <- minimize(bfgs(), quad, start_zeros(3), gr = quad_g)
  expect_identical(a@par, b@par)
  expect_identical(a@counts, b@counts)
})


test_that("set.seed reproduces a uniform start", {
  set.seed(11); a <- minimize(bfgs(), quad, start_runif(-3, 3, npar = 3))
  set.seed(11); b <- minimize(bfgs(), quad, start_runif(-3, 3, npar = 3))
  expect_identical(a@par, b@par)
})


test_that("the constructors refuse nonsense", {
  expect_error(start_zeros(0), "'npar'")
  expect_error(start_zeros(2.5), "'npar'")
  expect_error(start_zeros("3"), "'npar'")
  expect_error(start_runif(1, 0), "strictly below")
  expect_error(start_runif(-Inf, 1), "finite")
  expect_error(start_runif(NA, 1), "finite")
})
