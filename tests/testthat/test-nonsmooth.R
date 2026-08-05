# The three methods that do not need a smooth objective.
#
# The two tests that matter most are the ones where an unguarded implementation
# gives a plausible wrong answer: McKinnon's function, on which Nelder-Mead
# converges to a point that is not a minimizer while every value-based rule
# reports success, and a kinked objective, on which a descent method arrives at
# the answer and then reports that it did not.

quad <- function(p) sum((p - c(1, 2))^2)


test_that("all three solve a smooth problem, slowly", {
  set.seed(1)
  for (o in list(nelder_mead(), compass(), compass(directions = "mads"))) {
    r <- minimize(o, quad, c(0, 0))
    expect_true(r@converged, label = r@optimizer@name)
    expect_equal(r@par, c(1, 2), tolerance = 1e-6, label = r@optimizer@name)
  }
  r <- minimize(bundle(), quad, c(0, 0), gr = function(p) 2 * (p - c(1, 2)))
  expect_equal(r@par, c(1, 2), tolerance = 1e-6)
})


# --- McKinnon ---------------------------------------------------------------

# f(x, y) = theta phi |x|^tau + y + y^2 for x <= 0, theta x^tau + y + y^2
# otherwise. Strictly convex with continuous first derivatives at tau = 2, with
# its minimum at (0, -1/2) and value -1/4. From McKinnon's initial simplex,
# Nelder-Mead performs inside contractions for ever and converges to (0, 0),
# where the value is 0 and the gradient is not zero.
mckinnon <- function(p, tau = 2, theta = 6, phi = 60) {
  x <- p[1]; y <- p[2]
  (if (x <= 0) theta * phi * abs(x)^tau else theta * x^tau) + y + y^2
}
mck_simplex <- rbind(c(0, 0), c((1 + sqrt(33)) / 8, (1 - sqrt(33)) / 8), c(1, 1))


test_that("Nelder-Mead converges to a non-minimizer without the safeguard", {
  # The failure must be reproduced before the fix for it can be believed.
  r <- minimize(nelder_mead(simplex = mck_simplex, max_restarts = 0),
                mckinnon, par = c(0, 0))
  expect_equal(r@par, c(0, 0), tolerance = 1e-6)
  expect_equal(r@value, 0, tolerance = 1e-9)
  # and it says it converged, which is the point: nothing in the values is wrong
  expect_true(r@converged)
})


test_that("the degeneracy restart rescues it", {
  r <- minimize(nelder_mead(simplex = mck_simplex, max_restarts = 3,
                            keep_trace = TRUE),
                mckinnon, par = c(0, 0))
  expect_equal(r@par, c(0, -0.5), tolerance = 1e-6)
  expect_equal(r@value, -0.25, tolerance = 1e-9)
  expect_match(r@message, "simplex restarted")
  expect_true(any(r@trace$safeguard == "restart"))
})


test_that("the restart does not disturb a run that never degenerates", {
  a <- minimize(nelder_mead(max_restarts = 0), quad, c(0, 0))
  b <- minimize(nelder_mead(max_restarts = 3), quad, c(0, 0))
  expect_equal(a@par, b@par)
  expect_equal(a@counts[["f"]], b@counts[["f"]])
  expect_false(grepl("restarted", b@message))
})


# --- pattern search ---------------------------------------------------------

# |p1 + p2| + 0.1 ||p||^2: the kink runs along the anti-diagonal, so no
# coordinate direction descends it. The minimum is at the origin, value 0.
kink <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)


test_that("a coordinate poll stalls on a diagonal kink and a mads one does not", {
  # compass() defaults to mads, so the coordinate poll has to be asked for.
  co <- minimize(compass(directions = "coordinate"), kink, c(1, 0.5))
  set.seed(2)
  ra <- minimize(compass(directions = "mads"), kink, c(1, 0.5))

  # the coordinate poll stops on the ridge, at a point that is not stationary
  expect_gt(co@value, 1e-2)
  expect_equal(co@par[1] + co@par[2], 0, tolerance = 1e-6)
  # mads directions do better, by an order of magnitude
  expect_lt(ra@value, co@value / 5)
  # ...but neither reaches the minimum. That is not a defect in either: a poll
  # of finitely many directions is the wrong instrument for a kink, and it is
  # why bundle() exists. It solves this to 1e-7 in the test below.
  expect_gt(ra@value, 1e-9)
})


test_that("a random poll is reproducible under set.seed", {
  set.seed(4); a <- minimize(compass(directions = "mads"), quad, c(0, 0))
  set.seed(4); b <- minimize(compass(directions = "mads"), quad, c(0, 0))
  expect_identical(a@par, b@par)
})


test_that("complete polling needs fewer iterations than opportunistic", {
  a <- minimize(compass(opportunistic = TRUE), quad, c(0, 0))
  b <- minimize(compass(opportunistic = FALSE), quad, c(0, 0))
  expect_lte(b@iterations, a@iterations)
})


# --- the bundle method ------------------------------------------------------

test_that("bundle finds the median, where a descent method cannot certify it", {
  set.seed(1)
  y <- rnorm(101)                       # odd, so the minimizer is unique
  f <- function(p) sum(abs(y - p))
  g <- function(p) -sum(sign(y - p))

  r <- minimize(bundle(), f, par = 0, gr = g)
  expect_true(r@converged)
  expect_equal(r@par, median(y), tolerance = 1e-10)

  # The motivation, asserted rather than described: bfgs walks to the same
  # point and then reports failure, because the subgradient it evaluates there
  # has norm 1 and no line search can improve on it.
  b <- minimize(bfgs(), f, par = 0, gr = g)
  expect_false(b@converged)
  expect_equal(b@value, r@value, tolerance = 1e-10)
})


test_that("with an even sample the minimizer is a segment, and it lands in it", {
  # sum|y - mu| is piecewise linear with slope (#below - #above), which is
  # exactly zero between the two middle order statistics. Every point of that
  # segment is a minimizer; median() returns its midpoint by convention. Testing
  # against median() with a tolerance would therefore be testing a convention,
  # and would report an error where there is none.
  set.seed(5)
  y <- rnorm(200)
  s <- sort(y)
  f <- function(p) sum(abs(y - p))
  g <- function(p) -sum(sign(y - p))

  r <- minimize(bundle(), f, par = 0, gr = g)
  expect_gte(r@par, s[100])
  expect_lte(r@par, s[101])
  # the objective is flat across it, to the last digit
  expect_equal(f(s[100]), f(s[101]))
  expect_equal(r@value, f(s[100]), tolerance = 1e-9)
})


test_that("bundle solves the diagonal kink that defeats a poll", {
  g <- function(p) rep(sign(p[1] + p[2]), 2) + 0.2 * p
  r <- minimize(bundle(), kink, c(1, 0.5), gr = g)
  expect_true(r@converged)
  expect_lt(r@value, 1e-6)
  expect_lt(max(abs(r@par)), 1e-3)
})


test_that("serious and null steps are both counted and reported", {
  set.seed(1)
  y <- rnorm(51)
  r <- minimize(bundle(keep_trace = TRUE), function(p) sum(abs(y - p)),
                par = 0, gr = function(p) -sum(sign(y - p)))
  expect_match(r@message, "serious")
  expect_match(r@message, "null")
  expect_true(any(r@trace$safeguard == "null"))
  # the predicted decrease is what the rule watched, and it fell
  expect_true(all(r@trace$stationarity >= 0))
  expect_lt(tail(r@trace$stationarity, 1), r@trace$stationarity[1])
})


test_that("the reported gradient is the aggregate, and it is small", {
  set.seed(1)
  y <- rnorm(101)
  r <- minimize(bundle(), function(p) sum(abs(y - p)), par = 0,
                gr = function(p) -sum(sign(y - p)))
  # every individual subgradient there has |g| >= 1; the aggregate does not
  expect_lt(abs(r@gradient), 1e-3)
})


test_that("bundle keeps its bundle bounded and still converges", {
  set.seed(3)
  y <- rnorm(101)
  small <- minimize(bundle(bundle_size = 3), function(p) sum(abs(y - p)),
                    par = 0, gr = function(p) -sum(sign(y - p)))
  expect_true(small@converged)
  expect_equal(small@par, median(y), tolerance = 1e-8)
})


# --- what they refuse, and what they claim ----------------------------------

test_that("a gradient rule is refused by all three", {
  for (o in list(nelder_mead(criterion = crit_grad()),
                 compass(criterion = crit_grad()),
                 bundle(criterion = crit_grad()))) {
    expect_error(minimize(o, quad, c(0, 0)), "needs gradient")
  }
})


test_that("an objective rule is accepted by all three", {
  expect_silent(minimize(nelder_mead(criterion = crit_abs_obj(), maxit = 20),
                         quad, c(0, 0)))
  expect_silent(minimize(compass(criterion = crit_abs_obj(), maxit = 20),
                         quad, c(0, 0)))
  expect_silent(minimize(bundle(criterion = crit_abs_obj(), maxit = 20),
                         quad, c(0, 0), gr = function(p) 2 * (p - c(1, 2))))
})


test_that("a derivative-free method does not claim to have differenced one", {
  # It computes no gradient at all, so reporting that one was obtained by
  # finite differences would be a false statement about the run's exactness.
  r <- minimize(nelder_mead(), quad, c(0, 0))
  expect_false(grepl("finite difference", r@message))
  expect_null(r@gradient)

  r <- minimize(compass(), quad, c(0, 0))
  expect_false(grepl("finite difference", r@message))
  expect_null(r@gradient)

  # bundle does use one, so with no gr supplied it must say so
  r <- minimize(bundle(maxit = 20), quad, c(0, 0))
  expect_match(r@message, "finite difference")
})


test_that("bounds work for all three, on the user's scale", {
  f <- function(p) sum((p - c(1, 2))^2)
  for (o in list(nelder_mead(), compass(), bundle())) {
    r <- minimize(o, f, c(0.5, 0.5), gr = function(p) 2 * (p - c(1, 2)),
                  lower = c(0, 0), upper = c(5, 1))
    expect_gt(r@par[2], 0.9); expect_lt(r@par[2], 1)
    expect_equal(r@par[1], 1, tolerance = 1e-3, label = r@optimizer@name)
  }
})


test_that("the stationarity column appears only for these methods", {
  r <- minimize(nelder_mead(maxit = 5, keep_trace = TRUE), quad, c(0, 0))
  expect_true("stationarity" %in% names(r@trace))
  expect_false("gnorm" %in% names(r@trace))

  r <- minimize(bfgs(maxit = 5, keep_trace = TRUE), quad, c(0, 0))
  expect_true("gnorm" %in% names(r@trace))
  expect_false("stationarity" %in% names(r@trace))
})


test_that("crit_stationary reads what it says", {
  expect_identical(crit_needs(crit_stationary()), "stationarity")
  expect_true(crit_met(crit_stationary(1e-6), list(stationarity = 1e-9)))
  expect_false(crit_met(crit_stationary(1e-6), list(stationarity = 1e-3)))
  expect_false(crit_met(crit_stationary(1e-6), list(stationarity = NULL)))
})


test_that("the constructors refuse nonsense", {
  expect_error(nelder_mead(step = 0), "'step'")
  expect_error(nelder_mead(adaptive = 1), "'adaptive'")
  expect_error(nelder_mead(max_restarts = -1), "'max_restarts'")
  expect_error(nelder_mead(simplex = 1:3), "'simplex'")
  expect_error(minimize(nelder_mead(simplex = matrix(0, 2, 2)), quad, c(0, 0)),
               "3 rows")

  expect_error(compass(directions = "diagonal"), "should be one of")
  expect_error(compass(expand = 0.5), "'expand'")
  expect_error(compass(shrink = 1), "'shrink'")
  expect_error(compass(opportunistic = NA), "'opportunistic'")

  expect_error(bundle(m_serious = 0), "'m_serious'")
  expect_error(bundle(m_serious = 1), "'m_serious'")
  expect_error(bundle(bundle_size = 1), "'bundle_size'")
  expect_error(bundle(t_min = 10, t_max = 1), "strictly below")
  expect_error(bundle(t0 = 1e20), "between")
})
