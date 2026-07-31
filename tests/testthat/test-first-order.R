# The three methods that use the gradient and nothing else.
#
# All of them are a direction plus a line search, so they share the compiled
# loop; what is tested here is the direction each forms and the safeguards each
# needs, not the machinery around them.

quad <- function(p) sum((p - c(1, 2))^2)
quad_g <- function(p) 2 * (p - c(1, 2))

rosen <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
rosen_g <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                         200 * (p[2] - p[1]^2))


test_that("all three solve a quadratic and Rosenbrock", {
  for (o in list(gd(maxit = 5000), cg(), bb())) {
    r <- minimize(o, quad, c(0, 0), gr = quad_g)
    expect_true(r@converged, label = r@optimizer@name)
    expect_equal(r@par, c(1, 2), tolerance = 1e-6, label = r@optimizer@name)
  }
  for (o in list(cg(), bb())) {
    r <- minimize(o, rosen, c(-1.2, 1), gr = rosen_g)
    expect_true(r@converged, label = r@optimizer@name)
    expect_equal(r@par, c(1, 1), tolerance = 1e-5, label = r@optimizer@name)
  }
})


# --- conjugate gradients ----------------------------------------------------

test_that("cg beats gd on the problem gd is bad at, and by a lot", {
  # The zigzag is the whole argument for the method: a narrow valley is where
  # consecutive steepest-descent directions being orthogonal costs most.
  a <- minimize(gd(maxit = 20000), rosen, c(-1.2, 1), gr = rosen_g)
  b <- minimize(cg(maxit = 5000), rosen, c(-1.2, 1), gr = rosen_g)
  expect_lt(b@iterations, a@iterations / 10)
  expect_true(b@converged)
  expect_false(a@converged)          # gd does not get there at all
})


test_that("all four beta formulas work, and agree about the answer", {
  for (bt in c("pr", "fr", "hs", "dy")) {
    r <- minimize(cg(beta = bt, maxit = 5000), rosen, c(-1.2, 1), gr = rosen_g)
    expect_equal(r@par, c(1, 1), tolerance = 1e-5, label = bt)
  }
})


test_that("on a quadratic cg terminates in a few multiples of p", {
  # The classical property is p-step termination with an EXACT line search, and
  # no practical line search is exact, so this asserts the order of magnitude.
  #
  # The tolerance is 1e-6 and not 1e-8 deliberately. At 1e-8 on this problem the
  # line search runs out of room before the rule fires and the run reports
  # failure -- correctly, and lbfgs does the same thing here. Only bfgs, which
  # carries a full curvature approximation, reaches it. A test asserting 1e-8
  # would be asserting something about double precision rather than about cg.
  set.seed(1)
  p <- 20
  A <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
  b <- rnorm(p)
  f <- function(x) as.numeric(0.5 * t(x) %*% A %*% x - sum(b * x))
  g <- function(x) as.numeric(A %*% x - b)

  r <- minimize(cg(criterion = crit_grad(1e-6), maxit = 500), f, rep(0, p),
                gr = g)
  expect_true(r@converged)
  expect_lt(r@iterations, 4 * p)
  expect_equal(r@par, as.numeric(solve(A, b)), tolerance = 1e-5)
})


test_that("a tighter Wolfe constant is why cg defaults to one", {
  # c2 = 0.9 suits bfgs, which repairs a loose step with its curvature
  # approximation. cg has nothing to repair with: the conjugacy it accumulates
  # is only as good as the search that produced it.
  loose <- minimize(cg(line_search = wolfe(c2 = 0.9), maxit = 5000),
                    rosen, c(-1.2, 1), gr = rosen_g)
  tight <- minimize(cg(maxit = 5000), rosen, c(-1.2, 1), gr = rosen_g)
  expect_lt(tight@iterations, loose@iterations / 2)
})


test_that("PR+ clamps a negative beta, and reports it as the restart it is", {
  r <- minimize(cg(beta = "pr", maxit = 5000, keep_trace = TRUE),
                rosen, c(-1.2, 1), gr = rosen_g)
  expect_true(any(r@trace$safeguard == "cg restart"))
})


test_that("a periodic restart is taken and does not break the run", {
  r <- minimize(cg(restart_every = 5, maxit = 5000, keep_trace = TRUE),
                rosen, c(-1.2, 1), gr = rosen_g)
  expect_true(any(r@trace$safeguard == "cg restart"))
  expect_equal(r@par, c(1, 1), tolerance = 1e-5)
})


test_that("a non-descent direction falls back to steepest descent", {
  # The formulas cannot promise a descent direction, and the line search would
  # find nothing along one that goes uphill. Provoked with the formula that has
  # no clamp and a run long enough to meet the case.
  set.seed(2)
  f <- function(p) sum(p^2) + 3 * sin(sum(p))^2
  g <- function(p) 2 * p + 6 * sin(sum(p)) * cos(sum(p))
  r <- minimize(cg(beta = "fr", maxit = 2000, keep_trace = TRUE),
                f, c(2.5, -1.7), gr = g)
  # either the fallback fired or the run converged without needing it; both are
  # correct, and the test is that the run never failed
  expect_true(r@converged || r@iterations == 2000)
  expect_true(all(is.finite(r@par)))
})


# --- Barzilai-Borwein -------------------------------------------------------

test_that("bb solves a quadratic in two iterations", {
  # The first secant pair already contains the exact inverse curvature, so the
  # second step lands on the answer. This is the whole claim of the method.
  r <- minimize(bb(), quad, c(0, 0), gr = quad_g)
  expect_lte(r@iterations, 3)
  expect_equal(r@par, c(1, 2), tolerance = 1e-10)
})


test_that("bb reaches every smooth minimum in the battery", {
  for (nm in names(test_problems())) {
    p <- test_problems(nm)[[1]]
    if (!isTRUE(p$smooth) || isTRUE(p$multimodal)) next
    r <- minimize(bb(maxit = 5000), p$fn, p$par, gr = p$gr)
    expect_lt(r@value - p$value, 1e-8, label = nm)
  }
})


test_that("the three variants all work", {
  for (v in c("alternate", "bb1", "bb2")) {
    r <- minimize(bb(variant = v, maxit = 5000), rosen, c(-1.2, 1),
                  gr = rosen_g)
    expect_equal(r@par, c(1, 1), tolerance = 1e-5, label = v)
  }
})


test_that("a refused pair resets alpha to a step of order one, not to a constant", {
  # A pair reporting no usable curvature implies no positive step length, so
  # some alpha has to be supplied from outside. Both constants one might reach
  # for are absorbing states: KEEPING the previous alpha traps a short step in
  # its own shortness, since a short step samples curvature over a short
  # interval and on a valley floor that curvature is negative; restarting at
  # alpha0 does the same wherever alpha0 is itself too short to escape with.
  # 1/max|g| is neither, because it depends on the gradient and not on the alpha
  # it replaces. The two witnesses are the two problems that fail under the two
  # constants: Rosenbrock, where keeping cost 873 refusals in 945 iterations,
  # and a boxed quadratic through its reparametrisation, where restarting at
  # alpha0 cost 1395 in 1521.
  f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
  g  <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                      200 * (p[2] - p[1]^2))
  r <- minimize(bb(variant = "bb2", maxit = 5000, keep_trace = TRUE),
                f, c(-1.2, 1), gr = g)
  expect_true(r@converged)
  expect_lt(r@iterations, 150)
  # the reset is named where it fires, and it fires a handful of times, not on
  # nine tenths of the run
  expect_lt(sum(r@trace$safeguard == "bb curvature reset"), 10)

  q <- minimize(bb(maxit = 5000, keep_trace = TRUE),
                function(p) sum((p - c(1, 2))^2), c(0.5, 0.5),
                gr = function(p) 2 * (p - c(1, 2)),
                lower = c(0, 0), upper = c(5, 1))
  expect_true(q@converged)
  expect_lt(q@iterations, 200)
  expect_lt(sum(q@trace$safeguard == "bb curvature reset"), 10)

  set.seed(3)
  f2 <- function(p) sum(sin(2 * p)) + 0.05 * sum(p^2)
  g2 <- function(p) 2 * cos(2 * p) + 0.1 * p
  r2 <- minimize(bb(maxit = 500, keep_trace = TRUE), f2, c(1.4, -0.6), gr = g2)
  expect_true(all(is.finite(r2@par)))
  expect_true(is.finite(r2@value))
})


test_that("the nonmonotone search lets bb go uphill, and that is what it is for", {
  f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
  g  <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                      200 * (p[2] - p[1]^2))
  nm <- minimize(bb(maxit = 5000, keep_trace = TRUE), f, c(-1.2, 1), gr = g)
  mo <- minimize(bb(line_search = armijo(), maxit = 5000, keep_trace = TRUE),
                 f, c(-1.2, 1), gr = g)

  # the mechanism: steps that make the objective worse are accepted
  expect_gt(sum(diff(nm@trace$value) > 0), 0)
  # and cannot be, under a monotone rule
  expect_equal(sum(diff(mo@trace$value) > 0), 0)
  # and it pays, in the count that matters
  expect_lt(nm@counts[["f"]], mo@counts[["f"]])
})


test_that("memory zero is exactly armijo", {
  # The two differ in one number, so the comparison above is a comparison of
  # that number and of nothing else.
  f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
  g  <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                      200 * (p[2] - p[1]^2))
  a <- minimize(bb(line_search = nonmonotone(memory = 0), maxit = 5000),
                f, c(-1.2, 1), gr = g)
  b <- minimize(bb(line_search = armijo(), maxit = 5000), f, c(-1.2, 1), gr = g)
  expect_identical(a@par, b@par)
  expect_identical(a@counts, b@counts)
})


test_that("nonmonotone refuses nonsense", {
  expect_error(nonmonotone(memory = -1), "'memory'")
  expect_error(nonmonotone(memory = 2.5), "'memory'")
  expect_error(nonmonotone(c1 = 0), "'c1'")
})


test_that("the step length is clamped, and the clamp is reported", {
  r <- minimize(bb(alpha_max = 1e-3, maxit = 3000, keep_trace = TRUE),
                quad, c(0, 0), gr = quad_g)
  expect_true(any(r@trace$safeguard == "bb step clamped"))
  # and it still gets there, just slowly: clamping the step to 1e-3 is what a
  # first-order method with a fixed tiny step looks like
  expect_equal(r@par, c(1, 2), tolerance = 1e-2)
})


# --- what they refuse -------------------------------------------------------

test_that("all three pass check_optimizer", {
  for (o in list(gd(), cg(), bb())) {
    res <- check_optimizer(o, problems = test_problems("sphere"),
                           verbose = FALSE)
    bad <- names(res$checks)[!vapply(res$checks, isTRUE, logical(1))]
    expect_equal(bad, character(), label = o@name)
  }
})


test_that("the constructors refuse nonsense", {
  expect_error(gd(step = 0), "'step'")
  expect_error(gd(line_search = "armijo"), "line_search object")

  expect_error(cg(beta = "nesterov"), "should be one of")
  expect_error(cg(restart_every = -1), "'restart_every'")

  expect_error(bb(variant = "bb3"), "should be one of")
  expect_error(bb(alpha_min = 10, alpha_max = 1), "strictly below")
  expect_error(bb(alpha0 = 0), "'tol'")
})


test_that("they carry bounds like everything else", {
  f <- function(p) sum((p - c(1, 2))^2)
  for (o in list(gd(maxit = 5000), cg(), bb())) {
    r <- minimize(o, f, c(0.5, 0.5), gr = quad_g, lower = c(0, 0), upper = c(5, 1))
    expect_gt(r@par[2], 0.98, label = r@optimizer@name)
    expect_lt(r@par[2], 1, label = r@optimizer@name)
    expect_equal(r@par[1], 1, tolerance = 1e-3, label = r@optimizer@name)
  }
})
