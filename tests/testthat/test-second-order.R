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
    # the repair is recorded; WHICH of the two strings it is -- plain or
    # "(capped)" -- is the next test's business, and here the claim is only
    # that a repair was reported at all
    expect_true(any(startsWith(r@trace$safeguard, "hessian modified")),
                label = mod)
  }
})


test_that("a repaired Newton step is capped and an unrepaired one is not", {
  # WHY THIS EXISTS. At this start the repaired direction's first component is
  # g_1 / (floor * max|lambda|) = 0.1995 / 3.97e-08 = 5.03e+06: its length is
  # set by `floor` and by no curvature of the objective. Before the cap the
  # line search was the only thing that bounded it, at one function evaluation
  # per backtrack.
  H <- saddle_he(c(0.05, 0.6))
  ev <- eigen(H, symmetric = TRUE)$values
  expect_lt(min(ev), 0)
  lo <- max(1e-8, 1e-8 * max(abs(ev)))
  d1 <- abs(saddle_gr(c(0.05, 0.6))[1]) / lo
  expect_gt(d1, 1e6)

  for (mod in c("eigen", "ridge")) {
    r <- minimize(newton(hessian_mod = mod, keep_trace = TRUE),
                  saddle, c(0.05, 0.6), gr = saddle_gr, he = saddle_he)
    expect_true(r@converged, label = mod)
    expect_equal(r@value, -1, tolerance = 1e-8, label = mod)
    expect_true("hessian modified (capped)" %in% r@trace$safeguard,
                label = mod)
    # measured against optimizers7 0.4.0 on the same problem: 28 function
    # evaluations (eigen) and 30 (ridge) become 5 and 5, the backtracks
    # being what the uncapped length cost. Ten separates them with room.
    expect_lt(r@counts[["f"]], 10, label = mod)
  }

  # IT ONLY EVER SHORTENS, and the two halves of that are asserted apart.
  #
  # (i) a repaired direction already of order one is untouched: at floor = 1
  #     the floored eigenvalue is 3.97 and the direction is (0.05, -0.3), so
  #     the cap cannot fire and the guard stays the plain string. Measured
  #     identical either side of the change: 8 iterations, 9 evaluations.
  r1 <- minimize(newton(floor = 1, keep_trace = TRUE), saddle, c(0.05, 0.6),
                 gr = saddle_gr, he = saddle_he)
  expect_true(r1@converged)
  expect_true("hessian modified" %in% r1@trace$safeguard)
  expect_false("hessian modified (capped)" %in% r1@trace$safeguard)
  expect_identical(r1@iterations, 8L)

  # (ii) A GENUINE NEWTON STEP IS NOT TOUCHED AT ALL. Where the Cholesky
  #     succeeds the step's length is the objective's own curvature and
  #     capping it would throw away the quadratic convergence. Rosenbrock from
  #     this start repairs nothing, and its counts are what they were before
  #     the cap existed: 21 iterations, 29 f, 22 g, 21 h.
  r2 <- minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1), gr = rosen_gr,
                 he = rosen_he)
  expect_false(any(startsWith(r2@trace$safeguard, "hessian modified")))
  expect_identical(r2@iterations, 21L)
  expect_identical(r2@counts[["f"]], 29L)
  expect_identical(r2@counts[["g"]], 22L)
  expect_identical(r2@counts[["h"]], 21L)
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

  # crit_grad() rather than the default: the tolerances below are about how
  # close each METHOD gets, and the 0.6.0 default also stops on a stalled
  # objective. See crit_any().
  cr <- crit_grad()
  small <- minimize(lbfgs(memory = 3, criterion = cr, maxit = 5000),
                    f, rep(0, 30), gr = g)
  large <- minimize(lbfgs(memory = 30, criterion = cr, maxit = 5000),
                    f, rep(0, 30), gr = g)
  full  <- minimize(bfgs(criterion = cr, maxit = 5000), f, rep(0, 30), gr = g)

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


test_that("the quasi-Newton settings are validated at their own bounds", {
  # curv_tol multiplies |s||y|, so a negative value makes the comparison hold
  # for every pair and turns the skip protection off in silence.
  expect_error(bfgs(curv_tol = -1), "'curv_tol'")
  expect_error(lbfgs(curv_tol = -1e-16), "'curv_tol'")
  expect_error(bfgs(curv_tol = Inf), "'curv_tol'")

  # max_skip is compared with `>=`, so a negative value behaves as zero.
  expect_error(bfgs(max_skip = -3), "'max_skip'")
  expect_error(bfgs(max_skip = 2.5), "whole")

  # memory reaches the kernel through as.integer(), so 2.5 would run at 2.
  expect_error(lbfgs(memory = 2.5), "whole")

  # The negative controls: the values either side of each bound still build,
  # and zero is a setting rather than a mistake in the two places it means
  # something -- curv_tol = 0 is the textbook condition s'y > 0, and
  # max_skip = 0 resets on the first skip.
  expect_s7_class(bfgs(curv_tol = 0), optimizers7:::Bfgs)
  expect_s7_class(bfgs(max_skip = 0), optimizers7:::Bfgs)
  expect_s7_class(bfgs(curv_tol = 1e10), optimizers7:::Bfgs)
  expect_s7_class(lbfgs(memory = 1), optimizers7:::Lbfgs)
  expect_s7_class(lbfgs(curv_tol = 0), optimizers7:::Lbfgs)
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


test_that("the first quasi-Newton step is of order one in the parameters", {
  # A quasi-Newton direction is scaled so that step = 1 means the NEWTON step,
  # which is why bfgs() and lbfgs() default to it. On the first iteration
  # there is no curvature and the direction degenerates to -g, for which one
  # is not a natural unit: the trial displacement IS the gradient, so on a
  # badly scaled objective the first point tried is an arbitrary distance
  # away and the line search pays to backtrack all the way in.
  #
  # A quadratic 0.5*c*|x|^2 makes the arithmetic explicit. From x0 = 1 the
  # gradient is c and the minimizer along -g sits at alpha = 1/c, so an
  # unscaled first direction needs about log2(c) halvings to reach it, while
  # a direction of unit max-norm reaches it at alpha = 1 and the line search
  # accepts immediately.
  cs <- 1e4
  f <- function(p) 0.5 * cs * sum(p^2)
  g <- function(p) cs * p
  x0 <- c(1, 1)

  for (o in list(bfgs(), lbfgs())) {
    r <- minimize(o, f, x0, gr = g)
    lab <- o@name
    expect_true(r@converged, info = lab)
    # Measured with the scaling and without it, on this problem: one
    # iteration and 4 evaluations against two iterations and 19, for both
    # methods. Ten separates them with room on either side.
    expect_lt(sum(as.numeric(r@counts)), 10, label = lab)
  }

  # and it only ever SHORTENS: a problem whose gradient at the start is
  # already of order one takes the identical path, which is what makes this
  # safe to change under two packages that depend on the default. Here
  # max|g| = 0.6 at the start, so the scaling does not fire -- measured, 2
  # iterations and 6 evaluations either way.
  f2 <- function(p) 0.15 * sum(p^2)
  g2 <- function(p) 0.3 * p
  r2 <- minimize(lbfgs(), f2, c(2, -1), gr = g2)
  expect_true(r2@converged)
  expect_lt(max(abs(r2@par)), 1e-6)
  expect_identical(r2@iterations, 2L)
})


test_that("a Newton step longer than max_length is shortened to it", {
  # sqrt(1 + x^2) is convex and its curvature falls as |x| grows, so the
  # Newton step from x = 2 is -x (1 + x^2) = -10 and overshoots the minimum
  # at zero by eight. That is the shape of a criterion flattening towards an
  # asymptote, where a long step lands in the flat region.
  f  <- function(p) sqrt(1 + p^2)
  gr <- function(p) p / sqrt(1 + p^2)
  he <- function(p) matrix((1 + p^2)^(-3 / 2), 1, 1)
  expect_equal(-gr(2) / he(2)[1, 1], -10)

  free <- minimize(newton(keep_trace = TRUE), f, 2, gr = gr, he = he)
  capped <- minimize(newton(max_length = 1, keep_trace = TRUE), f, 2,
                     gr = gr, he = he)
  for (r in list(free, capped)) {
    expect_true(r@converged)
    expect_equal(r@par, 0, tolerance = 1e-6)
  }
  # the bounded run takes the steps -1 and -1 at full length and says so,
  # where the free one has to backtrack its first step to a quarter
  expect_true(all(capped@trace$safeguard == "newton step capped"))
  expect_true(all(capped@trace$step == 1))
  expect_false("newton step capped" %in% free@trace$safeguard)
  expect_lt(capped@counts[["f"]], free@counts[["f"]])

  # A BOUND THAT DOES NOT BIND CHANGES NOTHING. Rosenbrock's Newton steps from
  # the customary start are all shorter than 10, and the run is the one
  # newton() gives without the argument, to the iteration and the evaluation.
  a <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
  b <- minimize(newton(max_length = 10), rosen, c(-1.2, 1), gr = rosen_gr,
                he = rosen_he)
  expect_identical(b@par, a@par)
  expect_identical(b@counts, a@counts)
})


test_that("max_length is validated", {
  expect_error(newton(max_length = 0), "positive")
  expect_error(newton(max_length = -1), "positive")
  expect_error(newton(max_length = NA_real_), "positive")
  expect_error(newton(max_length = c(1, 2)), "positive")
  expect_identical(newton()@max_length, Inf)
})


test_that("typical sizes read every length in the parameters' own units", {
  expect_null(newton()@typical)
  expect_error(newton(typical = 0), "positive")
  expect_error(newton(typical = c(1, NA)), "positive")
  expect_error(newton(typical = Inf), "positive")
  expect_error(minimize(newton(typical = c(1, 2, 3)), rosen, c(-1.2, 1),
                        gr = rosen_gr, he = rosen_he), "one per parameter")
  # two separate quadratics whose minima lie 1000 and 2 away: the Newton step
  # is (1000, -2). A bound of 5 in the parameters' own units shrinks the step
  # by 200; with typical sizes (1000, 1) the same bound leaves it whole, and
  # with (100, 1) it is halved, the first component reaching 5 typical sizes.
  f  <- function(p) (p[1] - 1000)^2 + (p[2] + 2)^2
  gr <- function(p) 2 * c(p[1] - 1000, p[2] + 2)
  he <- function(p) diag(2, 2)
  one  <- minimize(newton(max_length = 5, keep_trace = TRUE, maxit = 1),
                   f, c(0, 0), gr = gr, he = he)
  wide <- minimize(newton(max_length = 5, typical = c(1000, 1),
                          keep_trace = TRUE, maxit = 1),
                   f, c(0, 0), gr = gr, he = he)
  half <- minimize(newton(max_length = 5, typical = c(100, 1),
                          keep_trace = TRUE, maxit = 1),
                   f, c(0, 0), gr = gr, he = he)
  expect_equal(one@par, c(5, -0.01))
  expect_equal(wide@par, c(1000, -2))
  expect_equal(half@par, c(500, -1))
  expect_false("newton step capped" %in% wide@trace$safeguard)
  expect_true("newton step capped" %in% half@trace$safeguard)

  # A REPAIRED direction is scaled to one typical size, not to one unit.
  # -100 cos(p / 100) is concave at p = 250, so the eigen floor fires and the
  # direction is scaled to length one: one unit, or one typical size.
  g2  <- function(p) -cos(p[1] / 100) * 100
  gr2 <- function(p) sin(p[1] / 100)
  he2 <- function(p) matrix(cos(p[1] / 100) / 100, 1, 1)
  # at p = 250 the curvature cos(2.5)/100 is negative and the gradient
  # sin(2.5) positive, so the step is a repaired one pointing down
  a <- minimize(newton(maxit = 1, keep_trace = TRUE), g2, 250, gr = gr2,
                he = he2)
  b <- minimize(newton(maxit = 1, keep_trace = TRUE, typical = 100), g2, 250,
                gr = gr2, he = he2)
  expect_equal(250 - a@par, 1)
  expect_equal(250 - b@par, 100)

  # AND NOTHING MOVES WITHOUT THEM, to the iteration and the evaluation
  r1 <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
  r2 <- minimize(newton(typical = 1), rosen, c(-1.2, 1), gr = rosen_gr,
                 he = rosen_he)
  expect_identical(r2@par, r1@par)
  expect_identical(r2@counts, r1@counts)
})
