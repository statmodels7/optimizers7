# The skeleton's job is the frame, not the algorithm, so this is what the frame
# has to get right: the three shapes an objective arrives in, the stopping rule
# as an object a user can replace, the budgets, the safeguards, and an honest
# account of what happened.

quad <- function(target) function(p) sum((p - target)^2)
quad_gr <- function(target) function(p) 2 * (p - target)

# Rosenbrock: gradient descent does not solve it quickly, which is exactly what
# a test of the iteration budget needs.
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) {
  c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
    200 * (p[2] - p[1]^2))
}


test_that("the three shapes of objective reach the same answer", {
  target <- c(1, 2)
  f <- quad(target)

  with_gr <- minimize(gd(), f, c(0, 0), gr = quad_gr(target))
  no_gr   <- minimize(gd(), f, c(0, 0))

  expect_equal(with_gr@par, target, tolerance = 1e-6)
  expect_equal(no_gr@par, target, tolerance = 1e-6)

  # the run says so when the gradient was differenced rather than given
  expect_match(no_gr@message, "finite differences")
  expect_false(grepl("finite differences", with_gr@message))
  expect_equal(with_gr@counts[["g"]] > 0, TRUE)
  expect_equal(no_gr@counts[["g"]], 0)
})


test_that("a resampling objective is just an objective", {
  # finite_sum() is gone. It declared that an objective was a sum over
  # observations so that a stochastic method could ask it for a subset, and it
  # had exactly one caller. An objective that resamples is a closure and needs
  # no class, which is what this asserts.
  set.seed(42)
  y <- rnorm(200, mean = 3)
  batch <- function(par) { i <- sample.int(200, 40); sum((y[i] - par)^2) / 2 }
  batch_gr <- function(par) { i <- sample.int(200, 40); -sum(y[i] - par) }

  res <- minimize(adam(alpha = 0.05, decay = 0.01, maxit = 2000),
                  batch, par = 0, gr = batch_gr)
  expect_equal(res@par, mean(y), tolerance = 5e-2)
})


test_that("criteria are objects, compose, and can be written by a user", {
  expect_s3_class(crit_grad(), "S7_object")
  expect_equal(crit_needs(crit_grad()), "gradient")
  # Every optimiser evaluates the objective, so a rule that reads it needs
  # nothing declared and can never be refused. There was a token for it while
  # Adam drew its own minibatches and the objective could be an estimate; with
  # that gone the token could never be unsatisfied, so it went too.
  expect_length(crit_needs(crit_rel_obj()), 0)
  expect_length(crit_needs(crit_abs_par()), 0)

  both <- crit_any(crit_grad(1e-8), crit_rel_obj(1e-12))
  expect_equal(crit_needs(both), "gradient")

  # x has genuinely moved here, so the parameter rule is not satisfied and the
  # conjunction must fail even though the gradient rule is satisfied.
  st <- list(iter = 1, f_new = 1, f_old = 1 + 1e-20,
             x_new = c(1, 1), x_old = c(1, 1.5), gradient = c(1e-12, 1e-12))
  expect_true(crit_met(crit_grad(1e-8), st))
  expect_true(crit_met(both, st))
  expect_false(crit_met(crit_abs_par(1e-8), st))
  expect_false(crit_met(crit_all(crit_grad(1e-8), crit_abs_par(1e-8)), st))
  expect_true(crit_met(crit_any(crit_grad(1e-8), crit_abs_par(1e-8)), st))

  # a criterion nobody in this package wrote
  CritTiny <- S7::new_class("CritTiny", parent = criterion,
                            properties = list(tol = S7::class_numeric))
  S7::method(crit_met, CritTiny) <- function(criterion, state) {
    state$f_new < criterion@tol
  }
  mine <- CritTiny(label = "objective below 1e-6", tol = 1e-6)
  res <- minimize(gd(criterion = mine), quad(c(1, 2)), c(0, 0),
                  gr = quad_gr(c(1, 2)))
  expect_true(res@converged)
  expect_equal(res@criterion_met, "objective below 1e-6")
})


test_that("a criterion the method cannot evaluate is refused, not ignored", {
  # A rule needing a gradient, given to a method that computes none, would test
  # NULL every iteration and never fire; the run would end on the budget and
  # report a reason nowhere near the truth.
  NoGrad <- S7::new_class("NoGrad", parent = optimizer)
  # It computes no gradient and says so. What an optimiser declares is what it
  # can honestly supply, so a rule reading a gradient is refused while one
  # reading the objective is not -- every optimiser evaluates that.
  S7::method(optimizer_provides, NoGrad) <- function(optimizer) character()
  # The signature must include every named argument of the generic, bounds
  # among them, or S7 refuses to register the method.
  S7::method(minimize, NoGrad) <-
    function(optimizer, fn, par, gr = NULL, he = NULL, lower = -Inf, upper = Inf, ...) {
      prepare_objective(optimizer, fn, par, gr)
      "unreachable"
    }
  ng <- NoGrad(name = "no-gradient method", criterion = crit_grad(1e-8),
               maxit = 10, max_eval = 100, verbose = FALSE, refresh = 1,
               keep_trace = FALSE)
  expect_error(minimize(ng, quad(c(0, 0)), c(1, 1)), "needs gradient")

  # and one it can evaluate is accepted
  ng2 <- NoGrad(name = "no-gradient method", criterion = crit_rel_obj(),
                maxit = 10, max_eval = 100, verbose = FALSE, refresh = 1,
                keep_trace = FALSE)
  expect_equal(minimize(ng2, quad(c(0, 0)), c(1, 1)), "unreachable")
})


test_that("converged is never TRUE merely because a budget ran out", {
  # The commonest defect in a hand-written loop, and the one that turns a
  # failure into a wrong answer that looks right.
  res <- minimize(gd(maxit = 5, criterion = crit_grad(1e-12)),
                  rosen, c(-1.2, 1), gr = rosen_gr)
  expect_false(res@converged)
  expect_equal(res@iterations, 5)
  expect_match(res@criterion_met, "iteration budget")

  ev <- minimize(gd(maxit = 1000, max_eval = 20,
                                  criterion = crit_grad(1e-12)),
                 rosen, c(-1.2, 1), gr = rosen_gr)
  expect_false(ev@converged)
  expect_match(ev@criterion_met, "evaluation budget")
})


test_that("a step to a non-finite objective is refused, not taken", {
  # Propagating one NaN contaminates every iterate after it, and the run then
  # fails somewhere far from the cause.
  f <- function(p) if (sum(p^2) > 4) Inf else sum((p - c(1, 1))^2)
  res <- minimize(gd(step = 10, keep_trace = TRUE), f, c(0, 0))

  expect_true(all(is.finite(res@par)))
  expect_true(is.finite(res@value))
  expect_equal(res@par, c(1, 1), tolerance = 1e-5)
  expect_true(any(res@trace$safeguard != "none"))
})


test_that("the trace is a data frame or NULL, never an empty list", {
  f <- quad(c(1, 2))
  bare <- minimize(gd(), f, c(0, 0), gr = quad_gr(c(1, 2)))
  expect_null(bare@trace)

  kept <- minimize(gd(keep_trace = TRUE, criterion = crit_grad(1e-10)),
                   rosen, c(-1.2, 1), gr = rosen_gr)
  expect_s3_class(kept@trace, "data.frame")
  expect_equal(names(kept@trace),
               c("iteration", "value", "gnorm", "step", "safeguard"))
  expect_equal(nrow(kept@trace), kept@iterations)
  # the objective never increases: the Armijo condition guarantees it
  expect_true(all(diff(kept@trace$value) <= 0))
})


test_that("sufficient decrease is required, not mere non-increase", {
  # On a quadratic with step 1 the update x - 2(x - t) = 2t - x reflects through
  # the minimum and leaves the objective exactly unchanged. Accepting it lets
  # the iterate oscillate forever while a criterion watching the objective sees
  # no change and reports convergence at a point that is not a minimum.
  target <- c(1, 2)
  res <- minimize(gd(step = 1), quad(target), c(0, 0),
                  gr = quad_gr(target))
  expect_equal(res@par, target, tolerance = 1e-8)
  expect_true(res@value < 1e-12)
})


test_that("maximize() restores the sign of everything it flipped", {
  f <- function(p) -sum((p - c(1, 2))^2)
  res <- maximize(gd(keep_trace = TRUE), f, c(0, 0))
  expect_equal(res@par, c(1, 2), tolerance = 1e-6)
  expect_equal(res@value, 0, tolerance = 1e-10)
  expect_true(all(res@trace$value <= 1e-10))
  expect_true(all(diff(res@trace$value) >= 0))   # ascent, after the flip
})


test_that("constructors refuse settings that are not settings", {
  expect_error(gd(maxit = 0), "positive")
  expect_error(gd(refresh = -1), "non-negative")
  expect_error(gd(verbose = "yes"), "TRUE or FALSE")
  expect_error(gd(step = -1), "positive")
  expect_error(gd(line_search = armijo(shrink = 1)),
               "between 0 and 1")
  expect_error(gd(criterion = "gradient"), "criterion object")
  expect_error(crit_grad(0), "positive")
  expect_error(crit_any(), "At least one")
  expect_error(crit_any(crit_grad(), "nonsense"), "criterion")
})


test_that("verbose reports on the refresh schedule and nowhere else", {
  out <- capture.output(
    minimize(gd(verbose = TRUE, refresh = 3, maxit = 10,
                              criterion = crit_grad(1e-14)),
             rosen, c(-1.2, 1), gr = rosen_gr)
  )
  # a header, the iterations divisible by three, and a closing line
  rows <- grep("^\\s+\\d+\\s", out, value = TRUE)
  its <- as.integer(sub("^\\s*(\\d+).*$", "\\1", rows))
  expect_true(all(its %% 3 == 0))
  expect_true(any(grepl("stopped after", out)))

  # invisible(), or the returned object auto-prints and the capture is not empty
  quiet <- capture.output(
    invisible(minimize(gd(verbose = FALSE), quad(c(1, 2)),
                       c(0, 0), gr = quad_gr(c(1, 2))))
  )
  expect_length(quiet, 0)
})


test_that("format_elapsed picks the unit from the magnitude", {
  expect_equal(format_elapsed(5e-6), "5 us")
  expect_equal(format_elapsed(0.446), "446 ms")
  expect_equal(format_elapsed(1.79), "1.79 s")
  expect_equal(format_elapsed(150), "2 min 30 s")
  expect_equal(format_elapsed(4380), "1 h 13 min")
  expect_true(is.na(format_elapsed(numeric(0))))
  expect_true(is.na(format_elapsed(NA_real_)))
})


test_that("print shows the leading parameters, rounded, and the elapsed time", {
  f <- function(p) sum((p - seq_len(10))^2)
  r <- minimize(bfgs(), f, rep(0, 10), gr = function(p) 2 * (p - seq_len(10)))
  out <- capture.output(print(r))
  expect_true(any(grepl("(6 of 10 shown)", out, fixed = TRUE)))
  expect_true(any(grepl("elapsed", out)))
  # max_par and digits are honoured
  out2 <- capture.output(print(r, digits = 2, max_par = 3))
  expect_true(any(grepl("(3 of 10 shown)", out2, fixed = TRUE)))
  # nothing is truncated when the vector is short
  r2 <- minimize(bfgs(), function(p) sum((p - 1:2)^2), c(0, 0))
  expect_false(any(grepl("shown", capture.output(print(r2)))))

  expect_error(print(r, digits = -1), "'digits'")
  expect_error(print(r, max_par = 0), "'max_par'")
})


test_that("newton without derivatives warns when the budget cannot carry it", {
  f <- function(p) sum((p - seq_len(30))^2)
  # 4 * 30^2 = 3600 per iteration against a 2000-evaluation budget
  expect_warning(minimize(newton(max_eval = 2000, maxit = 3), f, rep(0, 30)),
                 "evaluations per iteration")
  # an analytic gradient removes the p^2 term and the warning with it
  expect_no_warning(minimize(newton(max_eval = 2000, maxit = 3), f, rep(0, 30),
                             gr = function(p) 2 * (p - seq_len(30))))
})


test_that("the evaluation budget defaults to Inf and binds only when finite", {
  for (o in list(gd(), cg(), bb(), bfgs(), lbfgs(), newton(), adam(),
                 nelder_mead(), compass(), bundle())) {
    expect_identical(o@max_eval, Inf, label = o@name)
  }
  # unlimited: the run ends on the criterion, not on any budget
  f <- function(p) sum((p - c(1, 2))^2)
  r <- minimize(bfgs(), f, c(0, 0), gr = function(p) 2 * (p - c(1, 2)))
  expect_true(r@converged)
  # a finite budget still stops the run and is reported as the reason
  r2 <- suppressWarnings(
    minimize(newton(max_eval = 50, maxit = 100),
             function(p) sum((p - seq_len(10))^2), rep(0, 10)))
  expect_false(r2@converged)
  expect_match(r2@message, "budget")
  # maxit stays the backstop, so it cannot be Inf
  expect_error(gd(maxit = Inf), "finite")
})


test_that("a gradient that does not belong to the objective draws a warning", {
  # The exact shape of the user error this exists for: an objective constant in
  # the parameter next to a gradient that is not.
  f_const <- function(p) 42
  g_lin   <- function(p) 2 * (p - c(1, 2))
  expect_warning(minimize(bfgs(maxit = 2), f_const, c(0, 0), gr = g_lin),
                 "does not appear to be the gradient")

  # a correct pair passes in silence
  f <- function(p) sum((p - c(1, 2))^2)
  expect_no_warning(minimize(bfgs(), f, c(0, 0), gr = g_lin))

  # and so does a subgradient of a non-smooth objective, the tolerance being
  # deliberately loose
  z <- c(-2, -1, 0.5, 1, 3)
  sad <- function(p) sum(abs(z - p))
  sub <- function(p) -sum(sign(z - p))
  expect_no_warning(minimize(bundle(), sad, 0.2, gr = sub))

  # the option turns it off
  old <- options(optimizers7.check_gradient = FALSE)
  on.exit(options(old), add = TRUE)
  expect_no_warning(minimize(bfgs(maxit = 2), f_const, c(0, 0), gr = g_lin))
  options(old)

  # a zero gradient at the start is not decidable and is not warned about
  expect_no_warning(minimize(bfgs(maxit = 2), f, c(1, 2),
                             gr = function(p) c(0, 0)))

  # Nor is a gradient merely NEGLIGIBLE against the size of the objective. A
  # caller who starts at the optimum is doing the best possible thing, and the
  # difference of fn along a direction of no slope is its own truncation error,
  # so the comparison has no signal. This warned before the guard existed.
  big <- function(p) 400 + 1e4 * sum((p - c(1, 2))^4)
  big_gr <- function(p) 4e4 * (p - c(1, 2))^3
  expect_no_warning(minimize(bfgs(maxit = 2), big, c(1, 2), gr = big_gr))
  # and the same pair away from the optimum still gets checked, so the guard
  # has not simply switched the check off
  expect_no_warning(minimize(bfgs(maxit = 2), big, c(3, 4), gr = big_gr))
  expect_warning(
    minimize(bfgs(maxit = 2), big, c(3, 4), gr = function(p) -big_gr(p)),
    "does not appear to be the gradient"
  )

  # multistart warns once, not once per start
  w <- 0
  withCallingHandlers(
    minimize(multistart(bfgs(maxit = 2), n = 5, ncores = 1),
             f_const, c(0, 0), gr = g_lin),
    warning = function(x) { w <<- w + 1; invokeRestart("muffleWarning") })
  expect_equal(w, 1)
})
