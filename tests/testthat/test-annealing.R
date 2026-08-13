# Simulated annealing, and the composition that hands its answer on.
#
# What the checks below are about is the CONTRACT -- the point reported is the
# best seen, convergence is not inferred from the schedule, the run repeats --
# rather than the power, which is the battery's business and is printed rather
# than asserted.

rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
sphere <- function(p) sum(p^2)

test_that("it satisfies the optimizer contract", {
  set.seed(6)
  out <- check_optimizer(sa(), verbose = FALSE)
  expect_true(all(out$checks),
              info = paste(names(out$checks)[!out$checks], collapse = ", "))
  set.seed(7)
  out2 <- check_optimizer(sa(visiting = "cauchy"), verbose = FALSE)
  expect_true(all(out2$checks),
              info = paste(names(out2$checks)[!out2$checks], collapse = ", "))
})

test_that("the point reported is the best seen, not the last", {
  # A walk that accepts uphill moves ends wherever it happens to be, so the
  # two are different points and reporting the wrong one is invisible unless
  # the value is checked against the parameter.
  set.seed(1)
  r <- minimize(sa(), rastrigin, c(4.4, -3.6))
  expect_equal(r@value, rastrigin(r@par), tolerance = 1e-12)
  # and it is not worse than where it started
  expect_lt(r@value, rastrigin(c(4.4, -3.6)))
})

test_that("convergence is not inferred from the schedule finishing", {
  # One temperature level cannot have settled anything, and a run that merely
  # exhausted its budget must say so.
  set.seed(2)
  r <- minimize(sa(maxit = 1), sphere, c(2, 2))
  expect_false(r@converged)
  expect_equal(r@iterations, 1)
  expect_match(r@criterion_met, "budget")
})

test_that("the stopping rule is Corana's and reads the walk, not the best", {
  # The best-so-far stops improving long before the walk settles. Measured on
  # the version that read it: the run stopped at the seventh level whatever
  # its budget, and returned a point two tenths from the answer. The rule
  # reads the value the walk ends each level at, plus its distance from the
  # best, so a longer budget buys a tighter answer.
  set.seed(1)
  short <- minimize(sa(maxit = 15), sphere, c(3, 3, 3))
  set.seed(1)
  long <- minimize(sa(maxit = 100), sphere, c(3, 3, 3))
  expect_gt(long@counts[["f"]], short@counts[["f"]])
  expect_lt(long@value, short@value / 100)
  expect_lt(max(abs(long@par)), 1e-3)
})

test_that("the step adapts per coordinate", {
  # Coordinates on scales orders of magnitude apart is the case the adaptation
  # exists for: one step length cannot serve a parameter of size 1000 beside
  # one of size 0.001, and a method that used one would leave the second
  # untouched.
  f <- function(p) (p[1] - 1000)^2 / 1e6 + (p[2] - 0.001)^2 * 1e6
  set.seed(3)
  r <- minimize(sa(), f, c(0, 0))
  expect_lt(abs(r@par[1] - 1000), 5)
  expect_lt(abs(r@par[2] - 0.001), 1e-3)
})

test_that("the temperature is calibrated to the objective's own scale", {
  # The same problem scaled by a million must be solved as well, which a fixed
  # initial temperature cannot do: it would accept everything at one scale and
  # nothing at the other.
  for (s in c(1e-6, 1, 1e6)) {
    set.seed(4)
    r <- minimize(sa(), function(p) s * sum((p - 2)^2), c(-3, 5))
    expect_lt(max(abs(r@par - 2)), 1e-2, label = paste("scale", s))
  }
  # and a temperature the caller sets is used instead
  set.seed(4)
  expect_silent(minimize(sa(t0 = 10, maxit = 5), sphere, c(1, 1)))
  expect_error(sa(t0 = -1), "positive number")
})

test_that("a heavy-tailed proposal is available and is not the default", {
  expect_identical(sa()@visiting, "uniform")
  expect_identical(sa(visiting = "cauchy")@visiting, "cauchy")
  expect_match(sa(visiting = "cauchy")@name, "cauchy")
  # both reach the answer on a unimodal problem
  for (v in c("uniform", "cauchy")) {
    set.seed(5)
    r <- minimize(sa(visiting = v), sphere, c(4, -4))
    expect_lt(max(abs(r@par)), 1e-2, label = v)
  }
})

test_that("a rule it cannot evaluate is refused, and it says which", {
  expect_identical(optimizer_provides(sa()), "stationarity")
  expect_error(minimize(sa(criterion = crit_grad()), sphere, c(1, 1)),
               "does not provide")
})

test_that("the arguments are validated where they are written", {
  expect_error(sa(cooling = 1), "in \\(0, 1\\)")
  expect_error(sa(cooling = 0), "in \\(0, 1\\)")
  expect_error(sa(cycles = 0), "positive whole number")
  expect_error(sa(steps = 2.5), "positive whole number")
  expect_error(sa(target_accept = 0.95), "in \\(0.1, 0.9\\)")
  expect_error(sa(adjust = 0), "positive number")
  expect_error(sa(step = -1), "positive number")
})

test_that("the run repeats from the seed it recorded", {
  set.seed(11)
  a <- minimize(sa(maxit = 10), rastrigin, c(2, 2))
  expect_false(is.null(a@seed))
  assign(".Random.seed", a@seed, envir = globalenv())
  b <- minimize(sa(maxit = 10), rastrigin, c(2, 2))
  expect_equal(a@par, b@par)
  expect_equal(a@value, b@value)
})


test_that("the compiled loop and the R twin are the same run", {
  # They draw from R's generator in the same order, so from one seed this is
  # not a comparison at a tolerance but the same run twice: what it pins is
  # the order of the draws, and above all that the Metropolis test consumes a
  # uniform ONLY on an uphill proposal, which is where a transcription drifts.
  cheap <- function(p) sum(p^2) + 3 * sum(sin(p))
  for (cau in c(FALSE, TRUE)) {
    set.seed(21)
    a <- minimize(sa(maxit = 12, criterion = crit_stationary(1e-300),
                     visiting = if (cau) "cauchy" else "uniform"),
                  cheap, c(2, -1, 0.5))
    set.seed(21)
    b <- sa_run_r(cheap, c(2, -1, 0.5), cauchy = cau, maxit = 12)
    lab <- if (cau) "cauchy" else "uniform"
    expect_equal(as.numeric(a@par), b$par, info = lab)
    expect_equal(a@value, b$value, info = lab)
    # the same number of evaluations, which says neither took a branch the
    # other did not
    expect_identical(as.numeric(a@counts[["f"]]), as.numeric(b$n_value),
                     info = lab)
  }
})


# --- chain -----------------------------------------------------------------

test_that("a chain runs its stages in order, each from the last", {
  set.seed(3)
  a <- minimize(sa(maxit = 20), rastrigin, c(3.5, -2.5))
  set.seed(3)
  ch <- minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))
  # the second stage started where the first finished, so it can only have
  # improved on it
  expect_lte(ch@value, a@value)
  # and the work of both is reported
  expect_gt(ch@counts[["f"]], a@counts[["f"]])
  expect_match(ch@optimizer@name, "then")
})

test_that("a chain of one is that stage", {
  set.seed(4); a <- minimize(chain(bfgs()), rastrigin, c(0.4, 0.3))
  set.seed(4); b <- minimize(bfgs(), rastrigin, c(0.4, 0.3))
  expect_equal(a@par, b@par)
  expect_equal(a@value, b@value)
  expect_identical(a@converged, b@converged)
  expect_identical(a@counts, b@counts)
})

test_that("the verdict is the last stage's", {
  # A first stage that exhausts its budget is how a global search ordinarily
  # ends and says nothing about the whole; what decides is the method that
  # finished.
  set.seed(5)
  ch <- minimize(chain(sa(maxit = 1), bfgs()), sphere, c(2, 2))
  expect_true(ch@converged)
  expect_identical(ch@criterion_met, bfgs()@criterion@label)
  # and the other way round: a last stage starved reports failure
  set.seed(5)
  ch2 <- minimize(chain(bfgs(), sa(maxit = 1)), sphere, c(2, 2))
  expect_false(ch2@converged)
})

test_that("a chain reports what it can promise", {
  # provides is the last stage's, since the result comes from there
  expect_identical(optimizer_provides(chain(sa(), bfgs())), "gradient")
  expect_identical(optimizer_provides(chain(bfgs(), sa())), "stationarity")
  # bounded only when every stage is: they all receive the bounds, and a
  # proximal method takes its constraint inside the operator instead
  expect_true(optimizer_bounded(chain(sa(), bfgs())))
  soft <- prox_grad(prox = function(v, t) v, g = function(b) 0)
  expect_false(optimizer_bounded(chain(sa(), soft)))
  # a rule the last stage cannot evaluate is refused before anything runs
  expect_error(minimize(chain(bfgs(), sa(criterion = crit_grad())),
                        sphere, c(1, 1)), "does not provide")
})

test_that("a chain takes optimizers and says so when it does not", {
  expect_error(chain(), "at least one")
  expect_error(chain(bfgs(), "lbfgs"), "must be an optimizer")
  expect_error(chain(bfgs(), verbose = 1), "TRUE or FALSE")
})

test_that("a chain composes with multistart", {
  set.seed(8)
  r <- minimize(multistart(chain(sa(maxit = 5), bfgs()), n = 3),
                rastrigin, c(3.5, -2.5))
  expect_true(is.finite(r@value))
  expect_equal(r@value, rastrigin(r@par), tolerance = 1e-8)
})

test_that("the traces of the stages are stacked when they agree", {
  set.seed(9)
  r <- minimize(chain(sa(maxit = 5, keep_trace = TRUE),
                      sa(maxit = 5, keep_trace = TRUE), keep_trace = TRUE),
                sphere, c(2, 2))
  expect_true(is.data.frame(r@trace))
  expect_true("stage" %in% names(r@trace))
  expect_identical(sort(unique(r@trace$stage)), c(1L, 2L))
})
