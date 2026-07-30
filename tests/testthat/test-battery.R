# check_optimizer(), the battery, and the defects the validator must catch.
#
# A validator nobody has tried to fool proves nothing, so half of this file is
# optimisers that lie in one specific way each. The other half confirms the
# honest ones still pass, which is what stops the checks from being satisfied by
# failing everything.


# --- the problems themselves ------------------------------------------------

test_that("every stated solution really is one", {
  for (p in test_problems()) {
    expect_equal(p$fn(p$solution), p$value, tolerance = 1e-12,
                 label = paste(p$name, "value at the solution"))
    expect_lt(max(abs(p$gr(p$solution))), 1e-8, label = paste(p$name, "gradient"))
  }
})


test_that("every analytic gradient matches the function it belongs to", {
  # The battery is the reference the checks are measured against, so it needs a
  # reference of its own. A central difference at a point that is NOT the
  # solution, since at the solution a wrong gradient could still come out zero.
  for (p in test_problems()) {
    if (!isTRUE(p$smooth)) next        # a difference across a kink is not one
    x <- p$par + 0.37
    h <- 1e-5
    fd <- vapply(seq_along(x), function(j) {
      xp <- x; xm <- x; xp[j] <- x[j] + h; xm[j] <- x[j] - h
      (p$fn(xp) - p$fn(xm)) / (2 * h)
    }, numeric(1))
    expect_equal(p$gr(x), fd, tolerance = 1e-5, label = p$name)
  }
})


test_that("test_problems() selects and refuses by name", {
  expect_named(test_problems("beale"), "beale")
  expect_length(test_problems(c("sphere", "booth")), 2L)
  expect_error(test_problems("banana"), "No such test problem")
})


# --- every shipped optimiser keeps the contract -----------------------------

all_optimizers <- function() {
  list(gradient_descent(), newton(), bfgs(), lbfgs(),
       adam(alpha = 0.05, maxit = 4000),
       nelder_mead(), compass(), compass(directions = "coordinate"),
       bundle(), multistart(bfgs(), n = 6))
}

test_that("every optimiser in the package passes all twelve checks", {
  set.seed(1)
  for (o in all_optimizers()) {
    res <- check_optimizer(o, problems = test_problems("sphere"),
                           verbose = FALSE)
    bad <- names(res$checks)[!vapply(res$checks, isTRUE, logical(1))]
    expect_equal(bad, character(), label = o@name)
  }
})


test_that("the report prints, and says which checks failed", {
  out <- capture.output(check_optimizer(bfgs(),
                                        problems = test_problems("sphere")))
  expect_true(any(grepl("Checking optimizer: BFGS", out)))
  expect_true(any(grepl("All checks passed", out)))
  expect_true(any(grepl("battery", out)))
  expect_false(any(grepl("FAILED", out)))
})


test_that("the battery separates power from correctness", {
  # BFGS solves the smooth unimodal problems and is defeated by the kink; that
  # is a true statement about the method, and the point of reporting the gaps
  # as information rather than as a verdict is that it does not become a
  # failure.
  res <- check_optimizer(bfgs(), verbose = FALSE)
  b <- res$battery
  smooth <- b$problem %in% c("sphere", "rosenbrock", "booth", "beale", "powell")
  expect_true(all(b$gap[smooth] < 1e-10))
  expect_true(all(res$checks))

  # the kink defeats it, and bundle does not mind it
  expect_gt(b$gap[b$problem == "abs_sum"], 1e-6)
  bb <- check_optimizer(bundle(), problems = test_problems("abs_sum"),
                        verbose = FALSE)
  expect_lt(bb$battery$gap, 1e-6)
})


test_that("a local method finding another local minimum is not a failure", {
  res <- check_optimizer(bfgs(), problems = test_problems("rastrigin"),
                         verbose = FALSE)
  expect_gt(res$battery$gap, 1)        # it found a different minimum
  expect_true(all(res$checks))         # and behaved correctly throughout
})


test_that("bundle takes serious steps on a steep problem", {
  # A regression, and the battery is what found it. With t0 read as a bare
  # multiplier the first step was as long as the gradient was big -- 233 units on
  # rosenbrock -- landing where the objective is 1e11; the gradient then squared
  # faster than halving t could shrink it, and the run spent 3000 iterations
  # taking rejected steps and returned the point it started from.
  #
  # The signature of that failure is ZERO serious steps, so that is what is
  # asserted, alongside the value. Testing the value alone would let a version
  # that limped to a poor answer pass.
  for (nm in c("rosenbrock", "beale", "powell")) {
    p <- test_problems(nm)[[1]]
    r <- minimize(bundle(maxit = 2000), p$fn, p$par, gr = p$gr)
    expect_true(r@converged, label = nm)
    expect_lt(r@value - p$value, 1e-4)
    serious <- as.integer(sub(" serious.*", "", r@message))
    expect_gt(serious, 0)
  }
})


test_that("a diverging model is reported rather than run to the budget", {
  # The scaling stops the runaway from starting; it cannot stop an objective
  # from overflowing when told explicitly to take an enormous first step. The
  # trial at -1e70 has a finite value and a subgradient of 4e210, whose square
  # is not representable, so the subproblem's matrix holds an infinity and every
  # iteration after it is wasted. Saying so beats spending the budget.
  f <- function(p) p^4 + p
  g <- function(p) 4 * p^3 + 1
  r <- minimize(bundle(t0 = 1e70, t_max = 1e70, maxit = 500), f, par = 0.5,
                gr = g)
  expect_false(r@converged)
  expect_lt(r@iterations, 500)
  expect_match(r@message, "diverged")
})


test_that("a run that can no longer shrink its trust region stops saying so", {
  # An objective finite nowhere but the point it started from. Every trial is
  # rejected and t is halved, until t reaches its floor -- at which point the
  # only control the method has is exhausted and every further iteration would
  # re-propose the identical rejected step. The trigger needs no chosen
  # constant, which is why it is this one and not a count of retries.
  # Exact equality, not all.equal: its default tolerance is 1.5e-8, so once the
  # step shrank below that the objective came back finite and the run took null
  # steps instead of being refused -- testing the wrong branch entirely.
  x0 <- 0.7
  f <- function(p) if (p == x0) 1 else NaN
  g <- function(p) 1
  r <- minimize(bundle(maxit = 500), f, par = x0, gr = g)
  expect_false(r@converged)
  expect_lt(r@iterations, 500)
  expect_match(r@message, "smallest trust region")
  expect_equal(r@par, x0)
})


# --- optimisers that lie -----------------------------------------------------

# Each one runs BFGS honestly and then corrupts exactly one thing it reports,
# so that a check failing identifies the lie rather than a general malaise.
Liar <- S7::new_class("Liar", parent = optimizer,
                      properties = list(lie = S7::class_character))

S7::method(minimize, Liar) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    r <- minimize(bfgs(criterion = optimizer@criterion,
                       maxit = optimizer@maxit,
                       keep_trace = optimizer@keep_trace),
                  fn, par, gr = gr, he = he, bounds = bounds)
    switch(optimizer@lie,
      value     = r@value <- r@value + 1,
      converged = { r@converged <- TRUE
                    r@criterion_met <- optimizer@criterion@label },
      budget    = r@iterations <- optimizer@maxit + 5,
      bounds    = r@par <- rep(5, length(r@par)),
      gradient  = r@gradient <- r@gradient + 0.5,
      counts    = r@counts <- c(f = 0, g = 0, h = 0),
      stop("unknown lie")
    )
    r
  }

liar <- function(lie, ...) {
  Liar(name = paste0("liar (", lie, ")"), criterion = crit_grad(1e-8),
       maxit = 200, max_eval = 10000, verbose = FALSE, refresh = 10,
       keep_trace = FALSE, lie = lie, ...)
}

fails_only <- function(lie, which) {
  res <- check_optimizer(liar(lie), problems = test_problems("sphere"),
                         verbose = FALSE)
  bad <- names(res$checks)[!vapply(res$checks, isTRUE, logical(1))]
  expect_true(which %in% bad, label = paste0(lie, ": ", which, " must fail"))
}

test_that("a corrupted value is caught", {
  fails_only("value", "value agrees with par")
})

test_that("a gradient that is not the gradient at par is caught", {
  fails_only("gradient", "gradient agrees with par")
})

test_that("claiming convergence when the run was cut off is caught", {
  fails_only("converged", "convergence is not assumed")
})

test_that("exceeding the iteration budget is caught", {
  fails_only("budget", "budgets are respected")
})

test_that("returning a point ON its bound is caught", {
  fails_only("bounds", "bounds are respected strictly")
})

test_that("reporting no evaluations is caught", {
  fails_only("counts", "evaluations are counted")
})


test_that("an honest optimiser of the same shape passes, so the checks are not trivial", {
  # Without this, every check above could be satisfied by a validator that
  # simply always fails.
  # A class of its own rather than another method on Liar: re-registering one
  # would emit "Overwriting method" and, worse, leave the corrupting version
  # replaced for whatever ran after it.
  Honest <- S7::new_class("Honest", parent = optimizer)
  S7::method(minimize, Honest) <-
    function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
      minimize(bfgs(criterion = optimizer@criterion, maxit = optimizer@maxit,
                    keep_trace = optimizer@keep_trace),
               fn, par, gr = gr, he = he, bounds = bounds)
    }
  honest <- Honest(name = "honest", criterion = crit_grad(1e-8), maxit = 200,
                   max_eval = 10000, verbose = FALSE, refresh = 10,
                   keep_trace = FALSE)
  res <- check_optimizer(honest, problems = test_problems("sphere"),
                         verbose = FALSE)
  expect_true(all(res$checks))
})


test_that("check_optimizer refuses something that is not an optimiser", {
  expect_error(check_optimizer(bfgs), "must be an optimizer")
  expect_error(check_optimizer("bfgs"), "must be an optimizer")
})


test_that("with_criterion reaches the rule that is actually consulted", {
  # multistart carries a criterion only so that printing it tells the truth; the
  # rule evaluated belongs to the optimiser inside. Setting the outer one and
  # expecting a different run is how a check passes while testing nothing.
  m <- multistart(nelder_mead(), n = 3)
  m2 <- with_criterion(m, crit_grad())
  expect_identical(m2@optimizer@criterion@label, crit_grad()@label)
  expect_error(minimize(m2, function(p) sum(p^2), c(1, 1)), "needs gradient")
})
