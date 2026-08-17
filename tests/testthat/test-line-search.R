# The line search is shared machinery: Newton, BFGS and L-BFGS will all call it,
# so what it guarantees has to be true before any of them is written.

rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) {
  c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
    200 * (p[2] - p[1]^2))
}


test_that("both searches reach the same minimum of a hard function", {
  x0 <- c(-1.2, 1)
  a <- minimize(gd(line_search = armijo(), maxit = 50000,
                                 max_eval = 1e6, criterion = crit_grad(1e-8)),
                rosen, x0, gr = rosen_gr)
  w <- minimize(gd(line_search = wolfe(), maxit = 50000,
                                 max_eval = 1e6, criterion = crit_grad(1e-8)),
                rosen, x0, gr = rosen_gr)

  expect_true(a@converged)
  expect_true(w@converged)
  expect_equal(a@par, c(1, 1), tolerance = 1e-5)
  expect_equal(w@par, c(1, 1), tolerance = 1e-5)
})


test_that("the objective never increases, whichever search is used", {
  for (ls in list(armijo(), wolfe())) {
    res <- minimize(gd(line_search = ls, keep_trace = TRUE,
                                     maxit = 300, criterion = crit_grad(1e-10)),
                    rosen, c(-1.2, 1), gr = rosen_gr)
    expect_true(all(diff(res@trace$value) <= 0),
                label = paste("monotone under", ls@label))
  }
})


test_that("the two searches cost what they are supposed to cost", {
  x0 <- c(-1.2, 1)
  o <- function(ls) gd(line_search = ls, maxit = 50000,
                                     max_eval = 1e7, criterion = crit_grad(1e-8))
  a <- minimize(o(armijo()), rosen, x0, gr = rosen_gr)
  w <- minimize(o(wolfe()),  rosen, x0, gr = rosen_gr)

  # Armijo evaluates no gradient of its own, and the loop computes exactly one
  # per iteration: the gradient at the accepted point serves both the stopping
  # rule and the next iteration's direction. One more is spent before the loop
  # starts. That is an exact invariant, and it breaking means something began
  # recomputing a gradient it already had -- which is precisely the waste the
  # shared loop was written to remove.
  expect_equal(a@counts[["g"]], a@iterations + 1)

  # Wolfe evaluates gradients inside the search, so it costs slightly more per
  # iteration -- and needs fewer of them, which is the whole trade.
  expect_gt(w@counts[["g"]] / w@iterations, 1)
  expect_lt(w@counts[["g"]] / w@iterations, 2)
  expect_lt(w@iterations, a@iterations)

  # measured 0.86 on this problem; the assertion is only that the trade is
  # favorable at all, since the exact ratio is not a contract
  work <- function(r) r@counts[["f"]] + r@counts[["g"]]
  expect_lt(work(w), work(a))
})


test_that("Armijo demands sufficient decrease, not mere non-increase", {
  # The reflection through the minimum: x - 2(x - t) leaves f exactly unchanged.
  target <- c(1, 2)
  res <- minimize(gd(step = 1, line_search = armijo()),
                  function(p) sum((p - target)^2), c(0, 0),
                  gr = function(p) 2 * (p - target))
  expect_equal(res@par, target, tolerance = 1e-8)
})


test_that("a direction that is not one of descent is reported, not searched", {
  # No step along an ascent direction can decrease the objective, so a search
  # that merely shrank would burn its whole budget before admitting it.
  # Maximizing with a *minimizing* driver produces exactly that situation.
  f <- function(p) -sum(p^2)          # unbounded below along -g from anywhere
  res <- minimize(gd(step = 1, keep_trace = TRUE, maxit = 5),
                  f, c(1, 1), gr = function(p) -2 * p)
  # it must not claim to have converged
  expect_false(res@converged)
})


test_that("a non-finite region is stepped around, not into", {
  f <- function(p) if (sum(p^2) > 4) Inf else sum((p - c(1, 1))^2)
  for (ls in list(armijo(), wolfe())) {
    res <- minimize(gd(step = 50, line_search = ls,
                                     keep_trace = TRUE), f, c(0, 0))
    expect_true(is.finite(res@value), label = ls@label)
    expect_true(all(is.finite(res@par)), label = ls@label)
    expect_equal(res@par, c(1, 1), tolerance = 1e-4, label = ls@label)
  }
})


test_that("the search reports which safeguard it used", {
  f <- function(p) if (sum(p^2) > 4) Inf else sum((p - c(1, 1))^2)
  res <- minimize(gd(step = 50, keep_trace = TRUE), f, c(0, 0))
  expect_true(any(res@trace$safeguard != "none"))
  expect_true(all(res@trace$safeguard %in%
                    c("none", "step shortened", "objective non-finite",
                      "step adjusted", "curvature condition not met",
                      "not a descent direction", "no acceptable step")))
})


test_that("line searches are objects with validated settings", {
  expect_s3_class(armijo(), "S7_object")
  expect_s3_class(wolfe(), "S7_object")
  expect_match(armijo()@label, "Armijo")
  expect_match(wolfe()@label, "Wolfe")

  expect_error(armijo(c1 = 0), "between 0 and 1")
  expect_error(armijo(c1 = 1), "between 0 and 1")
  expect_error(armijo(shrink = 2), "between 0 and 1")
  expect_error(armijo(max_step = 0), "positive")
  expect_error(wolfe(c1 = 0.9, c2 = 0.1), "greater than")
  expect_error(gd(line_search = "armijo"), "line_search object")
})


test_that("the spec handed to C++ has one shape whichever search it describes", {
  a <- line_search_spec(armijo(c1 = 1e-3, shrink = 0.3, max_step = 12))
  w <- line_search_spec(wolfe(c1 = 1e-3, c2 = 0.4, max_step = 12))
  expect_equal(names(a), names(w))
  expect_equal(a$type, "armijo")
  expect_equal(w$type, "wolfe")
  expect_equal(a$c1, 1e-3)
  expect_equal(w$c2, 0.4)
  expect_type(a$max_step, "integer")
})


test_that("Wolfe's curvature condition actually holds where it claims", {
  # Checked directly rather than through the optimizer: at the accepted point,
  # |g'd| must have shrunk by the factor c2. This is the property BFGS will
  # depend on, so it is worth asserting on its own.
  x <- c(-1.2, 1)
  g <- rosen_gr(x)
  d <- -g
  c1 <- 1e-4; c2 <- 0.9
  res <- minimize(gd(line_search = wolfe(c1 = c1, c2 = c2),
                                   maxit = 1, criterion = crit_grad(1e-14)),
                  rosen, x, gr = rosen_gr)
  s <- res@trace$step[1]
  if (is.null(s)) s <- NA_real_

  res_tr <- minimize(gd(line_search = wolfe(c1 = c1, c2 = c2),
                                      maxit = 1, keep_trace = TRUE,
                                      criterion = crit_grad(1e-14)),
                     rosen, x, gr = rosen_gr)
  s <- res_tr@trace$step[1]
  x_new <- x + s * d
  expect_lte(rosen(x_new), rosen(x) + c1 * s * sum(g * d))          # Armijo
  expect_lte(abs(sum(rosen_gr(x_new) * d)), c2 * abs(sum(g * d)))   # curvature
})


test_that("a search that cannot move at a stationary point reports convergence", {
  # A caller that starts from a closed-form estimate hands the optimizer the
  # answer, and every line search then rejects every step. That is not a
  # failure: the stopping rule is asked before the run gives up, with no
  # previous iterate, so only the state at the point can end it.
  f  <- function(p) sum((p - c(1, 2))^2)
  gr <- function(p) 2 * (p - c(1, 2))

  for (o in list(bfgs(), lbfgs(), cg(), gd(), newton(), bb())) {
    r <- minimize(o, f, c(1, 2), gr = gr)
    expect_true(r@converged)
    expect_identical(r@message, "")
    expect_true(nzchar(r@criterion_met))
  }
})

test_that("a search that cannot move away from a stationary point still fails", {
  # The counterexample: a mis-stated gradient makes the direction ascend, so
  # no step is acceptable at a point that is nowhere near stationary. The
  # gradient rule cannot fire there and the failure must survive.
  f <- function(p) sum((p - c(1, 2))^2)
  r <- suppressWarnings(
    minimize(bfgs(), f, c(-3, 4), gr = function(p) -2 * (p - c(1, 2)))
  )
  expect_false(r@converged)
  expect_match(r@message, "no acceptable step")
})


test_that("a resolution stops a search that cannot verify a smaller decrease", {
  # An objective computed by a procedure returns slightly different values for
  # the same argument. Below that spread its values carry no information, and
  # the Armijo demand c1 * s * |g'd| shrinks with the step, so every remaining
  # backtrack compares numbers the objective cannot tell apart.
  mk <- function(jitter) {
    A <- diag(c(1, 12))
    k <- 0L
    list(
      fn = function(z) {
        k <<- k + 1L
        3e4 + 0.5 * drop(t(z) %*% A %*% z) +
          jitter * sin(k * 2.399963229728653)
      },
      gr = function(z) drop(A %*% z))
  }

  # with no resolution declared the run reports failure at a point it is not
  # leaving, which is the behaviour this exists to correct
  o <- mk(1e-3)
  bare <- suppressWarnings(minimize(lbfgs(), o$fn, c(3, -2), gr = o$gr))
  expect_false(bare@converged)

  o <- mk(1e-3)
  told <- suppressWarnings(
    minimize(lbfgs(line_search = wolfe(resolution = 1e-3)), o$fn, c(3, -2),
             gr = o$gr))
  expect_true(told@converged)
  expect_match(told@criterion_met, "resolution")
  # the same ANSWER, reached for less. Both runs end where the objective stops
  # carrying information, so they agree to the size of the jitter and not to
  # machine precision, which is the whole point of declaring one
  expect_equal(told@value, bare@value, tolerance = 1e-3)
  expect_true(max(abs(told@par)) < 1e-2)
  expect_lt(told@counts[["f"]], bare@counts[["f"]])
})


test_that("a declared resolution does not rescue a genuine failure", {
  # THE CONTROL THAT REFUTED THE FIRST DESIGN. A mis-stated gradient makes the
  # direction ascend at a point nowhere near stationary. With the resolution
  # asked inside the backtracking loop this run was promoted to converged,
  # because a vanishing step stops resolving the change whichever of the two
  # is wrong. Asked at the full step the predicted improvement is 10 against a
  # resolution of 1e-3, so the failure survives.
  f <- function(p) sum((p - c(1, 2))^2)
  for (res in c(1e-3, 1e-1)) {
    r <- suppressWarnings(
      minimize(bfgs(line_search = wolfe(resolution = res)), f, c(-3, 4),
               gr = function(p) -2 * (p - c(1, 2)))
    )
    expect_false(r@converged)
    expect_match(r@message, "no acceptable step")
  }
})


test_that("the resolution is off by default and validated", {
  expect_identical(armijo()@resolution, 0)
  expect_identical(wolfe()@resolution, 0)
  expect_identical(nonmonotone()@resolution, 0)
  expect_identical(line_search_spec(armijo(resolution = 2))$resolution, 2)
  expect_identical(line_search_spec(wolfe(resolution = 2))$resolution, 2)
  expect_identical(line_search_spec(nonmonotone(resolution = 2))$resolution, 2)
  for (bad in list(-1, c(1, 2), NA_real_, Inf, "x")) {
    expect_error(armijo(resolution = bad), "non-negative finite")
    expect_error(wolfe(resolution = bad), "non-negative finite")
    expect_error(nonmonotone(resolution = bad), "non-negative finite")
  }
  # and with it at zero every method reaches the same answer it always did
  f <- function(p) sum((p - c(1, 2))^2)
  gr <- function(p) 2 * (p - c(1, 2))
  a <- minimize(bfgs(line_search = wolfe()), f, c(-3, 4), gr = gr)
  b <- minimize(bfgs(line_search = wolfe(resolution = 0)), f, c(-3, 4), gr = gr)
  expect_identical(a@par, b@par)
  expect_identical(a@counts, b@counts)
})


test_that("armijo takes a resolution on its own backtracking", {
  # bb() runs a nonmonotone Armijo search, so the armijo branch of the kernel
  # is the one exercised here rather than wolfe's zoom.
  mk <- function(jitter) {
    k <- 0L
    list(fn = function(z) {
           k <<- k + 1L
           3e4 + 0.5 * sum(c(1, 12) * z^2) + jitter * sin(k * 2.3999632297)
         },
         gr = function(z) c(1, 12) * z)
  }
  o <- mk(1e-3)
  bare <- suppressWarnings(
    minimize(gd(line_search = armijo(), maxit = 200L), o$fn, c(3, -2),
             gr = o$gr))
  o <- mk(1e-3)
  told <- suppressWarnings(
    minimize(gd(line_search = armijo(resolution = 1e-3), maxit = 200L),
             o$fn, c(3, -2), gr = o$gr))
  expect_true(told@converged)
  expect_match(told@criterion_met, "resolution")
  expect_lt(told@counts[["f"]], bare@counts[["f"]])
})


test_that("a resolution that moves is asked again at every iteration", {
  # An objective that settles as the run goes -- a fit warm-started from the
  # previous evaluation locates its answer better each time -- has a
  # resolution that improves with it, and the reading taken once at the start
  # is the one from the worst point of the whole run.
  mk <- function() {
    A <- diag(c(1, 12)); k <- 0L
    list(fn = function(z) {
           k <<- k + 1L
           3e4 + 0.5 * drop(t(z) %*% A %*% z) + 1e-3 * sin(k * 2.3999632297)
         },
         gr = function(z) drop(A %*% z))
  }

  asked <- 0L
  moving <- function() { asked <<- asked + 1L; 1e-3 }
  o <- mk()
  r <- suppressWarnings(
    minimize(lbfgs(line_search = wolfe(resolution = moving)), o$fn, c(3, -2),
             gr = o$gr))
  expect_true(r@converged)
  expect_match(r@criterion_met, "resolution")
  # once per invocation of the search, so once an iteration and not once per
  # trial: the count tracks the iterations rather than the evaluations
  expect_gt(asked, 0L)
  expect_lte(asked, r@iterations + 1L)

  # a constant function and the number itself are the same run
  o <- mk()
  a <- suppressWarnings(
    minimize(lbfgs(line_search = wolfe(resolution = 1e-3)), o$fn, c(3, -2),
             gr = o$gr))
  o <- mk()
  b <- suppressWarnings(
    minimize(lbfgs(line_search = wolfe(resolution = function() 1e-3)), o$fn,
             c(3, -2), gr = o$gr))
  expect_identical(a@par, b@par)
  expect_identical(a@counts, b@counts)
  expect_identical(a@converged, b@converged)

  # what is not finite and positive asks nothing, which is how a caller with
  # no reading yet says so without a branch of its own
  o <- mk()
  none <- suppressWarnings(
    minimize(lbfgs(line_search = wolfe(resolution = function() NA_real_)),
             o$fn, c(3, -2), gr = o$gr))
  expect_false(none@converged)

  # a function that takes arguments is refused at construction, where the
  # caller can see it, rather than at the first iteration
  expect_error(wolfe(resolution = function(x) 1), "no arguments")
  expect_error(armijo(resolution = function(x) 1), "no arguments")
  expect_error(nonmonotone(resolution = function(x) 1), "no arguments")
})
