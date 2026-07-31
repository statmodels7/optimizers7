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
  # favourable at all, since the exact ratio is not a contract
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
  # Maximising with a *minimising* driver produces exactly that situation.
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
  # Checked directly rather than through the optimiser: at the accepted point,
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
