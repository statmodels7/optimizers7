rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                          200 * (p[2] - p[1]^2))
rosen_he <- function(p) matrix(c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
                                 -400 * p[1], 200), 2, 2)
rosen_t3 <- function(p, d) c(2400 * p[1] * d[1]^2 - 800 * d[1] * d[2],
                             -400 * d[1]^2)

test_that("the transcribed third derivative is a difference of the Hessian", {
  x <- c(0.3, -0.7); d <- c(0.4, 1.1); h <- 1e-5
  ref <- as.numeric((rosen_he(x + h * d) - rosen_he(x - h * d)) %*% d) / (2 * h)
  expect_lt(max(abs(rosen_t3(x, d) - ref)), 1e-6)
})

test_that("chebyshev() reaches the minimum in fewer iterations than newton()", {
  for (st in list(c(-1.2, 1), c(2, 2), c(0.8, 0.6))) {
    a <- minimize(newton(), rosen, st, gr = rosen_gr, he = rosen_he)
    b <- minimize(chebyshev(), rosen, st, gr = rosen_gr, he = rosen_he,
                  t3 = rosen_t3)
    expect_true(b@converged)
    expect_equal(b@par, c(1, 1), tolerance = 1e-8)
    expect_lt(b@iterations, a@iterations)
  }
})

test_that("t3 is called, and a second difference of the gradient serves without it", {
  calls <- 0L
  t3c <- function(p, d) { calls <<- calls + 1L; rosen_t3(p, d) }
  b <- minimize(chebyshev(), rosen, c(2, 2), gr = rosen_gr, he = rosen_he,
                t3 = t3c)
  expect_gt(calls, 0L)
  f <- minimize(chebyshev(), rosen, c(2, 2), gr = rosen_gr, he = rosen_he)
  expect_identical(f@iterations, b@iterations)
  expect_equal(f@par, b@par, tolerance = 1e-10)
  # the difference costs two gradients an iteration
  expect_gt(f@counts[["g"]], b@counts[["g"]])
})

test_that("a correction that is too long is refused and recorded", {
  b <- minimize(chebyshev(keep_trace = TRUE), rosen, c(-1.2, 1),
                gr = rosen_gr, he = rosen_he, t3 = rosen_t3)
  expect_true(any(b@trace$safeguard == "cubic correction refused"))
  # a tiny ratio refuses every correction and reproduces newton()
  n <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
  t <- minimize(chebyshev(ratio = 1e-12), rosen, c(-1.2, 1), gr = rosen_gr,
                he = rosen_he, t3 = rosen_t3)
  expect_identical(t@iterations, n@iterations)
  expect_identical(t@par, n@par)
})

test_that("arguments are validated", {
  expect_error(chebyshev(ratio = 0), "'ratio'")
  expect_error(chebyshev(ratio = c(1, 2)), "'ratio'")
  expect_error(minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr,
                        he = rosen_he, t3 = 1), "'t3'")
  expect_error(minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr,
                        he = rosen_he, t3 = rosen_t3, lower = -5),
               "unconstrained scale")
  # bounds without t3 run on the difference
  r <- minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he,
                lower = -5, upper = 5)
  expect_true(r@converged)
  expect_equal(r@par, c(1, 1), tolerance = 1e-6)
  expect_error(minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr,
                        he = rosen_he, t3 = function(p, d) 1),
               "as long as the parameters")
})

test_that("chebyshev() passes check_optimizer()", {
  r <- check_optimizer(chebyshev(), verbose = FALSE)
  expect_true(all(r$checks))
})
