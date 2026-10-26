# The Hessian modification objects of newton() and chebyshev().

saddle    <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
saddle_gr <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
saddle_he <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)

run_saddle <- function(opt, start = c(0.05, 0.6)) {
  minimize(opt, saddle, start, gr = saddle_gr, he = saddle_he)
}


test_that("the constructors build objects and reject invalid floors", {
  expect_s7_class(modify_eigen(), optimizers7:::ModifyEigen)
  expect_s7_class(modify_ridge(floor = 1e-6), optimizers7:::ModifyRidge)
  expect_s7_class(modify_eigen(), hessian_modification)
  expect_equal(modify_eigen()@floor, 1e-8)
  expect_equal(modify_ridge(floor = 1e-6)@floor, 1e-6)
  for (bad in list(0, -1, NA_real_, c(1, 2), "a", Inf * 0)) {
    expect_error(modify_eigen(floor = bad), "'floor' must be a single positive")
    expect_error(modify_ridge(floor = bad), "'floor' must be a single positive")
  }
  expect_error(hessian_modification(label = "x"), "abstract")
})


test_that("the objects print", {
  expect_output(print(modify_eigen()), "<hessian_modification> eigenvalue floor")
  expect_output(print(modify_ridge(floor = 1e-6)), "ridge \\(floor = 1e-06\\)")
  expect_output(print(newton(modify = modify_ridge())), "modify += ridge")
})


test_that("a string for modify is rejected with the object to use instead", {
  expect_error(newton(modify = "ridge"),
               "no longer accepted. Use modify_ridge\\(floor = \\.\\.\\.\\)")
  expect_error(newton(modify = "eigen"), "Use modify_eigen\\(floor = \\.\\.\\.\\)")
  expect_error(chebyshev(modify = "ridge"), "Use modify_ridge")
  expect_error(newton(modify = "other"), "modify_eigen\\(\\) or modify_ridge\\(\\)")
  expect_error(newton(modify = 1), "must be a Hessian modification object")
  expect_error(newton(modify = NULL), "must be a Hessian modification object")
  # the former arguments no longer exist
  expect_error(newton(hessian_mod = "ridge"), "unused argument")
  expect_error(newton(floor = 1e-6), "unused argument")
  expect_error(chebyshev(hessian_mod = "ridge"), "unused argument")
})


test_that("the default modification is the eigenvalue floor at 1e-8", {
  expect_s7_class(newton()@modify, optimizers7:::ModifyEigen)
  expect_equal(newton()@modify@floor, 1e-8)
  expect_s7_class(chebyshev()@modify, optimizers7:::ModifyEigen)
  expect_s7_class(chebyshev(modify = modify_ridge())@modify,
                  optimizers7:::ModifyRidge)
})


test_that("both shipped modifications reproduce the numbers of the string form", {
  # Values recorded with the code at the base of this change, in which the
  # modification was newton(hessian_mod = , floor = ). The comparison is
  # exact for the parameters up to 1e-12 and for the counts.
  ref <- function(opt, par, it, counts, start = c(0.05, 0.6)) {
    r <- run_saddle(opt, start)
    expect_equal(r@par, par, tolerance = 1e-12)
    expect_identical(r@iterations, as.integer(it))
    expect_equal(unname(unlist(r@counts[c("f", "g", "h")])), counts)
  }
  ref(newton(), c(1.0000000004227825, 0), 4, c(5, 5, 4))
  ref(newton(modify = modify_eigen()), c(1.0000000004227825, 0), 4, c(5, 5, 4))
  ref(newton(modify = modify_ridge()), c(1.0000000004227825, 0), 4, c(5, 5, 4))
  ref(chebyshev(modify = modify_eigen(floor = 1e-6)),
      c(1.0000000000072631, 0), 4, c(5, 11, 4), start = c(0.3, -2))
  ref(newton(modify = modify_ridge(floor = 1)),
      c(1.0000000134451885, 0), 5, c(6, 6, 5), start = c(0.3, -2))

  k <- test_problems("himmelblau")[[1]]
  r <- minimize(newton(), k$fn, k$par, gr = k$gr, he = k$he)
  expect_equal(r@par, c(3, 2), tolerance = 1e-12)
  expect_identical(r@iterations, 8L)
})


test_that("the floor changes the run", {
  a <- run_saddle(newton(modify = modify_eigen(floor = 1e-8), keep_trace = TRUE))
  b <- run_saddle(newton(modify = modify_eigen(floor = 1), keep_trace = TRUE))
  expect_false(identical(a@iterations, b@iterations))
  expect_identical(b@iterations, 8L)
})


test_that("the R methods of hessian_modify() match the compiled modifications", {
  # A subclass of a shipped class is called through R, so it runs the R method
  # of the parent; the run must agree with the compiled one.
  Eig2 <- S7::new_class("Eig2", parent = optimizers7:::ModifyEigen)
  Rdg2 <- S7::new_class("Rdg2", parent = optimizers7:::ModifyRidge)
  for (fl in c(1e-8, 1e-3, 1)) {
    for (pair in list(list(modify_eigen(floor = fl),
                           Eig2(label = "e2", floor = fl)),
                      list(modify_ridge(floor = fl),
                           Rdg2(label = "r2", floor = fl)))) {
      for (opt in c("newton", "chebyshev")) {
        make <- get(opt)
        a <- run_saddle(make(modify = pair[[1]], keep_trace = TRUE))
        b <- run_saddle(make(modify = pair[[2]], keep_trace = TRUE))
        expect_equal(b@par, a@par, tolerance = 1e-10)
        expect_identical(b@iterations, a@iterations)
        expect_identical(b@trace$safeguard, a@trace$safeguard)
      }
    }
  }
})


test_that("hessian_modify() returns a positive definite matrix", {
  H <- matrix(c(0.5, 1.5, 1.5, 0.5), 2, 2)
  expect_lt(min(eigen(H)$values), 0)
  for (m in list(modify_eigen(), modify_ridge(), modify_eigen(floor = 0.1))) {
    Hm <- hessian_modify(m, H)
    expect_equal(dim(Hm), dim(H))
    expect_true(isSymmetric(Hm))
    expect_gt(min(eigen(Hm)$values), 0)
  }
  # the eigenvalue floor leaves the eigenvalues above it unchanged; the
  # bound is floor times the largest absolute eigenvalue (2) when that is
  # above one
  Hm <- hessian_modify(modify_eigen(floor = 0.1), H)
  expect_equal(sort(eigen(Hm)$values), c(0.2, 2), tolerance = 1e-12)
  # the ridge shifts every eigenvalue by the same amount
  Hr <- hessian_modify(modify_ridge(floor = 1e-8), H)
  d <- eigen(Hr)$values - eigen(H)$values
  expect_equal(d[1], d[2], tolerance = 1e-12)
})


test_that("a user-defined modification is passed to newton() and chebyshev()", {
  ModifyAbs <- S7::new_class("ModifyAbs", parent = hessian_modification,
    properties = list(floor = S7::class_numeric))
  calls <- 0L
  S7::method(hessian_modify, ModifyAbs) <- function(modify, H) {
    calls <<- calls + 1L
    e <- eigen(H, symmetric = TRUE)
    e$vectors %*% (pmax(abs(e$values), modify@floor) * t(e$vectors))
  }
  m <- ModifyAbs(label = "absolute eigenvalues", floor = 1e-8)

  r <- run_saddle(newton(modify = m, keep_trace = TRUE))
  expect_true(r@converged)
  expect_equal(r@value, -1, tolerance = 1e-8)
  expect_gt(calls, 0L)
  expect_true(any(startsWith(r@trace$safeguard, "hessian modified")))
  expect_output(print(newton(modify = m)), "absolute eigenvalues")

  r2 <- run_saddle(chebyshev(modify = m, keep_trace = TRUE))
  expect_true(r2@converged)
  expect_equal(r2@value, -1, tolerance = 1e-8)

  # a Hessian that is positive definite never reaches the method
  calls <- 0L
  k <- test_problems("rosenbrock")[[1]]
  rr <- minimize(newton(modify = m), k$fn, k$par, gr = k$gr, he = k$he)
  expect_true(rr@converged)
  expect_identical(calls, 0L)
})


test_that("a modification that returns an unusable matrix falls back to the gradient", {
  Bad <- S7::new_class("Bad", parent = hessian_modification)
  S7::method(hessian_modify, Bad) <- function(modify, H) diag(1, 3)
  r <- run_saddle(newton(modify = Bad(label = "bad"), keep_trace = TRUE,
                         maxit = 5))
  expect_true("hessian solve failed" %in% r@trace$safeguard)

  # a matrix with the wrong sign is not a descent direction either
  Neg <- S7::new_class("Neg", parent = hessian_modification)
  S7::method(hessian_modify, Neg) <- function(modify, H) -diag(nrow(H))
  r2 <- run_saddle(newton(modify = Neg(label = "neg"), keep_trace = TRUE,
                          maxit = 5))
  expect_true("hessian solve failed" %in% r2@trace$safeguard)
})


test_that("a subclass without a method signals an error when it is used", {
  NoMethod <- S7::new_class("NoMethod", parent = hessian_modification)
  expect_error(run_saddle(newton(modify = NoMethod(label = "none"))),
               "hessian_modify")
})
