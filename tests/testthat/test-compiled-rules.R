# The compiled stopping rules against the R ones: every shipped rule, on the
# smooth test problems, with every descent direction, must end the run at the
# same iteration and point as crit_met() in R.

run_both <- function(opt_fn, q, par, ...) {
  old <- options(optimizers7.compiled_rules = FALSE,
                 optimizers7.check_gradient = FALSE)
  on.exit(options(old))
  a <- minimize(opt_fn(), q$fn, par, gr = q$gr, hess = q$hess, ...)
  options(optimizers7.compiled_rules = TRUE)
  b <- minimize(opt_fn(), q$fn, par, gr = q$gr, hess = q$hess, ...)
  list(a = a, b = b)
}

same_run <- function(a, b) {
  identical(a@par, b@par) && identical(a@value, b@value) &&
    identical(a@gradient, b@gradient) && identical(a@counts, b@counts) &&
    identical(a@iterations, b@iterations) &&
    identical(a@converged, b@converged) &&
    identical(a@criterion_met, b@criterion_met) &&
    identical(a@message, b@message) && identical(a@trace, b@trace)
}

test_that("compiled rules stop every descent method where the R rules do", {
  P <- test_problems()
  P <- P[vapply(P, function(q) isTRUE(q$smooth), logical(1))]
  opts <- list(
    function() newton(),
    function() newton(max_length = 5, line_search = armijo(max_step = 12)),
    function() bfgs(criterion = crit_any(crit_grad(1e-7, norm = "2"),
                                         crit_rel_par(1e-10))),
    function() lbfgs(),
    function() lbfgs(criterion = crit_any(crit_rel_obj(1e-10), crit_rel_par(1e-8))),
    function() lbfgs(criterion = crit_all(crit_grad(1e-4), crit_abs_obj(1e-6))),
    function() lbfgs(criterion = crit_never(), maxit = 40),
    function() cg(),
    function() bb(),
    function() gd(maxit = 200)
  )
  n <- 0L
  for (q in P) for (o in opts) for (box in c(FALSE, TRUE)) {
    r <- if (box) {
      run_both(o, q, q$par, lower = q$par - 3, upper = q$par + 3)
    } else {
      run_both(o, q, q$par)
    }
    expect_true(same_run(r$a, r$b), info = paste(q$name, r$a@optimizer@name, box))
    n <- n + 1L
  }
  expect_equal(n, length(P) * length(opts) * 2L)
})

test_that("criterion_spec() describes the shipped rules and only those", {
  s <- criterion_spec(crit_any(crit_grad(1e-7, norm = "2"), crit_abs_obj(1e-8),
                               crit_all(crit_rel_par(1e-6), crit_never())))
  expect_identical(s$type, "any")
  expect_identical(s$criteria[[1L]], list(type = "grad", tol = 1e-7, norm = "2"))
  expect_identical(s$criteria[[2L]], list(type = "abs_obj", tol = 1e-8))
  expect_identical(s$criteria[[3L]]$type, "all")
  expect_identical(s$criteria[[3L]]$criteria[[2L]], list(type = "never"))
  Tiny <- S7::new_class("Tiny", parent = criterion,
                        properties = list(tol = S7::class_numeric))
  S7::method(crit_met, Tiny) <- function(criterion, state) abs(state$f_new) < criterion@tol
  expect_null(criterion_spec(Tiny(label = "tiny", tol = 1e-3)))
  expect_null(criterion_spec(crit_any(crit_grad(), Tiny(label = "tiny", tol = 1e-3))))
})

test_that("a rule without a compiled form is asked in R", {
  Tiny <- S7::new_class("Tiny2", parent = criterion,
                        properties = list(tol = S7::class_numeric))
  S7::method(crit_met, Tiny) <- function(criterion, state) abs(state$f_new) < criterion@tol
  q <- test_problems("sphere")[[1L]]
  old <- options(optimizers7.compiled_rules = TRUE)
  on.exit(options(old))
  r <- minimize(lbfgs(criterion = Tiny(label = "tiny", tol = 1e-12)), q$fn,
                q$par, gr = q$gr)
  expect_true(r@converged)
  expect_identical(r@criterion_met, "tiny")
  expect_lt(abs(r@value), 1e-12)
})
