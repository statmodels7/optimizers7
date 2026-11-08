test_that("descent_settings() returns the arguments minimize() passes to the loop", {
  q <- function(p) sum((p - c(1, 2))^2)
  qg <- function(p) 2 * (p - c(1, 2))
  qh <- function(p) diag(2, 2)
  opts <- list(newton(), newton(max_length = 5), bfgs(), lbfgs(memory = 4L),
               gd(), cg(), bb(),
               lbfgs(criterion = crit_any(crit_grad(1e-7), crit_abs_obj(1e-9)),
                     line_search = wolfe(max_step = 12)),
               newton(typical_size = c(1, 3)))
  orig <- descent_run
  seen <- NULL
  local_mocked_bindings(descent_run = function(...) {
    seen <<- list(...)
    orig(...)
  })
  for (o in opts) {
    seen <- NULL
    minimize(o, q, c(0, 0), gr = qg, hess = qh)
    s <- descent_settings(o, 2)
    for (k in c("method", "line_search", "maxit", "max_eval", "verbose",
                "refresh", "keep_trace", "step", "rule_spec")) {
      expect_identical(s[[k]], seen[[k]], info = paste(o@name, k))
    }
  }
})

test_that("descent_settings() is NULL where the loop does not run", {
  expect_null(descent_settings(nelder_mead(), 2))
  expect_null(descent_settings(chebyshev(), 2))
})
