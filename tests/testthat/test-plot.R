rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                          200 * (p[2] - p[1]^2))

test_that("the graph of a result draws the trace and names its elements", {
  skip_if_not_installed("ggplot2")
  res <- minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1), gr = rosen_gr)
  pdf(file = NULL)
  on.exit(dev.off())
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, res@optimizer@name)
  b <- ggplot2::ggplot_build(p)
  # the line is the objective at every iteration of the trace
  expect_equal(b$data[[1]]$y, res@trace$value)
  # the points are the iterations at which a safeguard fired
  hit <- res@trace$safeguard != "none"
  expect_true(any(hit))
  expect_equal(b$data[[2]]$x, res@trace$iteration[hit])
  expect_identical(b$plot$scales$get_scales("colour")$get_labels(),
                   c("objective", "iteration at which a safeguard fired"))
  expect_identical(p$theme$legend.position, "bottom")
  expect_identical(plot(res, legend = FALSE)$theme$legend.position, "none")
  expect_s3_class(ggplot2::autoplot(res), "ggplot")
})

test_that("a run without a safeguard has no entry for it", {
  skip_if_not_installed("ggplot2")
  res <- minimize(newton(keep_trace = TRUE), function(p) sum((p - 1:2)^2),
                  c(0, 0))
  p <- ggplot2::autoplot(res)
  expect_length(p$layers, 1L)
  b <- ggplot2::ggplot_build(p)
  expect_identical(b$plot$scales$get_scales("colour")$get_labels(), "objective")
})

test_that("the graph of a result rejects what it cannot draw", {
  skip_if_not_installed("ggplot2")
  res <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr)
  expect_error(plot(res), "keep_trace = TRUE")
  res <- minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1), gr = rosen_gr)
  expect_error(plot(res, log = "y"), "no further arguments")
  expect_error(plot(res, legend = NA), "'legend' must be TRUE or FALSE")
})
