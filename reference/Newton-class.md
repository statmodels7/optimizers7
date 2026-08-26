# S7 Class for Newton's Method

An optimizer holding the initial step length, the line search, and the
two settings that decide how an indefinite Hessian is repaired. Built by
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).
It is the only shipped method that reads the `he` argument of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Usage

``` r
Newton(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  hessian_mod = character(0),
  floor = integer(0)
)
```

## Arguments

- step, line_search:

  The initial step length and the
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object, as in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- hessian_mod:

  How an indefinite Hessian is repaired, `"eigen"` or `"ridge"`.

- floor:

  The smallest eigenvalue the repaired Hessian may have.

## Value

An S7 object of class `Newton` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the four properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `Newton` carries
four: `step` and `line_search`, shared with the other line-search
methods, and `hessian_mod` and `floor`, which are its own.

## See also

[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
for the constructor,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the method that needs no Hessian.
