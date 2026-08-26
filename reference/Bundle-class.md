# S7 Class for the Proximal Bundle Method

An optimizer holding the proximity weight and its bounds, the acceptance
fraction that separates a serious step from a null one, the size of the
bundle of linearizations, and the effort spent on the subproblem. Built
by
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md).

## Usage

``` r
Bundle(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  t0 = integer(0),
  t_min = integer(0),
  t_max = integer(0),
  m_serious = integer(0),
  bundle_size = integer(0),
  qp_iters = integer(0),
  qp_tol = integer(0)
)
```

## Arguments

- t0:

  Initial proximity weight, as a step length on the parameter scale.

- t_min, t_max:

  Bounds on the proximity weight.

- m_serious:

  Fraction of the predicted decrease a serious step must achieve.

- bundle_size:

  Largest number of linearizations kept.

- qp_iters, qp_tol:

  Effort spent on the subproblem.

## Value

An S7 object of class `Bundle` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the seven properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `Bundle` carries
seven of its own, in three groups: the proximity weight (`t0`, `t_min`,
`t_max`), the acceptance test (`m_serious`), and the model and its
subproblem (`bundle_size`, `qp_iters`, `qp_tol`).

## See also

[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
for the constructor,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the rule it reads.
