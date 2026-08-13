# What Simulated Annealing Can Offer a Stopping Rule

The objective and a stationarity measure, but no gradient.

## Arguments

- optimizer:

  An `Sa` object.

## Value

A character vector.

## Details

The measure is Corana's termination rule, by how much the best value has
moved over the last `n_eps` temperature levels. A rule reading a
gradient is rejected at construction rather than left testing `NULL` and
never firing.
