# optimizers7: An S7 Framework for Optimisation Algorithms

Optimisation algorithms as objects, built on the S7 object system with
the numerical work in C++. An optimiser carries its own settings,
safeguards and stopping rule, so swapping one for another changes a
single word; the stopping rule is itself a composable object, so what a
run means by convergence is something the caller sets rather than
something fixed when the package was written. Eleven methods are
provided, spanning second-order, first-order, derivative-free and
non-smooth problems, along with multi-start. Box constraints are removed
by reparameterisation instead of enforced, so every method accepts them
without knowing they exist, and every method reports which of its
safeguards fired.

## See also

Useful links:

- <https://statmodels7.github.io/optimizers7/>

- <https://github.com/statmodels7/optimizers7>

- Report bugs at <https://github.com/statmodels7/optimizers7/issues>

## Author

**Maintainer**: Giovanni Tinervia <giovannitinervia9@gmail.com>
