# optimizers7: An S7 Framework for Optimization Algorithms

Optimization algorithms as objects, built on the S7 object system with
the numerical work in C++. An optimizer carries its own settings, its
safeguards and its stopping rule, so swapping one for another changes a
single word. The stopping rule is a composable object of its own, so the
caller decides what a run means by convergence. Twelve methods cover
second-order, first-order, derivative-free, non-smooth and global
search, and two wrappers compose them: multi-start and a chain of
stages. Box constraints are removed by reparametrization, so every
method accepts them without knowing they exist, and every method reports
which of its safeguards fired.

## See also

Useful links:

- <https://statmodels7.github.io/optimizers7/>

- <https://github.com/statmodels7/optimizers7>

- Report bugs at <https://github.com/statmodels7/optimizers7/issues>

## Author

**Maintainer**: Giovanni Tinervia <giovannitinervia9@gmail.com>
