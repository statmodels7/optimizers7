#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

# Starting values that the caller does not have to write down.
#
# Every optimizer needs a point to start from, and for a model with a dozen
# parameters writing that point out is busywork of the least interesting kind --
# the more so because the natural thing to write, a vector of zeros, is wrong on
# the parameter scale as soon as a variance is involved. A starter is an object
# that stands in for the vector: it says how the values are to be produced, and
# it is turned into an actual vector at the top of minimize(), on the
# unconstrained scale, and mapped back through the box.
#
# Resolution happens in the BODY of the minimize() generic, before dispatch, in
# the same way distributions7 intercepts its link scale. Every optimizer
# therefore accepts a starter, user-written ones included, and no method has to
# know that starters exist.


#' @title S7 Class for a Starting-Value Generator
#'
#' @description
#' The abstract parent of [start_zeros()] and [start_runif()]. A starter
#' stands in for the vector of starting values: it says how the values are to
#' be produced and carries how many of them there are, and [minimize()]
#' turns it into an actual vector before dispatch. Every optimizer therefore
#' accepts one, including an optimizer written outside the package, and no
#' method needs to know that starters exist.
#'
#' @details
#' The class is abstract and carries one property, `npar`, an integer or
#' `NULL`. A subclass needs a method for [starting_values()] and nothing
#' else. The class is not exported as of version 0.6.0, so subclassing it is
#' available inside the package and not outside it; [is_starter()] is what
#' [minimize()] tests, and a class parented elsewhere fails that test.
#'
#' @param npar The number of parameters, an integer, or `NULL` to have
#'   [minimize()] work it out from the bounds or from the objective.
#'
#' @return An S7 object. The class is abstract, so every value is an object of
#'   one of its subclasses.
#'
#' @seealso [start_zeros()], [start_runif()], [starting_values()].
#' @name starter-class
#' @aliases starter
#' @keywords internal
starter <- S7::new_class("starter", abstract = TRUE,
  properties = list(npar = S7::class_any))


#' @title Is This a Starter?
#'
#' @description
#' `TRUE` when `x` inherits from the abstract [starter] class, `FALSE`
#' otherwise. [minimize()] asks it to decide whether `par` is a vector to be
#' used as given or an object to be resolved into one, and the question has a
#' name so that the test is written once.
#'
#' @param x Any object.
#'
#' @return A single logical.
#'
#' @examples
#' is_starter(start_zeros())
#' is_starter(start_runif(-2, 2))
#' is_starter(c(0, 0))
#' is_starter(bfgs())
#'
#' @keywords internal
is_starter <- function(x) S7::S7_inherits(x, starter)


#' @title Produce a Vector of Starting Values
#'
#' @description
#' Turns a starter into an actual numeric vector of length `npar`, on the
#' **unconstrained** scale.
#'
#' @param starter A [start_zeros()] or [start_runif()]
#'   object, or a user-defined starter.
#' @param npar The number of parameters wanted.
#'
#' @details
#' The unconstrained scale is the one the optimizer works on when there are
#' bounds, and that is where a starter is entitled to be simple: zero means
#' the middle of an interval, one for a variance, one half for a probability,
#' and no value can fall outside a bound. [minimize()] maps the result back
#' through [bounded_transform()] before any method sees it.
#'
#' `npar` is passed by [minimize()], which has already settled it from the
#' starter's own `npar`, from the length of the bounds, or by probing the
#' objective with [infer_npar()]. A method for this generic therefore reads
#' the argument and not `starter@npar`.
#'
#' The two shipped starters are [start_zeros()] and [start_runif()], and they
#' are what a caller has available. A starter of a third kind is a subclass of
#' the abstract [starter] class with a method for this generic and nothing
#' else, but that class is not exported as of version 0.6.0, and
#' [minimize()] accepts as `par` only a numeric vector or an object
#' inheriting from it. A class of your own carrying a `starting_values()`
#' method is therefore refused with `'par' must be a numeric vector of
#' starting values`; resolve the vector yourself and pass it.
#'
#' @return A numeric vector of length `npar`, on the unconstrained scale.
#'
#' @examples
#' starting_values(start_zeros(), 3)
#'
#' set.seed(1)
#' starting_values(start_runif(-2, 2), 3)
#'
#' # Resolving by hand gives exactly the vector minimize() would have used,
#' # so passing a starter and passing its values are the same run.
#' f <- function(p) sum((p - 1:3)^2)
#' set.seed(1); a <- minimize(bfgs(), f, start_runif(-2, 2, npar = 3))
#' set.seed(1); b <- minimize(bfgs(), f, starting_values(start_runif(-2, 2), 3))
#' identical(a@par, b@par)
#'
#' @seealso [start_zeros()] and [start_runif()] for the two shipped starters,
#'   [minimize()] for where the resolution happens, [infer_npar()] for how
#'   `npar` is worked out when nothing declares it.
#' @export
starting_values <- S7::new_generic("starting_values", "starter",
  function(starter, npar) S7::S7_dispatch())


# --- zeros ------------------------------------------------------------------

#' @title S7 Class for the Zero Starter
#'
#' @description
#' A starter that produces a vector of zeros on the unconstrained scale.
#' Built by [start_zeros()]. It carries `npar` alone, adding no property of
#' its own to the abstract [starter] class.
#'
#' @param npar The number of parameters, an integer, or `NULL`.
#'
#' @return An S7 object of class `ZeroStart` inheriting from [starter].
#'
#' @seealso [start_zeros()] for the constructor,
#'   [starting_values.ZeroStart()] for what it produces.
#' @name ZeroStart-class
#' @aliases ZeroStart
#' @keywords internal
ZeroStart <- S7::new_class("ZeroStart", parent = starter)


#' @title Start From Zero on the Unconstrained Scale
#'
#' @description
#' A starting point of all zeros, which after the bound transform is the
#' sensible middle of every parameter's domain: one for a positive parameter,
#' one half for a probability, zero for an unbounded one.
#'
#' @param npar The number of parameters. Defaults to `NULL`, meaning work
#'   it out from the bounds or from the objective; see [minimize()].
#'
#' @details
#' Zero is the right constant only because it is applied on the unconstrained
#' scale. What it becomes after the bound transform depends on the box:
#'
#' \tabular{ll}{
#'   **bounds** \tab **starting value** \cr
#'   \eqn{(-\infty, \infty)} \tab 0 \cr
#'   \eqn{(0, \infty)} \tab 1 \cr
#'   \eqn{(0, 1)} \tab 0.5 \cr
#'   \eqn{(2, 6)} \tab 4 \cr
#'   \eqn{(-\infty, 5)} \tab 4 \cr
#' }
#'
#' A vector of zeros on the *parameter* scale is no starting point at all for
#' a model with a scale parameter in it: it sits exactly on the boundary,
#' where the log-likelihood is usually infinite and the gradient certainly is.
#' [minimize()] refuses such a point for that reason.
#'
#' @return An S7 object of class [ZeroStart], inheriting from [starter], to be
#'   passed as `par`.
#'
#' @examples
#' f <- function(p) sum((p - c(1, 2, 3))^2)
#' minimize(bfgs(), f, start_zeros(3))@par
#'
#' # Zero on the unconstrained scale is one on a positive parameter's scale,
#' # and the midpoint of a two-sided box.
#' starting_values(start_zeros(), 3)
#' bounded_transform(c(0, Inf), 0)$h
#' bounded_transform(c(2, 6), 0)$h
#'
#' # The count need not be given when the bounds already say it.
#' minimize(bfgs(), f, start_zeros(), lower = c(0, 0, 0))@par
#'
#' @seealso [start_runif()] for a random start, [starting_values()] for the
#'   generic that resolves it, [minimize()] for how `npar` is settled.
#' @export
start_zeros <- function(npar = NULL) ZeroStart(npar = check_npar(npar))


#' @title Zero Starting Values
#' @name starting_values.ZeroStart
#'
#' @description
#' Returns `npar` zeros. Read on the unconstrained scale, so [minimize()]
#' maps them through the bounds before any method sees them: a positive
#' parameter starts at 1, a probability at 0.5, an unbounded one at 0.
#'
#' @param starter A `ZeroStart` object. Its own `npar` is not read here;
#'   [minimize()] has already settled the count and passes it.
#' @param npar The number of parameters wanted, a positive whole number.
#'
#' @return A numeric vector of `npar` zeros.
#'
#' @examples
#' starting_values(start_zeros(), 4)
#'
#' @keywords internal
S7::method(starting_values, ZeroStart) <- function(starter, npar) {
  numeric(npar)
}


# --- uniform ----------------------------------------------------------------

#' @title S7 Class for the Uniform Starter
#'
#' @description
#' A starter that draws each coordinate independently from a uniform on the
#' unconstrained scale. Built by [start_runif()]. It adds `min` and `max` to
#' the `npar` the abstract [starter] class carries; either may be one number
#' or one per parameter.
#'
#' @param min,max The range drawn from, in unconstrained units.
#' @param npar The number of parameters, an integer, or `NULL`.
#'
#' @return An S7 object of class `UniformStart` inheriting from [starter].
#'
#' @seealso [start_runif()] for the constructor,
#'   [starting_values.UniformStart()] for the draw.
#' @name UniformStart-class
#' @aliases UniformStart
#' @keywords internal
UniformStart <- S7::new_class("UniformStart", parent = starter,
  properties = list(min = S7::class_numeric, max = S7::class_numeric))


#' @title Start From a Uniform Draw on the Unconstrained Scale
#'
#' @description
#' Each coordinate is drawn independently from `runif(min, max)` on the
#' unconstrained scale and mapped back through the bounds, so no draw is ever
#' rejected for being outside the box.
#'
#' @param min,max The range to draw from, in unconstrained units. Both default
#'   to a width of one either side of zero, and both may be given per parameter
#'   rather than as a single number.
#' @param npar The number of parameters. Defaults to `NULL`, meaning work
#'   it out from the bounds or from the objective; see [minimize()].
#'
#' @details
#' The range is in unconstrained units, and that is what makes a single
#' default workable. A draw in \eqn{(-1, 1)} becomes a variance between
#' `0.368` and `2.72`, a probability between `0.269` and `0.731`, and a
#' parameter bounded on both sides lands well inside its interval. The same
#' numbers on the parameter scale would mean quite different things and would
#' sometimes be inadmissible.
#'
#' Widen it when the scale of the problem is unknown. `start_runif(-5, 5)`
#' spans `0.0067` to `148` for a positive parameter, four orders of
#' magnitude, and is still a range no draw can fall out of.
#'
#' The draw uses \R's ordinary generator, so [set.seed()] reproduces it, and
#' the state is recorded in the result's `seed`.
#'
#' @return An S7 object of class [UniformStart], inheriting from [starter], to
#'   be passed as `par`.
#'
#' @examples
#' f <- function(p) sum((p - c(1, 2, 3))^2)
#' set.seed(1)
#' minimize(bfgs(), f, start_runif(npar = 3))@par
#'
#' # What the range means on the parameter scale, for a positive parameter.
#' range(bounded_transform(c(0, Inf), c(-1, 1))$h)
#' range(bounded_transform(c(0, Inf), c(-5, 5))$h)
#'
#' # A wider net on a positive parameter whose scale is unknown.
#' set.seed(1)
#' minimize(bfgs(), function(p) (log(p) - 1)^2, start_runif(-5, 5, npar = 1),
#'          lower = 0)@par
#'
#' # One range per parameter is allowed; a length that is neither 1 nor npar
#' # is refused when the draw is made.
#' set.seed(1)
#' starting_values(start_runif(c(-1, -10), c(1, 10)), 2)
#' try(starting_values(start_runif(c(-1, -2, -3)), 2))
#'
#' @seealso [start_zeros()] for the deterministic start, [starting_values()]
#'   for the generic that resolves it, [multistart()], which uses a starter to
#'   generate its own starts.
#' @export
start_runif <- function(min = -1, max = 1, npar = NULL) {
  if (!is.numeric(min) || !length(min) || anyNA(min) || !all(is.finite(min))) {
    stop("'min' must be finite and numeric.", call. = FALSE)
  }
  if (!is.numeric(max) || !length(max) || anyNA(max) || !all(is.finite(max))) {
    stop("'max' must be finite and numeric.", call. = FALSE)
  }
  if (any(min >= max)) {
    stop("'min' must be strictly below 'max'.", call. = FALSE)
  }
  UniformStart(npar = check_npar(npar), min = as.numeric(min),
               max = as.numeric(max))
}


#' @title Uniform Starting Values
#' @name starting_values.UniformStart
#'
#' @description
#' Draws `npar` values, coordinate by coordinate, from the starter's `min`
#' and `max` on the unconstrained scale. A `min` or `max` of length one is
#' used for every coordinate; one of length `npar` gives each its own range,
#' and any other length raises an error naming both lengths.
#'
#' @param starter A `UniformStart` object, read for `min` and `max`.
#' @param npar The number of parameters wanted, a positive whole number.
#'
#' @return A numeric vector of length `npar`, drawn with `stats::runif()`, so
#'   `set.seed()` governs it.
#'
#' @examples
#' set.seed(1)
#' starting_values(start_runif(-2, 2), 4)
#'
#' # A range per coordinate.
#' set.seed(1)
#' starting_values(start_runif(c(-1, -10), c(1, 10)), 2)
#'
#' # A length that is neither 1 nor npar is a mistake, not a request.
#' try(starting_values(start_runif(c(-1, -2, -3)), 2))
#'
#' @keywords internal
S7::method(starting_values, UniformStart) <- function(starter, npar) {
  lo <- recycle_to(starter@min, npar, "min")
  up <- recycle_to(starter@max, npar, "max")
  stats::runif(npar, lo, up)
}


# --- resolution -------------------------------------------------------------

#' Validate a Declared Parameter Count
#'
#' @description
#' Checks that `npar` is a single positive whole number, or `NULL`, and
#' returns it as an integer. Called by both starter constructors, so
#' `start_zeros(0)` and `start_runif(npar = 2.5)` are refused in the same
#' words at the point they are written.
#'
#' @param npar `NULL` or a positive whole number.
#'
#' @return `NULL` when `npar` is `NULL`, otherwise the value as an integer.
#'   Raises an error otherwise.
#'
#' @keywords internal
check_npar <- function(npar) {
  if (is.null(npar)) return(NULL)
  if (!is.numeric(npar) || length(npar) != 1L || is.na(npar) ||
      npar < 1 || npar != round(npar)) {
    stop("'npar' must be a single positive whole number, or NULL.",
         call. = FALSE)
  }
  as.integer(npar)
}


#' Recycle a Length-One Vector, and Reject Any Other Mismatch
#'
#' @description
#' Returns `v` at length `n`: a single value is repeated, a value already of
#' length `n` is passed through as a double, and anything else raises an
#' error naming the argument and both lengths. \R's own recycling is
#' deliberately not used, since it is silent whenever the shorter length
#' divides the longer, and a partial range is far likelier to be a mistake
#' than a request.
#'
#' @param v A numeric vector.
#' @param n The length wanted.
#' @param nm The argument's name, for the message.
#'
#' @return A numeric vector of length `n`.
#'
#' @keywords internal
recycle_to <- function(v, n, nm) {
  if (length(v) == 1L) return(rep(as.numeric(v), n))
  if (length(v) != n) {
    stop("'", nm, "' must have length 1 or ", n, ", one per parameter; it has ",
         "length ", length(v), ".", call. = FALSE)
  }
  as.numeric(v)
}


#' @title How Many Parameters the Objective Takes
#'
#' @description
#' Works out the length of the parameter vector by trying lengths and seeing
#' which the objective accepts. Called by [minimize()] when a starter
#' was given without `npar` and the bounds do not say.
#'
#' @param fn The objective.
#' @param gr Its gradient, or `NULL`. Supplying one makes the answer much
#'   more likely to be unique; see Details.
#' @param probe A function of one integer returning a candidate parameter vector
#'   of that length, so that the objective is probed where it will be used.
#' @param npar_max The largest length tried. Defaults to `50`.
#'
#' @details
#' A length is accepted when `fn` returns a single finite number for it and
#' raises neither an error nor a warning, and, if `gr` was supplied, when
#' the gradient comes back with the same length as its argument.
#'
#' What decides it is whether the objective genuinely rejects the wrong
#' length, and the objective a modeling package hands over usually does:
#' `X %*% beta` with a parameter of the wrong length is an error rather
#' than a number, so a regression of any kind is settled at once. A gradient
#' helps when it spells its components out, since such a gradient returns a fixed
#' number of them whatever it is handed.
#'
#' A vectorized objective written in terms of the parameter alone is another
#' matter: both plausible
#' guesses about \R are wrong. Recycling warns only when the shorter length is
#' not a *divisor* of the longer, so `sum((p - c(1, 2, 3))^2)` accepts
#' a length-one vector in silence and returns a perfectly finite 14; and its
#' gradient `2 * (p - c(1, 2, 3))` returns six components for a length-six
#' argument, so the gradient rule passes it too. Rosenbrock's
#' `100 (p_2 - p_1^2)^2 + (1 - p_1)^2` accepts every length from two upwards
#' for the same kind of reason. None of this is a defect and none of it can be
#' guessed at.
#'
#' The probe therefore settles the objectives that have a
#' fixed width built into them and rejects the ones that do not, naming the two
#' lengths it found. When it rejects, `npar` or a vector of bounds is one
#' word.
#'
#' The search stops as soon as a *second* length is accepted: the answer is
#' then already known to be ambiguous and there is no reason to keep probing.
#' The cost is therefore two evaluations when the objective accepts any
#' length, and `npar_max` when it accepts exactly one. Either way it happens
#' once, before the run.
#'
#' @return A single integer, the one length accepted. Raises an error naming
#'   the two lengths it found when more than one is accepted, and an error
#'   when none is.
#'
#' @examples
#' # A hand-written gradient pins it exactly.
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' infer_npar(f, gr, function(k) numeric(k))
#'
#' # Without one, the same objective is happy with any length from two
#' # upwards, and the refusal names the two lengths that decided it.
#' try(infer_npar(f, NULL, function(k) numeric(k)))
#'
#' # An objective with a width built in is settled without a gradient, and
#' # this is the ordinary case for a model.
#' set.seed(1)
#' X <- matrix(rnorm(40 * 3), 40, 3)
#' y <- as.numeric(X %*% c(1, -2, 0.5) + rnorm(40))
#' infer_npar(function(b) sum((y - X %*% b)^2), NULL, function(k) numeric(k))
#'
#' @seealso [start_zeros()] and [start_runif()], the starters that make this
#'   question arise, and [minimize()], which asks it.
#' @export
infer_npar <- function(fn, gr, probe, npar_max = 50) {
  accepts <- function(k) {
    x <- probe(k)
    bad <- FALSE
    v <- withCallingHandlers(
      tryCatch(fn(x), error = function(e) NULL),
      warning = function(w) {
        bad <<- TRUE
        invokeRestart("muffleWarning")
      })
    if (bad || !is.numeric(v) || length(v) != 1L || !is.finite(v)) return(FALSE)
    if (is.null(gr)) return(TRUE)
    bad <- FALSE
    g <- withCallingHandlers(
      tryCatch(gr(x), error = function(e) NULL),
      warning = function(w) {
        bad <<- TRUE
        invokeRestart("muffleWarning")
      })
    !bad && is.numeric(g) && length(g) == k && all(is.finite(g))
  }

  found <- integer(0)
  for (k in seq_len(npar_max)) {
    if (accepts(k)) {
      found <- c(found, k)
      if (length(found) == 2L) break
    }
  }

  if (!length(found)) {
    stop("Could not work out how many parameters the objective takes: no ",
         "length\n  from 1 to ", npar_max, " gave a single finite value ",
         "without a warning.\n  Say how many, as in start_zeros(npar = 3), or ",
         "pass a numeric starting value.", call. = FALSE)
  }
  if (length(found) > 1L) {
    stop("The objective accepts more than one length of parameter vector (",
         found[1], " and ", found[2], "),\n  so the number of parameters ",
         "cannot be worked out from it. Say how many,\n  as in start_zeros(",
         "npar = ", found[2], ").", call. = FALSE)
  }
  found
}


#' Turn a Starter Into a Starting Value
#'
#' @description
#' The whole of what [minimize()] does with a starter, in one place:
#' settle the number of parameters, draw the values on the unconstrained scale,
#' and map them back through the box.
#'
#' @details
#' The number of parameters is looked for in three places, in order of how much
#' the caller was willing to say. `npar` on the starter itself is taken as
#' given. Failing that, a `lower` or `upper` of length greater than one
#' answers the question, since bounds are one per parameter. Failing both,
#' [infer_npar()] probes the objective.
#'
#' A numeric `par` passes through untouched, so this costs nothing at all
#' for the ordinary call.
#'
#' @param par Whatever was passed as `par`.
#' @param fn,gr The objective and its gradient.
#' @param lower,upper The bounds.
#'
#' @return A numeric vector on the parameter scale.
#'
#' @keywords internal
resolve_start <- function(par, fn, gr, lower, upper) {
  if (!is_starter(par)) return(par)
  if (!is.function(fn)) {
    stop("A starter needs an objective it can probe; 'fn' is not a function.",
         call. = FALSE)
  }

  p <- par@npar
  if (is.null(p)) {
    lens <- c(length(lower), length(upper))
    if (max(lens) > 1L) {
      p <- max(lens)
      if (all(lens > 1L) && lens[1] != lens[2]) {
        stop("'lower' and 'upper' have lengths ", lens[1], " and ", lens[2],
             ", so they disagree\n  about how many parameters there are.",
             call. = FALSE)
      }
    }
  }
  if (is.null(p)) {
    p <- infer_npar(fn, gr, function(k) starting_values(par, k))
  }

  eta <- starting_values(par, p)
  if (!is.numeric(eta) || length(eta) != p || anyNA(eta)) {
    stop("starting_values() must return ", p, " numbers, and not NA.",
         call. = FALSE)
  }

  # Onto the parameter scale, coordinate by coordinate: the starter works in
  # eta, the objective is a function of theta.
  lo <- recycle_to(if (is.null(lower)) -Inf else lower, p, "lower")
  up <- recycle_to(if (is.null(upper)) Inf else upper, p, "upper")
  if (all(!is.finite(lo)) && all(!is.finite(up))) return(eta)

  theta <- eta
  for (j in seq_len(p)) {
    if (lo[j] >= up[j]) {
      stop("For parameter ", j, " the lower bound must be strictly below the ",
           "upper one; they are ", format(lo[j]), " and ", format(up[j]), ".",
           call. = FALSE)
    }
    theta[j] <- bounded_transform(c(lo[j], up[j]), eta[j])$h
  }
  theta
}
