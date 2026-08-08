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
#' @description The abstract parent of \code{\link{start_zeros}} and
#'   \code{\link{start_runif}}.
#' @param npar The number of parameters, or \code{NULL} to work it out.
#' @return An S7 object.
#' @seealso \code{\link{start_zeros}}, \code{\link{start_runif}}
#' @name starter-class
#' @aliases starter
#' @keywords internal
starter <- S7::new_class("starter", abstract = TRUE,
  properties = list(npar = S7::class_any))


#' @title The Number of Parameters a Starter Was Given
#' @description The class of an optimizer is not the place to look for this, so
#'   it has a name.
#' @param x Any object.
#' @return \code{TRUE} for a starter.
#' @keywords internal
is_starter <- function(x) S7::S7_inherits(x, starter)


#' @title Produce a Vector of Starting Values
#'
#' @description
#' Turns a starter into an actual numeric vector of length \code{npar}, on the
#' \strong{unconstrained} scale.
#'
#' @param starter A \code{\link{start_zeros}} or \code{\link{start_runif}}
#'   object, or a user-defined starter.
#' @param npar The number of parameters wanted.
#'
#' @details
#' The unconstrained scale is the one the optimizer actually works on when there
#' are bounds, so this is where a starter is entitled to be simple: zero means
#' the middle of an interval, one for a variance, one half for a probability,
#' and there is no way for it to fall outside a bound. \code{\link{minimize}}
#' maps the result back through \code{\link{bounded_transform}} before any
#' method sees it.
#'
#' A user-defined starter is a subclass of \code{starter} with a method for
#' this generic; nothing else is required.
#'
#' @return A numeric vector of length \code{npar}.
#'
#' @examples
#' starting_values(start_zeros(), 3)
#'
#' set.seed(1)
#' starting_values(start_runif(-2, 2), 3)
#'
#' @seealso \code{\link{start_zeros}}, \code{\link{start_runif}},
#'   \code{\link{minimize}}
#' @export
starting_values <- S7::new_generic("starting_values", "starter",
  function(starter, npar) S7::S7_dispatch())


# --- zeros ------------------------------------------------------------------

#' @title S7 Class for the Zero Starter
#' @description The class \code{\link{start_zeros}} instantiates.
#' @param npar The number of parameters, or \code{NULL}.
#' @return An S7 object inheriting from \code{starter}.
#' @seealso \code{\link{start_zeros}}
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
#' @param npar The number of parameters. Defaults to \code{NULL}, meaning work
#'   it out from the bounds or from the objective; see \code{\link{minimize}}.
#'
#' @details
#' Zero is the right constant only because it is applied on the unconstrained
#' scale. A vector of zeros on the \emph{parameter} scale is not a starting
#' point at all for a model with a scale parameter in it: it sits exactly on the
#' boundary, where the log-likelihood is usually infinite and the gradient
#' certainly is.
#'
#' @return An S7 object of class \code{ZeroStart}.
#'
#' @examples
#' f <- function(p) sum((p - c(1, 2, 3))^2)
#' minimize(bfgs(), f, start_zeros(3))@par
#'
#' # zero on the unconstrained scale is one on a positive parameter's scale
#' minimize(bfgs(), f, start_zeros(3), lower = 0)@par
#'
#' @seealso \code{\link{start_runif}}, \code{\link{starting_values}}
#' @export
start_zeros <- function(npar = NULL) ZeroStart(npar = check_npar(npar))


#' @title Zero Starting Values
#' @name starting_values.ZeroStart
#' @description Returns \code{npar} zeros.
#' @param starter A \code{ZeroStart} object.
#' @param npar The number of parameters.
#' @return A numeric vector of zeros.
#' @keywords internal
S7::method(starting_values, ZeroStart) <- function(starter, npar) {
  numeric(npar)
}


# --- uniform ----------------------------------------------------------------

#' @title S7 Class for the Uniform Starter
#' @description The class \code{\link{start_runif}} instantiates.
#' @param min,max The range drawn from.
#' @param npar The number of parameters, or \code{NULL}.
#' @return An S7 object inheriting from \code{starter}.
#' @seealso \code{\link{start_runif}}
#' @name UniformStart-class
#' @aliases UniformStart
#' @keywords internal
UniformStart <- S7::new_class("UniformStart", parent = starter,
  properties = list(min = S7::class_numeric, max = S7::class_numeric))


#' @title Start From a Uniform Draw on the Unconstrained Scale
#'
#' @description
#' Each coordinate is drawn independently from \code{runif(min, max)} on the
#' unconstrained scale and mapped back through the bounds, so no draw is ever
#' rejected for being outside the box.
#'
#' @param min,max The range to draw from, in unconstrained units. Both default
#'   to a width of one either side of zero, and both may be given per parameter
#'   rather than as a single number.
#' @param npar The number of parameters. Defaults to \code{NULL}, meaning work
#'   it out from the bounds or from the objective; see \code{\link{minimize}}.
#'
#' @details
#' The range is in unconstrained units, which is what makes a single default
#' workable. A draw in \eqn{(-1, 1)} becomes a variance between \eqn{0.37} and
#' \eqn{2.7}, a probability between \eqn{0.27} and \eqn{0.73}, and a parameter
#' bounded on both sides lands well inside its interval; the same numbers on the
#' parameter scale would mean quite different things and would sometimes be
#' inadmissible.
#'
#' Widen it when the scale of the problem is unknown. \code{start_runif(-5, 5)}
#' spans four orders of magnitude for a positive parameter, which is usually more
#' than enough and is still a range no draw can fall out of.
#'
#' The draw uses \R's ordinary generator, so \code{\link{set.seed}} reproduces
#' it, and the seed is recorded in the result.
#'
#' @return An S7 object of class \code{UniformStart}.
#'
#' @examples
#' f <- function(p) sum((p - c(1, 2, 3))^2)
#' set.seed(1)
#' minimize(bfgs(), f, start_runif(npar = 3))@par
#'
#' # a wider net, and a positive parameter
#' set.seed(1)
#' minimize(bfgs(), function(p) (log(p) - 1)^2, start_runif(-5, 5, npar = 1),
#'          lower = 0)@par
#'
#' @seealso \code{\link{start_zeros}}, \code{\link{starting_values}}
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
#' @description Draws \code{npar} values from the starter's range.
#' @param starter A \code{UniformStart} object.
#' @param npar The number of parameters.
#' @return A numeric vector.
#' @keywords internal
S7::method(starting_values, UniformStart) <- function(starter, npar) {
  lo <- recycle_to(starter@min, npar, "min")
  up <- recycle_to(starter@max, npar, "max")
  stats::runif(npar, lo, up)
}


# --- resolution -------------------------------------------------------------

#' Validate a Declared Parameter Count
#'
#' @param npar \code{NULL} or a positive whole number.
#' @return \code{NULL}, or the value as an integer.
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
#' @param v A numeric vector.
#' @param n The length wanted.
#' @param nm The argument's name, for the message.
#' @return A numeric vector of length \code{n}.
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
#' which the objective accepts. Called by \code{\link{minimize}} when a starter
#' was given without \code{npar} and the bounds do not say.
#'
#' @param fn The objective.
#' @param gr Its gradient, or \code{NULL}. Supplying one makes the answer much
#'   more likely to be unique; see Details.
#' @param probe A function of one integer returning a candidate parameter vector
#'   of that length, so that the objective is probed where it will be used.
#' @param npar_max The largest length tried. Defaults to \code{50}.
#'
#' @details
#' A length is accepted when \code{fn} returns a single finite number for it and
#' raises neither an error nor a warning, and, if \code{gr} was supplied, when
#' the gradient comes back with the same length as its argument.
#'
#' What decides it is whether the objective genuinely rejects the wrong
#' length, and the objective a modeling package hands over usually does:
#' \code{X \%*\% beta} with a parameter of the wrong length is an error rather
#' than a number, so a regression of any kind is settled at once. A gradient
#' helps when it spells its components out, since such a gradient returns a fixed
#' number of them whatever it is handed.
#'
#' A vectorized objective written in terms of the parameter alone is another
#' matter: both plausible
#' guesses about \R are wrong. Recycling warns only when the shorter length is
#' not a \emph{divisor} of the longer, so \code{sum((p - c(1, 2, 3))^2)} accepts
#' a length-one vector in silence and returns a perfectly finite 14; and its
#' gradient \code{2 * (p - c(1, 2, 3))} returns six components for a length-six
#' argument, so the gradient rule passes it too. Rosenbrock's
#' \code{100 (p_2 - p_1^2)^2 + (1 - p_1)^2} accepts every length from two upwards
#' for the same kind of reason. None of this is a defect and none of it can be
#' guessed at.
#'
#' The probe therefore settles the objectives that have a
#' fixed width built into them and rejects the ones that do not, naming the two
#' lengths it found. When it rejects, \code{npar} or a vector of bounds is one
#' word.
#'
#' The search stops as soon as a \emph{second} length is accepted, because at
#' that point the answer is already known to be ambiguous and there is no reason
#' to keep probing. So the cost is two evaluations when the objective accepts any
#' length, and at most \code{npar_max} when it accepts exactly one. Either way it
#' happens once, before the run.
#'
#' @return A single integer.
#'
#' @examples
#' # a hand-written gradient pins it exactly
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' infer_npar(f, gr, function(k) numeric(k))
#'
#' # without one, the same objective is happy with any length from two upwards
#' try(infer_npar(f, NULL, function(k) numeric(k)))
#'
#' @seealso \code{\link{start_zeros}}, \code{\link{minimize}}
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
#' The whole of what \code{\link{minimize}} does with a starter, in one place:
#' settle the number of parameters, draw the values on the unconstrained scale,
#' and map them back through the box.
#'
#' @details
#' The number of parameters is looked for in three places, in order of how much
#' the caller was willing to say. \code{npar} on the starter itself is taken as
#' given. Failing that, a \code{lower} or \code{upper} of length greater than one
#' answers the question, since bounds are one per parameter. Failing both,
#' \code{\link{infer_npar}} probes the objective.
#'
#' A numeric \code{par} passes through untouched, so this costs nothing at all
#' for the ordinary call.
#'
#' @param par Whatever was passed as \code{par}.
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
