#' @include optimizer_class.R
NULL

# Modifications of an indefinite Hessian, as objects. newton() and chebyshev()
# hold one; the compiled loop runs the two shipped modifications in C++ and
# calls back into R for any other class, through hessian_modify().


#' @title S7 Class for Hessian Modifications
#'
#' @description
#' This class is the abstract parent of the objects that repair an indefinite
#' or nearly singular Hessian before a Newton step. [newton()] and
#' [chebyshev()] hold one in their `modify` argument. The two shipped
#' subclasses are built by [modify_eigen()] and [modify_ridge()]; a new
#' modification is a subclass of this class with a method for
#' [hessian_modify()].
#'
#' @details
#' A modification is used only at an iterate where the Cholesky factorization
#' of the Hessian fails. Where the factorization succeeds, the Newton step is
#' computed from the Hessian itself and no modification is involved. See
#' [hessian_modify()] for the requirements on a subclass and for an example of
#' a user-defined modification, and [modifications] for a comparison of the
#' shipped ones.
#'
#' @param label A short character label, shown when the optimizer that
#'   contains the modification is printed.
#'
#' @return An S7 object of class `hessian_modification`, carrying `label`.
#'   The class is abstract, so every value is an object of a subclass.
#'
#' @examples
#' # The class is abstract: use a constructor or a subclass.
#' try(hessian_modification(label = "mine"))
#'
#' modify_eigen()
#' modify_ridge(floor = 1e-6)
#'
#' @seealso [modify_eigen()], [modify_ridge()], [hessian_modify()],
#'   [modifications], [newton()], [chebyshev()]
#' @export
hessian_modification <- S7::new_class(
  "hessian_modification",
  properties = list(label = S7::class_character),
  abstract = TRUE
)


#' @title Modifications of the Hessian
#'
#' @description
#' Overview of the objects that repair a Hessian that is not positive
#' definite, passed to [newton()] and [chebyshev()] as `modify`.
#'
#' @details
#' The Newton direction \eqn{d = -H^{-1}g} is a descent direction only when
#' \eqn{H} is positive definite. When the Cholesky factorization of \eqn{H}
#' fails, the method replaces \eqn{H} by a positive definite matrix
#' \eqn{\tilde H} and solves \eqn{\tilde H d = -g}. The object in `modify`
#' selects how \eqn{\tilde H} is formed. The default is [modify_eigen()].
#'
#' \describe{
#'   \item{[modify_eigen()]}{Decomposes \eqn{H = Q \Lambda Q^\top} and raises
#'     every eigenvalue below
#'     \eqn{\max(\texttt{floor}, \texttt{floor}\,\max_i\lvert\lambda_i\rvert)}
#'     to that value. It costs a symmetric eigendecomposition and leaves the
#'     eigenvectors and the eigenvalues above the floor unchanged. The step
#'     component along an eigenvector \eqn{q_i} is
#'     \eqn{-q_i^\top g / \max(\lambda_i, \lambda_{\min})}, the Newton
#'     component where \eqn{\lambda_i \ge \lambda_{\min}} and a component
#'     of length \eqn{\lvert q_i^\top g\rvert/\lambda_{\min}} where the
#'     eigenvalue was raised, which is large when the floor is small.}
#'   \item{[modify_ridge()]}{Adds \eqn{\tau I} to \eqn{H}, starting from the
#'     smallest \eqn{\tau \ge \texttt{floor}} that makes the diagonal
#'     positive and doubling \eqn{\tau} until the Cholesky factorization
#'     succeeds. It costs a few factorizations and shifts every eigenvalue
#'     by the same amount.}
#'   \item{A user-defined subclass of [hessian_modification()]}{Any other
#'     rule, supplied as a method for [hessian_modify()]. See that page for
#'     an executable example.}
#' }
#'
#' Both shipped modifications run in compiled code. A user-defined class is
#' called from the compiled loop through R at each iteration at which the
#' Hessian is not positive definite.
#'
#' A repair is recorded in the trace as `hessian modified` (or
#' `hessian modified (capped)`) whichever modification was used; see
#' [newton()].
#'
#' @return Not applicable: this page documents the family. Each constructor
#'   returns an object inheriting from [hessian_modification()].
#'
#' @examples
#' f <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
#' g <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
#' h <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)
#'
#' # The Hessian at the start has a negative eigenvalue.
#' eigen(h(c(0.05, 0.6)), only.values = TRUE)$values
#'
#' sapply(list(eigen = modify_eigen(), ridge = modify_ridge()), function(m) {
#'   r <- minimize(newton(modify = m), f, c(0.05, 0.6), gr = g, he = h)
#'   c(value = r@value, iterations = r@iterations)
#' })
#'
#' @seealso [modify_eigen()], [modify_ridge()], [hessian_modification()],
#'   [hessian_modify()], [newton()], [chebyshev()]
#' @name modifications
NULL


#' @title S7 Class for the Eigenvalue Modification
#'
#' @description
#' This class holds the lower bound of the eigenvalues of the modified
#' Hessian. [modify_eigen()] builds it.
#'
#' @param floor The lower bound, as in [modify_eigen()].
#' @param label The label of every [hessian_modification()].
#'
#' @return An S7 object of class `ModifyEigen` inheriting from
#'   [hessian_modification()], with the property `floor` beside `label`.
#'
#' @seealso [modify_eigen()] for the constructor.
#' @name ModifyEigen-class
#' @aliases ModifyEigen
#' @keywords internal
ModifyEigen <- S7::new_class("ModifyEigen", parent = hessian_modification,
  properties = list(floor = S7::class_numeric))


#' @title S7 Class for the Ridge Modification
#'
#' @description
#' This class holds the smallest shift tried by the ridge modification.
#' [modify_ridge()] builds it.
#'
#' @param floor The smallest shift, as in [modify_ridge()].
#' @param label The label of every [hessian_modification()].
#'
#' @return An S7 object of class `ModifyRidge` inheriting from
#'   [hessian_modification()], with the property `floor` beside `label`.
#'
#' @seealso [modify_ridge()] for the constructor.
#' @name ModifyRidge-class
#' @aliases ModifyRidge
#' @keywords internal
ModifyRidge <- S7::new_class("ModifyRidge", parent = hessian_modification,
  properties = list(floor = S7::class_numeric))


#' @title Eigenvalue Modification of the Hessian
#'
#' @description
#' Builds the modification that decomposes the Hessian,
#' \eqn{H = Q \Lambda Q^\top}, and raises every eigenvalue below a floor to
#' the floor:
#' \deqn{\tilde H = Q \max(\Lambda, \lambda_{\min}) Q^\top, \qquad
#'   \lambda_{\min} = \max(\texttt{floor},\;
#'   \texttt{floor}\,\max_i\lvert\lambda_i\rvert).}
#' The eigenvectors and the eigenvalues above the floor are unchanged. The
#' step component along an eigenvector \eqn{q_i} is
#' \eqn{-q_i^\top g / \max(\lambda_i, \lambda_{\min})}: the Newton component
#' where \eqn{\lambda_i \ge \lambda_{\min}}, and, where the eigenvalue was
#' raised, a component of length \eqn{\lvert q_i^\top g\rvert/\lambda_{\min}}
#' that is large when the floor is small. The cost is a symmetric
#' eigendecomposition. It is the default of
#' [newton()] and [chebyshev()]. If the eigendecomposition fails, the ridge
#' modification is used.
#'
#' @param floor The lower bound of the eigenvalues, a single positive number.
#'   Defaults to `1e-8`. The bound is relative to the largest eigenvalue in
#'   absolute value when that value exceeds one.
#'
#' @return An S7 object of class `ModifyEigen`, inheriting from
#'   [hessian_modification()], with the properties `floor` and `label`.
#'
#' @examples
#' m <- modify_eigen(floor = 1e-6)
#' m
#'
#' # The method for hessian_modify() applied to an indefinite matrix.
#' H <- matrix(c(2, 1, 1, -3), 2, 2)
#' eigen(H, only.values = TRUE)$values
#' eigen(hessian_modify(m, H), only.values = TRUE)$values
#'
#' # Passed to the optimizer.
#' f <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
#' g <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
#' h <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)
#' minimize(newton(modify = modify_eigen(floor = 1e-6)), f, c(0.05, 0.6),
#'          gr = g, he = h)
#'
#' @seealso [modify_ridge()], [modifications] for the comparison,
#'   [hessian_modify()], [newton()], [chebyshev()]
#' @export
modify_eigen <- function(floor = 1e-8) {
  check_modify_floor(floor)
  ModifyEigen(label = paste0("eigenvalue floor (floor = ", format(floor), ")"),
              floor = floor)
}


#' @title Ridge Modification of the Hessian
#'
#' @description
#' Builds the modification that adds a multiple of the identity,
#' \eqn{\tilde H = H + \tau I}. The shift starts at
#' \eqn{\tau = \max(\texttt{floor},\; \texttt{floor} - \min_i H_{ii})}, which
#' makes the diagonal positive, and is doubled until the Cholesky
#' factorization of \eqn{\tilde H} succeeds (at most 60 times; the identity
#' matrix is used if none succeeds). This is the idea of Levenberg: the
#' direction moves from the Newton direction at \eqn{\tau = 0} to a scaled
#' steepest-descent direction for large \eqn{\tau}. It costs less than
#' [modify_eigen()] and shifts every eigenvalue, not only the small ones.
#'
#' @param floor The smallest shift, a single positive number. Defaults to
#'   `1e-8`.
#'
#' @return An S7 object of class `ModifyRidge`, inheriting from
#'   [hessian_modification()], with the properties `floor` and `label`.
#'
#' @examples
#' m <- modify_ridge(floor = 1e-6)
#' m
#'
#' # The method for hessian_modify() applied to an indefinite matrix.
#' H <- matrix(c(2, 1, 1, -3), 2, 2)
#' eigen(hessian_modify(m, H), only.values = TRUE)$values
#'
#' # Passed to the optimizer.
#' f <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
#' g <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
#' h <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)
#' minimize(newton(modify = modify_ridge(floor = 1e-6)), f, c(0.05, 0.6),
#'          gr = g, he = h)
#'
#' @seealso [modify_eigen()], [modifications] for the comparison,
#'   [hessian_modify()], [newton()], [chebyshev()]
#' @export
modify_ridge <- function(floor = 1e-8) {
  check_modify_floor(floor)
  ModifyRidge(label = paste0("ridge (floor = ", format(floor), ")"),
              floor = floor)
}


#' Validate the Floor of a Modification
#'
#' @param floor The value supplied.
#'
#' @return Invisibly `TRUE`. Raises an error naming `floor` otherwise.
#'
#' @keywords internal
check_modify_floor <- function(floor) {
  if (length(floor) != 1L || !is.numeric(floor) || is.na(floor) || floor <= 0) {
    stop("'floor' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate the Modification of a Newton Method
#'
#' @description
#' Checks that `x` inherits from [hessian_modification()]. A character string
#' signals an error that names the constructor to use.
#'
#' @param x The value supplied.
#'
#' @return Invisibly `TRUE`. Raises an error naming `modify` otherwise.
#'
#' @keywords internal
check_modify <- function(x) {
  if (S7::S7_inherits(x, hessian_modification)) return(invisible(TRUE))
  stop("'modify' must be a hessian_modification object such as ",
       "modify_eigen().", call. = FALSE)
}


#' @title Modify a Hessian
#'
#' @description
#' Generic that returns a positive definite matrix \eqn{\tilde H} in place of
#' a Hessian \eqn{H} for which the Cholesky factorization failed. [newton()]
#' and [chebyshev()] solve \eqn{\tilde H d = -g} for the direction. There is
#' one method for each [hessian_modification()] subclass; a new modification
#' is a subclass and a method, and [newton()] needs no change.
#'
#' @details
#' # Requirements on a method
#'
#' The method receives the modification object and the Hessian \eqn{H}, a
#' symmetric numeric matrix with finite entries, and returns a numeric matrix
#' of the same dimension. The returned matrix should be positive definite:
#' the method of the optimizer solves with it, and a direction that is not a
#' descent direction is replaced by a scaled gradient step, which is recorded
#' in the trace as `hessian solve failed`. A returned matrix of another
#' dimension or with non-finite entries is treated in the same way. The
#' method is called only at an iterate at which \eqn{H} is not positive
#' definite, and is called from compiled code through R, once per such
#' iterate.
#'
#' # Writing a modification
#'
#' A modification is defined by (1) a subclass of [hessian_modification()]
#' with the parameters of the rule as properties, (2) a constructor that sets
#' the `label` shown when the optimizer is printed, and (3) a method for
#' this generic, registered with `S7::method()`. The shipped classes
#' [modify_eigen()] and [modify_ridge()] are handled in compiled code; any
#' other class, including a subclass of one of them, is called through this
#' generic.
#'
#' @param modify A [hessian_modification()] object.
#' @param H A symmetric numeric matrix.
#'
#' @return A numeric matrix with the dimension of `H`, intended to be positive
#'   definite.
#'
#' @examples
#' # The shipped methods, on an indefinite matrix with eigenvalues 2 and -1.
#' H <- matrix(c(0.5, 1.5, 1.5, 0.5), 2, 2)
#' eigen(H, only.values = TRUE)$values
#' hessian_modify(modify_eigen(floor = 0.1), H)
#' hessian_modify(modify_ridge(), H)
#'
#' # A user-defined modification: replace each eigenvalue by its absolute
#' # value, with a lower bound (the saddle-free Newton method of Dauphin et
#' # al., 2014). Unlike the eigenvalue floor, it uses the magnitude of a
#' # negative curvature as the curvature of the quadratic model.
#' ModifyAbs <- S7::new_class("ModifyAbs", parent = hessian_modification,
#'   properties = list(floor = S7::class_numeric))
#'
#' modify_abs <- function(floor = 1e-8) {
#'   ModifyAbs(label = paste0("absolute eigenvalues (floor = ", floor, ")"),
#'             floor = floor)
#' }
#'
#' S7::method(hessian_modify, ModifyAbs) <- function(modify, H) {
#'   e <- eigen(H, symmetric = TRUE)
#'   e$vectors %*% (pmax(abs(e$values), modify@floor) * t(e$vectors))
#' }
#'
#' hessian_modify(modify_abs(), H)
#'
#' # The new object is passed to newton() like a shipped one, here on the
#' # Himmelblau function, started where its Hessian is not positive definite.
#' k <- test_problems("himmelblau")[[1]]
#' eigen(k$he(c(0.5, 0.5)), only.values = TRUE)$values
#' r <- minimize(newton(modify = modify_abs(), keep_trace = TRUE), k$fn,
#'               c(0.5, 0.5), gr = k$gr, he = k$he)
#' c(value = r@value, iterations = r@iterations, converged = r@converged)
#' unique(r@trace$safeguard)
#'
#' @seealso [hessian_modification()], [modify_eigen()], [modify_ridge()],
#'   [modifications], [newton()], [chebyshev()]
#' @references
#' Dauphin, Y. N., Pascanu, R., Gulcehre, C., Cho, K., Ganguli, S. and
#' Bengio, Y. (2014). Identifying and attacking the saddle point problem in
#' high-dimensional non-convex optimization. *Advances in Neural Information
#' Processing Systems* **27**.
#' @export
hessian_modify <- S7::new_generic("hessian_modify", "modify",
                                  function(modify, H) S7::S7_dispatch())

#' @rdname hessian_modify
#' @name hessian_modify.ModifyEigen
#' @keywords internal
S7::method(hessian_modify, ModifyEigen) <- function(modify, H) {
  H <- 0.5 * (H + t(H))
  e <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e) || anyNA(e$values)) {
    return(hessian_modify(ModifyRidge(label = "ridge", floor = modify@floor), H))
  }
  lo <- max(modify@floor, modify@floor * max(abs(e$values)))
  e$vectors %*% (pmax(e$values, lo) * t(e$vectors))
}

#' @rdname hessian_modify
#' @name hessian_modify.ModifyRidge
#' @keywords internal
S7::method(hessian_modify, ModifyRidge) <- function(modify, H) {
  H <- 0.5 * (H + t(H))
  p <- nrow(H)
  tau <- max(modify@floor, -min(diag(H)) + modify@floor)
  for (k in seq_len(60L)) {
    Hm <- H + tau * diag(p)
    ok <- tryCatch({chol(Hm); TRUE}, error = function(e) FALSE)
    if (ok) return(Hm)
    tau <- 2 * tau
  }
  diag(p)
}


#' Describe a Modification to the Compiled Loop
#'
#' @description
#' Returns the part of the `method` list of [newton()] and [chebyshev()] that
#' describes the Hessian modification: `modify_type` (`"eigen"`, `"ridge"` or
#' `"custom"`), `floor`, and `modify_fn`, `NULL` for the two shipped
#' modifications and otherwise a function of the Hessian that calls
#' [hessian_modify()].
#'
#' @param modify A [hessian_modification()] object.
#'
#' @return A list of three.
#'
#' @keywords internal
modify_spec <- function(modify) {
  cls <- S7::S7_class(modify)
  if (identical(cls, ModifyEigen)) {
    list(modify_type = "eigen", floor = modify@floor, modify_fn = NULL)
  } else if (identical(cls, ModifyRidge)) {
    list(modify_type = "ridge", floor = modify@floor, modify_fn = NULL)
  } else {
    list(modify_type = "custom", floor = 0, modify_fn = function(H) {
      hessian_modify(modify, H)
    })
  }
}


#' @title Print Method for Hessian Modifications
#' @name print.hessian_modification
#' @param x A [hessian_modification()] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(modify_ridge())
#' @keywords internal
S7::method(print, hessian_modification) <- function(x, ...) {
  cat("<hessian_modification> ", x@label, "\n", sep = "")
  invisible(x)
}
