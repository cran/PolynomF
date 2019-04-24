
#' Simpl orthogonal polynomials
#'
#' Generate a list of polynomials up to a specified degree,
#' orthogonal with respect to the natural inner product on
#' a discrete, finite set of x-values with equal weights.
#'
#' @param x A numeric vector
#' @param degree The desired maximum degree
#' @param norm Logical: should polynomials be normalised to length one?
#' @param ... Arguments passed on to the non-deprecated function.
#'
#' @return A list of orthogonal polynomials as a polylist object
#' @export
#'
#' @examples
#' x <- c(0:3, 5)
#' P <- poly_orth(x)
#' plot(P, lty = "solid")
#' Pf <- as.function(P)
#' zap(crossprod(Pf(x)))
poly_orth <- function(x, degree = length(unique(x)) - 1, norm = TRUE) {
  at <- attr(poly(x, degree), "coefs")
  a <- at$alpha
  N <- at$norm2
  x <- polynom()
  p <- polylist(polynom(0), polynom(1))
  for(j in 1:degree) {
    p[[j + 2]] <- (x - a[j]) * p[[j + 1]] - N[j + 1]/N[j] * p[[j]]
  }
  p <- p[-1]
  if(norm) {
    sqrtN <- sqrt(N[-1])
    for(j in 1 + 0:degree) {
      p[[j]] <- p[[j]]/sqrtN[j]
    }
  }
  setNames(p, paste0("P", 0:degree))
}

#' @rdname poly_orth
#' @export
poly.orth <- function(...) {
  .Deprecated("poly_orth")
  poly_orth(...)
}

# #' @rdname poly_orth
# #' @param weights A numeric vector of non-negative weights (default: all 1).
# #' @param ... Arguments transferred to \code{poly_orth}.
# #' @export
# poly_orth_weighted <- local({
#   x <- 0
#   w <- 1
#   ip <- function(p1, p2 = p1) {
#     sum(w*p1(x)*p2(x))
#   }
#   normalize <- function(p) {
#     p/sqrt(ip(p))
#   }
#   function(x, degree = length(unique(x[w > 0])) - 1, weights = 1,
#            norm = TRUE) {
#
#     weights <- rep_len(weights, length.out = length(x))
#     x <<- x
#     stopifnot(is.numeric(x) &&
#                 (length(unique(x)) > 0))
#     stopifnot(is.numeric(weights) &&
#                 all(weights >= 0))
#     w <<- weights/mean(weights)
#     stopifnot(is.numeric(degree) &&
#                 (length(degree) == 1) &&
#                 ((degree %% 1) == 0) &&
#                 (degree >= 0))
#     OP <- polylist(polynom(0), polynom(1))
#     if(degree == 0) {
#       return(OP[-1])
#     }
#     z <- polynom()
#     pn1 <- 1
#     pnn <- ip(OP[[2]])
#     for(j in 1:degree) {
#       an <- ip(OP[[j+1]], z*OP[[j+1]])/pnn
#       bn <- pnn/pn1
#       OP[[j+2]] <- (z - an)*OP[[j+1]] - bn*OP[[j]]
#       pn1 <- pnn
#       pnn <- ip(OP[[j+2]])
#     }
#     OP <- OP[-1]
#     if(norm) {
#       OP <- as_polylist(lapply(OP, normalize))
#     }
#     OP
#   }
# })

#' General Orthogonal Polynomials
#'
#' Generate sets of polynomials orthogonal with respect to a general
#' inner product.  The inner product is specified by an R function of
#' (at least) two polynomial arguments.
#'
#' Discrete orthogonal polynomials, equally or unequally weighted,
#' are included as special cases.  See the \code{Discrete} inner
#' product function.
#'
#' Computations are done using
#' the recurrence relation with computed coefficients.  If the algebraic
#' expressions for these recurrence relation coefficients are known
#' the computation can be made much more efficient.
#'
#' @param inner_product An R function of two \code{"polynom"} arguments
#'        with the second polynomial having a default value equal to the
#'        first.  Additional arguments may be specified.  See examples
#' @param degree A non-negative integer specifying the maximum degree
#' @param norm  Logical: should the polynomials be normalized?
#' @param ... additional arguments passed on to the inner product function
#' @param p,q Polynomials
#' @param alpha,beta Family parameters for the Jacobi polynomials
#' @param x numeric vector defining discrete orthogonal polynomials
#' @param w a weight function for discrete orthogonal polynomials
#'
#' @return A \code{"polylist"} object containing the orthogonal set
#' @export
#'
#' @examples
#' (P0 <- poly_orth(0:5, norm = FALSE))
#' (P1 <- poly_orth_general(Discrete, degree = 5, x = 0:5, norm = FALSE))
#' sapply(P0-P1, function(x) max(abs(coef(x))))  ## visual check for equality
#' (P0 <- poly_orth_general(Legendre, 5))
#'    ### should be same as P0, up to roundoff
#' (P1 <- poly_orth_general(Jacobi, 5, alpha = 0, beta = 0))
#'                 ### check
#' sapply(P0-P1, function(x) max(abs(coef(x))))
poly_orth_general <- function(inner_product, degree, norm = FALSE, ...) {
  stopifnot(is.numeric(degree) &&
              (length(degree) == 1) &&
              ((degree %% 1) == 0) &&
              (degree >= 0))
  OP <- polylist(polynom(0), polynom(1))
  if(degree == 0) {
    return(OP[-1])
  }
  z <- polynom()
  pn1 <- 1
  pnn <- inner_product(OP[[2]], ...)
  for(j in 1:degree) {
    an <- inner_product(OP[[j+1]], z*OP[[j+1]], ...)/pnn
    bn <- pnn/pn1
    OP[[j+2]] <- (z - an)*OP[[j+1]] - bn*OP[[j]]
    pn1 <- pnn
    pnn <- inner_product(OP[[j+2]], ...)
  }
  OP <- OP[-1]
  if(norm) {
    OP <- as_polylist(lapply(OP, function(p) p/sqrt(inner_product(p, p, ...))))
  }
  setNames(OP, paste0("P", 0:degree))
}

#' @rdname poly_orth_general
#' @export
Hermite <- function(p, q = p) {
  integrate(function(x) dnorm(x)*p(x)*q(x), lower = -Inf, upper = Inf)$value
}

#' @rdname poly_orth_general
#' @export
Legendre <- function(p, q = p) {
  integral(p*q, limits = c(-1,1))
}

#' @rdname poly_orth_general
#' @export
Chebyshev <- function(p, q = p) {
  integrate(function(x) sqrt(1 - x^2)*p(x)*q(x), lower = -1, upper = 1)$value
}

#' @rdname poly_orth_general
#' @export
Jacobi <- function(p, q = p, alpha = 0.5, beta = 0.5) {
  integrate(function(x) (1-x)^alpha*(1 + x)^beta*p(x)*q(x), lower = -1, upper = 1)$value
}

#' @rdname poly_orth_general
#' @export
Discrete <- function(p, q = p, x, w = function(x) rep(1, length(x))) {
  sum(w(x)*p(x)*q(x))
}


#' Remove minuscule coefficients
#'
#' A convenience function for setting polynomial coefficients likely
#' to be entirely round-off error to zero.  The decision is relegated
#' to the function \code{base::zapsmall}, to which this is a front-end.
#'
#' @param x A polynomial or polylist object
#' @param digits As for \code{base::zapsmall}
#' @param ... Passed on to \code{base::zapsmall}
#'
#' @return A polynomial or polylist object with minuscule coefficients set to zero.
#' @export
#'
#' @examples
#' (P <- poly_orth(-2:2, norm = FALSE))
#' zap(P)
zap <- function(x, digits = getOption("digits")) {
  UseMethod("zap")
}

#' @rdname zap
#' @export
zap.default <- base::zapsmall

#' @rdname zap
#' @export
zap.polynom <- function(x, digits = getOption("digits")) {
  polynom(base::zapsmall(coef(x), digits))
}

#' @rdname zap
#' @export
zap.polylist <- function(x, digits = getOption("digits")) {
  as_polylist(lapply(x, zap.polynom, digits))
}
