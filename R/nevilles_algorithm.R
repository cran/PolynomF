

#' Lagrange Interpolation Polynomials
#'
#' Compute the Lagrange Interpolation Polynomial from a given set of x- and y-values,
#' or, alterntively, compute the interpolated values at a set of given x-values.
#' Two algorithms are provided, namely Neville's algorithm, or a more direct version
#' based on the usual Lagrange formula.  The latter is generally faster but the former
#' can be more accurate numerically.
#'
#' @param x A numeric vector of x-values
#' @param y A numeric values of y-values corresponding to the x-values
#' @param x0 Either a polynomial object or a vector of x-values for which
#'           interpolated y-values are required.
#'
#' @return Either an interpolation polynomial object or a vector of interpolated y-values
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' x <- 1:5
#' y <- rnorm(x)
#' xout <- 0.5 + 1:4
#'
#' p1 <- neville(x, y)
#' plot(p1, xlim = range(x), ylim = extendrange(y, f = 1), panel.first = grid())
#' points(x, y, col = 4)
#' points(xout, lagrange(x, y, xout), col = 2)
neville <- function(x, y, x0 = polynomial()) {
  stopifnot("bad x" = is.numeric(x) && length(x) > 1,
            "bad y" = is.numeric(y) && length(y) == length(x))
  ox <- order(x, decreasing = TRUE)
  x <- x[ox]
  y <- y[ox]
  dup <- duplicated(x)
  if(any(dup)) {
    x <- x[!dup]
    y <- y[!dup]
    warning("data from some duplicated x-values discarded")
  }
  n <- length(x)
  stopifnot("x is too short" = n > 1)
  p0 <- as.list(y)
  for(gap in 1:(n-1)) {
    p <- list()
    for(i in 1:(n - gap)) {
      j <- i+gap
      p[[i]] <- ((x[i] - x0)*p0[[i + 1]] + (x0 - x[j])*p0[[i]])/(x[i] - x[j])
    }
    p0 <- p
  }
  p[[1]]
}

#' @rdname neville
#' @export
lagrange <- function(x, y, x0 = polynomial()) {
  stopifnot("bad x" = is.numeric(x) && length(x) > 1,
            "bad y" = is.numeric(y) && length(y) == length(x))
  ox <- order(x, decreasing = TRUE)
  x <- x[ox]
  y <- y[ox]
  dup <- duplicated(x)
  if(any(dup)) {
    x <- x[!dup]
    y <- y[!dup]
    warning("data from some duplicated x-values discarded")
  }
  stopifnot("x is too short" = length(x) > 1)
  P <- 1
  for(xi in x) {
    P <- P*(x0 - xi)
  }
  Q <- 0
  for(i in seq_along(y)) {
    Q <- Q + y[i]*(P/(x0 - x[i]))/prod(x[i] - x[-i])
  }
  Q
}
