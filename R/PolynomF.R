## revision of the polynomial class with a different representation
#' @useDynLib PolynomF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @importFrom stats deriv poly predict coef
#' @importFrom graphics lines par plot points grid
NULL

#' Polynomial construction
#'
#' Functions to construct polynomial objects and check class membership
#'
#' @param a  A \code{polynom} object, or a numeric vector of coefficients
#'           (in "power series" order) or a vector object which can be
#'           coerced to one.
#' @param x  An object of class \code{"polylist"}, at least potentially.
#' @param ... Additional arguments, currently ignored.
#' @param eps A small non-negative tolerance to check for zero components.
#'
#' @return A polynomial object.
#' @export
#'
#' @examples
#' (x <- polynom())
#' (p <- polynom(c(1, 5, 2, 2)/10))
#' plot(p, xlim = 0:1, ylim = 0:1)
#' P <- p
#' for(j in 1:7) {
#'   lines(P, col = j)
#'   P <- p(P)
#' }
#' (r <- solve(p-x))
#' segments(r, 0, r, p(r), lty = "dashed")
polynom <- function(a = c(0,1), ..., eps = 0) {
  ### constructor function
  a <- as.numeric(a)
  if(any(is.na(a)))
    stop("illegal coefficient vector")
  a[abs(a) < eps] <- 0
  while((k <- length(a)) > 1 &&
    abs(a[k]) <= eps) a <- a[-k]
  structure(function(x) {
    p <- 0
    for(a0 in rev(a)) p <- a0 + x*p
    p
  }, class = "polynom")
}

#' @rdname polynom
#' @export
as.polynom <- function(a) {
  ### coercion to polynom
  if(is.polynom(a)) a else polynom(as.vector(a))
}

#' @rdname polynom
#' @export
is.polynom <- function(a) {
  ### predicate function
  inherits(a, "polynom")
}

.tangent <- function(x0, p) {
  x <- polynom()
  p(x0) + deriv(p)(x0)*(x - x0)
}

#' Tangent lines
#'
#' Find the tangent line to a polynomial at one or more x-points
#'
#' @param p  A polynomial object
#' @param x0 A numeric vector of values at which the tangent line(s) are required
#'
#' @return A linear polynomial giving the tangent line, or a list of such polynomials
#' @export
#'
#' @examples
#' p <- poly.from.zeros(c(0, 0:5, 4))
#' plot(p, xlab = expression(italic(x)), ylab = expression(italic(P(x))),
#'   main = parse(text = paste("italic(P(x) ==",
#'                              as.character(p, decreasing = TRUE),")")))
#' x0 <- solve(deriv(p))        ## stationary points
#' lines(tangent(p, x0), col = "dark green", lty = "solid",
#'       limits = cbind(x0-1/4, x0+1/4))
#' points(x0, p(x0), col = "dark green")
#'
#' x0 <- solve(deriv(deriv(p))) ## points of inflexion
#' lines(tangent(p, x0), col = "red", lty = "solid", lwd = 2,
#'       limits = cbind(x0-1/4, x0+1/4))
#' points(x0, p(x0), col = "red")
#' legend("bottomleft", c("Stationary points", "Points of inflexion"),
#'        pch = 19, col = c("dark green", "red"), lty = "solid",
#'        cex = 0.7, bg = "beige", box.lwd = 0.25)
tangent <- function(p, x0) {   ## the tangent line(s) to p(x) at x = x0
  stopifnot(is.polynom(p), is.numeric(x0) && length(x0) > 0)
  if(length(x0) == 1) {
    .tangent(x0, p)
  } else {
    as.polylist(lapply(x0, .tangent, p = p))
  }
}

.poly.mult <- function(e1, e2) {
  poly_product(e1, e2)
}

.poly.quo <- function(e1, e2) {### quotient
  poly_divide(e1, e2)$quotient
}

.poly.rem <- function(e1, e2) { ### remainder
  poly_divide(e1, e2)$remainder
}

#' Polynomial arithmetic
#'
#' Group generic function to implement arithmetic operations on polynomial objects
#'
#' @param e1,e2 A numeric vector of a polynomial object. At least one of \code{e1}
#'              or \code{e2} must be an object of class \code{"polynom"} or
#'              \code{"polylist"}.
#'
#' @return A polynomial or polylist object representing the result of the operation.
#' @export
#'
#' @examples
#' x <- polynom()
#' (p <- (x-1)^5 - 1)
#' (p1 <- (p + 1)/(x - 1)^2 - 1)
#' for(i in 0:10) cat(coef((x+1)^i), "\n")
Ops.polynom <- function(e1, e2) {
  if(missing(e2))   ### unary operations
    return(switch(.Generic,
            "+" = e1,
            "-" = polynom(-coef(e1)),
            stop("unsupported unary operation")))
  e1 <- if(is.polynom(e1)) coef(e1) else as.numeric(e1)
  e2 <- if(is.polynom(e2)) coef(e2) else as.numeric(e2)
  l1 <- length(e1)
  l2 <- length(e2)
  e1.op.e2 <-
    switch(.Generic,
             ### addition and subtraction
         "+" = ,
         "-" = {
           e1 <- c(e1, rep.int(0, max(0, l2 - l1)))
           e2 <- c(e2, rep.int(0, max(0, l1 - l2)))
           NextMethod(.Generic)
         },
             ### product of two polynomials
         "*" = .poly.mult(e1, e2),
             ### quotient
         "/" =,
         "%/%" = .poly.quo(e1, e2),
             ### remainder
         "%%" = .poly.rem(e1, e2),
             ### non-negative integer powers
         "^" = {
           if(l2 != 1 || e2 < 0 || e2 %% 1 != 0)
             stop("unsupported polynom power")
           p <- 1
           while(e2 > 0) {
             if (e2 %% 2 == 1) {
               p <- .poly.mult(p, e1)
               e2 <- e2 - 1
             }
             e1 <- .poly.mult(e1, e1)
             l1 <- length(e1)
             e2 <- e2 / 2
           }
           p
         },
             ### equality and inequality
         "==" = return(l1 == l2 && all(e1 == e2)),
         "!=" = return(l1 != l2 || any(e1 != e2)),
         stop("unsupported operation on polynoms"))
  polynom(e1.op.e2)
}

#' @rdname Ops.polynom
#' @export
Ops.polylist <- function(e1, e2) {
  if(missing(e2))
    return(switch(.Generic,
      "+" = e1,
      "-" = as.polylist(lapply(e1, "-")),
      stop("unknown unary operator!")))
  switch(.Generic,
  "+" =,
  "-" =,
  "*" =,
  "/" =,
  "%/%" =,
  "%%" = as.polylist(mapply(.Generic, lapply(e1, as.polynom), lapply(e2, as.polynom))),
  "^" = as.polylist(mapply(.Generic, lapply(e1, as.polynom), e2)),
  "==" =,
  "!=" = unlist(mapply(.Generic, lapply(e1, as.polynom), lapply(e2, as.polynom))),
  stop("unsupported operation on polynoms"))
}

.accumulate <- function(f, init, x, right = TRUE) {
  ## .accumulate a la Abelson and Sussman.
  if(length(x) == 0)
    return(init)
  f <- match.fun(f)
  if(right)
    f(x[[1]], Recall(f, init, x[-1], right = TRUE))
  else
    Recall(f, f(init, x[[1]]), x[-1], right = FALSE)
}

#' Summary and Math methods for polynomials
#'
#' These provide methods for the generic function \code{Summary}
#' and \code{Math} for polynomial and polylist objects.  For \code{Summary}
#' only \code{sum} and \code{prod} members are implemented
#'
#' @name GroupGenerics
#' @param x a \code{"polynom"} or \code{"polylist"} objects.
#' @param ... Additional arguments
#' @param na.rm Logical: should missing values be removed?
#'
#' @return The result of the group generic operation
#' @export
#'
#' @examples
#' lis <- as.polylist(lapply(-2:3, function(x) polynom() - x))
#' prod(lis)
#' sum(lis)
#' solve(prod(lis))
#' solve(sum(lis))
Summary.polynom <- function(..., na.rm = FALSE) {
  ok <- switch(.Generic,
               sum = , prod = TRUE,
               FALSE)
  if(!ok)
    stop(gettextf("Generic '%s' not defined for '%s' objects.",
                  .Generic, .Class))
  switch(.Generic,
         "sum" = .accumulate("+", polynom(0), polylist(...)),
         "prod" = .accumulate("*", polynom(1), polylist(...)))
}

#' @rdname GroupGenerics
#' @export
Summary.polylist <- function(..., na.rm = FALSE) {
  ok <- switch(.Generic,
               sum = , prod = TRUE,
               FALSE)
  if(!ok)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class))
  switch(.Generic,
         "sum" = .accumulate("+", polynom(0), c(...)),
         "prod" = .accumulate("*", polynom(1), c(...)))
}

#' @rdname GroupGenerics
#' @export
Math.polynom <- function(x, ...) {
  x <- coef(x)
  switch(.Generic,
         round = ,
         signif = ,
         floor = ,
         ceiling = ,
         trunc = polynom(NextMethod(.Generic)),
         stop(paste(.Generic, "unsupported for polynoms")))
}

#' @rdname GroupGenerics
#' @export
Math.polylist <- function(x, ...) {
  sapply(x, .Generic, ...)
}

#' Polynomial coercion to character
#'
#' Produce a text representation of a polynomial object
#'
#' @param x The polynomial object in question
#' @param variable Character string: what variable name should be used?
#' @param decreasing Logical: in decreasing powers or increasing powers?
#' @param ... Additional arguments (ignored as yet)
#'
#' @return A character string representation of the polynomial
#' @export
#'
#' @examples
#' p <- poly.from.zeros(-2:3)
#' as.character(p, "z", FALSE)
#' as.character(p, "z", TRUE)
#' parse(text = as.character(p, "z", TRUE))[[1]]
as.character.polynom <- function(x, variable = "x", decreasing = FALSE, ...) {
  if(is.polynom(x)) p <- coef(x) else p <- unclass(x)
  lp <- length(p) - 1
  names(p) <- 0:lp
  p <- p[p != 0]

  if(length(p) == 0) return("0")

  if(decreasing) p <- rev(p)

  signs <- ifelse(p < 0, "- ", "+ ")
  signs[1] <- if(signs[1] == "- ") "-" else ""

  np <- names(p)
  p <- as.character(abs(p))
  p[p == "1" & np != "0"] <- ""

  pow <- paste(variable, "^", np, sep = "")
  pow[np == "0"] <- ""
  pow[np == "1"] <- variable
  stars <- rep.int("*", length(p))
  stars[p == "" | pow == ""] <- ""
  paste(signs, p, stars, pow, sep = "", collapse = " ")
}

#' Print method for polynomial objects
#'
#' Standard method for printing polynomial objects
#'
#' @param x A polynomial object
#' @param variable Character string: what variable name should be given?
#' @param digits Integer: how many decimal degits to use?
#' @param decreasing Logical: in descending powers, or ascending?
#' @param ... Additional arguments
#'
#' @return The original object \code{x}, invisibly
#' @export
print.polynom <- function(x, variable = "x",
                          digits = getOption("digits"),
                          decreasing = FALSE, ...) {
    x <- coef(x)
    p <- as.character.polynom(signif(x, digits = digits),
                              variable = variable, decreasing = decreasing, ...)
    pc <- nchar(p)
    ow <- max(35, getOption("width"))
    m2 <- 0
    while(m2 < pc) {
      m1 <- m2 + 1
      m2 <- min(pc, m2 + ow)
      if(m2 < pc)
        while(substring(p, m2, m2) != " " && m2 > m1 + 1)
          m2 <- m2 - 1
      cat(substring(p, m1, m2), "\n")
    }
    invisible(x)
  }

#' Coercion to function
#'
#' PolynomF objects ARE functions, but this coercion method
#' creates from a polynomial object a pure function with the
#' coefficients fully exposed in the code and which evaluates
#' the polynomial more efficiently.
#'
#' @param x A polynomial object
#' @param variable Character string: what variable name should be used?
#' @param ... Additional arguments
#'
#' @return An explicit R function evaluating the polynomial
#' @export
#'
#' @examples
#' p <- poly.from.zeros(-2:3)
#' p
#' as.function(p)
as.function.polynom <- function (x, variable = "x", ...) {
  a <- rev(coef(x))
  w <- as.name(".w")
  v <- as.name(variable)
  ex <- call("{", call("<-", w, 0))
  for (i in seq_along(a)) {
    ex[[i + 2]] <- call("<-", w, call("+", a[1], call("*", v, w)))
    a <- a[-1]
  }
  ex[[length(ex) + 1]] <- w
  f <- quote(function(x) NULL)
  names(f[[2]]) <- variable
  f <- eval(f)
  body(f) <- ex
  f
}

#' @rdname as.function.polynom
#' @export
as.function.polylist <- function(x, ...) {
  x <- lapply(x, as.function.polynom)
  function(z, ...) sapply(x, function(p) p(z), ...)
}

#' Polynomial coefficients
#'
#' Extract polynomial coefficients
#'
#' @param object A polynomial object or list thereof
#' @param ... Ignored
#'
#' @return A numeric vector of coefficients
#' @export
#'
#' @examples
#' p <- polynom(1:3)*polynom(5:1)
#' coef(p)
coef.polynom <- function(object,...) {
  get("a", envir = environment(object))
}

#' @rdname coef.polynom
#' @export
coef.polylist <- function(object, ...) {
  sapply(object, coef.polynom, ...)
}

#' Polynomial Calculus
#'
#' Find the derivative or indefinite integral of a polynomial object, or list thereof.
#'
#' @param expr A polynomial object, or list thereof
#' @param limits Real limits of a definite integral
#' @param ... Unused as yet
#'
#' @return A coeffieient vector, or list thereof
#' @export
#'
#' @examples
#' p <- poly.from.roots(-2:3)
#' p
#' deriv(p)
#' integral(p)
deriv.polynom <- function(expr, ...) {
  expr <- coef(expr)
  if(length(expr) == 1)
    return(polynom(0))
  expr <- expr[-1]
  polynom(expr * seq(along = expr))
}

#' @rdname deriv.polynom
#' @export
integral <- function(expr, ...) {
  UseMethod("integral")
}

#' @rdname deriv.polynom
#' @export
integral.default <- function(expr, ...) {
  stop(gettextf("No 'integral' method for objects of class '%s'",
                class(expr)))
}

#' @rdname deriv.polynom
#' @export
integral.polynom <- function(expr, limits = NULL, ...) {
  expr <- coef(expr)
  p <- polynom(c(0, expr/seq(along = expr)))
  if(is.null(limits))
    p
  else
    diff(p(limits))
}

#' Orthogonal polynomials
#'
#' Generate a list of polynomials up to a specified degree,
#' orthogonal with respect to the natural inner product
#'
#' @param x A numeric vector
#' @param degree The desired maximum degree
#' @param norm Logical: should polynomials be normalised to length one?
#'
#' @return A list of orthogonal polynomials
#' @export
#'
#' @examples
#' x <- c(0:3, 5)
#' P <- poly.orth(x)
#' plot(P)
#' Pf <- as.function(P)
#' zapsmall(crossprod(Pf(x)))
poly.orth <- function(x, degree = length(unique(x)) - 1, norm = TRUE) {
  at <- attr(poly(x, degree), "coefs")
  a <- at$alpha
  N <- at$norm2
  x <- polynom()
  p <- list(polynom(0), polynom(1))
  for(j in 1:degree)
    p[[j + 2]] <- (x - a[j]) * p[[j + 1]] - N[j + 1]/N[j] * p[[j]]
  p <- p[-1]
  if(norm) {
    sqrtN <- sqrt(N[-1])
    for(j in 1 + 0:degree) p[[j]] <- p[[j]]/sqrtN[j]
  }
  class(p) <- "polylist"
  p
}

.polylist_from_list <- function(x) {
  structure(lapply(x, as.polynom), class = "polylist")
}

#' @rdname polynom
#' @export
polylist <- function(...) {
  .polylist_from_list(list(...))
}

#' @rdname polynom
#' @export
is.polylist <- function(x) {
  inherits(x, "polylist")
}

#' @rdname polynom
#' @export
as.polylist <- function(x) {
  if(is.polylist(x)) x
  else if(is.list(x)) .polylist_from_list(x)
  else polylist(x)
}

#' @rdname deriv.polynom
#' @export
deriv.polylist <- function(expr, ...) {
  structure(lapply(expr, deriv), class = class(expr))
}

#' @rdname deriv.polynom
#' @export
integral.polylist <- function(expr, ...) {
  result <- lapply(expr, integral, ...)
  if (length(result) > 0 && is.polynom(result[[1]]))
    class(result) <- class(expr)
  result
}

#' @rdname plot.polynom
#' @export
plot.polylist <- function(x, xlim = 0:1, ylim = range(Px),
                          type = "l", xlab = "x", ylab = "P(x)",
                          ..., len = 1000) {
  p <- x
  if(missing(xlim)) {
    ## try to cover the "interesting" region
    xlim <- range(Re(unlist(lapply(p, summary.polynom))))
  }
  if(any(is.na(xlim))) {
    warning("summary of polynom fails. Using nominal xlim")
    xlim <- 0:1
  }
  if(diff(xlim) == 0)
    xlim <- xlim + c(-1, 1)/2
  if(length(xlim) > 2)
    x <- xlim
  else {
    eps <- diff(xlim)/100
    xlim <- xlim + c( - eps, eps)
    x <- seq(xlim[1], xlim[2], len = len)
  }
  Px <- unlist(lapply(p, function(Pi, x) Pi(x), x))
  if(!missing(ylim))
    Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
  plot(cbind(x, Px), xlab = xlab, ylab = ylab, type = "n",
       xlim = xlim, ylim = ylim, ...)
  grid(lty = "dashed")
  for(i in seq(along = p))
    lines(p[[i]], lty = i, col = i, ...)
  invisible()
}

#' Print method for polynomial objects
#'
#' @param x  A polynomial object or list thereof
#' @param ... Additional arguments passed on to methods
#'
#' @return The original object, invisibly.
#' @export
print.polylist <- function(x, ...) {
  y <- x
  x <- unclass(x)
  cat("List of polynomials:\n")
  if(length(x) > 0) {
    nam <- names(x)
    if(is.null(nam)) {
      for(i in 1:length(x)) {
        cat(paste("[[", i, "]]\n", sep=""))
        print(x[[i]], ...)
        cat("\n")
      }
    } else {
      for(n in nam) {
        cat(paste("$\"", n, "\"\n", sep=""))
        print(x[[n]], ...)
        cat("\n")
      }
    }
  } else {
    NextMethod("print", x, ...)
  }
  invisible(y)
}

#' Concatenation of polynomial objects into lists
#'
#' @param ... Polynomial or polylist objects
#' @param recursive Logical, should the concatenation flatten all component lists?
#'
#' @return A polylist object with all argumets included
#' @export
c.polynom <- function(..., recursive = FALSE) {
  .polylist_from_list(unlist(lapply(list(...), as.polylist),
                             recursive = FALSE))
}

#' @rdname c.polynom
#' @export
c.polylist <- c.polynom

#' Extract components of a list of polynomials
#'
#' @param x A polylist object
#' @param i An index vector of any crongruent form
#'
#' @return A polylist object of the components
#' @export
"[.polylist" <- function(x, i) {
  .polylist_from_list(NextMethod("["))
}


#' Component repition
#'
#' Repeat components of a polylist object
#'
#' @param x A single polynom or polylist object
#' @param times,... As for the base package function \code{rep}.
#'
#' @return The resulting polylist object.
#' @export
rep.polylist <- function(x, times, ...) {
  .polylist_from_list(NextMethod("rep"))
}

#' @rdname rep.polylist
#' @export
rep.polynom <- function(x, times, ...) {
  rep.polylist(polylist(x), times, ...)
}

#' Unique components
#'
#' Remove duplicated polynomials in a polylist object
#'
#' @param x A polylist object
#' @param incomparables Logical: as for the base function \code{unique}
#' @param ... As for the base function \code{unique}
#'
#' @return A polylist object with no duplicated components
#' @export
unique.polylist <- function(x, incomparables = FALSE, ...) {
  .polylist_from_list(NextMethod("unique"))
}

#' Chenge origin of a polynomial
#'
#' Given a polynomial P(x) and a new origin \code{o}, find
#' the polynomial Q(x) = P(x + o).  I.e. Q(0) = P(o)
#'
#' @param p A polynom or polylist object
#' @param o A single numeric quantity specifying the new x-origin
#' @param ... currently not used
#'
#' @return A polynom or polylist object with x measured from the new origin
#' @export
change.origin <- function(p, o, ...) {
  UseMethod("change.origin")
}

#' @rdname change.origin
#' @export
change.origin.default <- function(p, o, ...) {
  stop("unimplemented method")
}

#' @rdname change.origin
#' @export
change.origin.polynom <- function(p, o, ...) {
  p(polynom() + as.numeric(o)[1])
}


#' @rdname change.origin
#' @export
change.origin.polylist <- function(p, o, ...) {
  structure(lapply(p, change.origin, o = o), class = "polylist")
}

#' Plot method for polynomials
#'
#' Plot methods for polynom or polylist objects
#'
#' @param x A polynom or polylist object to be plotted
#' @param xlim,ylim as for graphics::plot
#' @param type as for graphics::plot
#' @param xlab,ylab as for graphics::plot
#' @param ... additional arguments passed on t methods
#' @param len positive integer defining the point or curve resolution
#' @param limits x-limits for the polynomial, default: the entire plot.
#'        For polylist objects this may be a two column matrix.
#' @param col,lty Colour(s) and line type(s) as for graphics::plot
#'
#' @return Nothing of interest, invisibly
#' @export
#'
#' @examples
#' p <- poly.from.zeros((-3):4)
#' plot(p)
#' lines(deriv(p), col = "red")
plot.polynom <- function(x, xlim = 0:1, ylim = range(Px),
                         type = "l", xlab = "x", ylab = "p(x)",
                         ..., len = 1000, limits = pu[1:2]) {
  p <- x
  if(missing(xlim))
    xlim <- range(c(0, Re(unlist(summary(p)))))
  if(any(is.na(xlim))) {
    warning("summary of polynom fails. Using nominal xlim")
    xlim <- 0:1
  }
  if(diff(xlim) == 0)
    xlim <- xlim + c(-1, 1)/2
  if(length(xlim) > 2)
    x <- xlim
  else {
    eps <- diff(xlim)/100
    xlim <- xlim + c(- eps, eps)
    x <- seq(xlim[1], xlim[2], len = len)
  }
  Px <- p(x)
  if(!missing(ylim))
    Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
  plot(x, Px, xlim = xlim, ylim = ylim, type = "n",
       xlab = xlab, ylab = ylab, ...)
  grid(lty = "dashed")
  pu <- par("usr")
  x <- seq(limits[1], limits[2], len = len)
  lines(x, p(x), type = type, ...)
}

#' @rdname plot.polynom
#' @export
lines.polynom <- function(x, ..., len = 1000, limits = pu[1:2])  {
  p <- x
  pu <- par("usr")
  x <- seq(limits[1], limits[2], len = len)
  lines(x, p(x), ...)
}

#' @rdname plot.polynom
#' @export
points.polynom <- function(x, ..., len = 100, limits = pu[1:2])  {
  p <- x
  pu <- par("usr")
  at <- seq(limits[1], limits[2], len = len)
  points(at, p(at), ...)
}

#' @rdname plot.polynom
#' @export
lines.polylist <- function(x, ..., len = 1000, limits = pu[1:2],
                           col = seq_along(x), lty = seq_along(x)) {
  n <- length(x)
  col <- rep_len(col, length.out = n)
  lty <- rep_len(lty, length.out = n)
  pu <- par("usr")
  if(!is.matrix(limits)) {
    limits <- rbind(limits[1:2])
  }
  limits <- cbind(rep_len(limits[,1], length.out = n),
                  rep_len(limits[,2], length.out = n))
  for(i in seq_along(x)) {
    lines(x[[i]], col = col[i], lty = lty[i], len = len,
          limits = limits[i, ],...)
  }
}

#' @rdname plot.polynom
#' @export
points.polylist <- function(x, ..., len = 100) {
  for(i in seq(along = x)) {
    points(x[[i]], pch = i, col = i, len = len, ...)
  }
}

#' Lagrange interpolation polynomial
#'
#' Calculate the Lagrange interpolation polynomial, or list of polynomials, given
#' a set of (x, y) points to fit
#'
#' @param x A numeric vector of x-points at which the y-values are specified.
#' @param y Either a numeric vector of the same length as \code{x} or a numeric
#'          matrix with rows matching the length of \code{x}.  If \code{y} is
#'          missing (not specified) then a polynomial with zero at \code{x} is
#'          returned.
#' @param tol A numeric tolerance for duplicated \code{x} values.
#' @param lab A character string vector of names for the list result when \code{y} is a matrix.
#' @param ... A list of specified zeros (for subsidiary functions)
#'
#' @return An interpolation polynomial, or list of interpolating polynomials.
#' @export
#'
#' @examples
#' (p <- poly.calc(0:5)) ## same as poly.from.zeros(0:5)
#' (p <- poly.calc(0:5, exp(0:5)))
#' plot(p)
#' curve(exp, add = TRUE, col = "red")
poly.calc <- function(x, y,
                      tol = sqrt(.Machine$double.eps),
                      lab = dimnames(y)[[2]]) {
  if(missing(y) || all(y == 0)) { ## case 1: polynomial from zeros
    p <- 1
    for(xi in x)
      p <- c(0, p) - c(xi * p, 0)
    return(polynom(p))
  }                ## case 2: Lagrange interpolating polynomial
  if(is.matrix(y)) {
    if(length(x) != nrow(y))
      stop("x and y are inconsistent in size")
    lis <- list()
    if(is.null(lab))
      lab <- paste("p", 1:(dim(y)[2]), sep = "")
    for(i in 1:dim(y)[2])
      lis[[lab[i]]] <- Recall(x, y[, i], tol)
    return(structure(lis, class = "polylist"))
  }
  if(any(toss <- duplicated(x))) {
    crit <- max(tapply(y, x, function(x) diff(range(x))))
    if(crit > tol)
      warning("some duplicated x-points have inconsistent y-values")
    keep <- !toss
    y <- y[keep]
    x <- x[keep]
  }
  if((m <- length(x)) != length(y))
    stop("x and y(x) do not match in length!")
  if(m <= 1)
    return(polynom(y))
  r <- 0
  for(i in 1:m)
    r <- r + (y[i] * coef(Recall(x[ - i])))/prod(x[i] - x[ - i])
  r[abs(r) < tol] <- 0
  polynom(r)
}

#' @rdname poly.calc
#' @export
poly.from.zeros <- function(...) {
  poly.calc(unlist(list(...)))
}

#' @rdname poly.calc
#' @export
poly.from.roots <- poly.from.zeros

#' @rdname poly.calc
#' @export
poly.from.values <- poly.calc

#' Evaluate a polynomial
#'
#' Evaluate a polynomial, or polylist object components.
#'
#' @param object A polynomial or polylist object
#' @param newdata A target object at which to evaluate.
#' @param ... Not used
#'
#' @return If \code{newdata} is a numeric vector, a numeric vector of
#'         results.  If \code{newdata} is a polynomial, then the composition
#'         is returned as a polynomial, or polylist object.
#' @export
predict.polynom <- function(object, newdata, ...) {
  object(newdata)
}

#' @rdname predict.polynom
#' @export
predict.polylist <- function(object, newdata, ...) {
  sapply(object, function(x, .n) x(.n), .n = newdata)
}

#' Find Polynomial Zeros
#'
#' Solve polynomial equations, a(x) = b(x), or alternatively
#' find the zeros of the polynomial a(x) - b(x)
#'
#' @param a,b Polynomials for the LHS and RHS respectively
#' @param ... Currently unused
#'
#' @return A vector of roots, usually complex
#' @export
#'
#' @examples
#' p <- poly.calc(0:5)
#' solve(p)
#' solve(p, 1)
solve.polynom <- function(a, b, ...) {
  if(!missing(b))
    a <- a - b
  a <- coef(a)
  if(a[1] == 0) {
    z <- rle(a)$lengths[1]
    a <- a[-(1:z)]
    r <- rep(0, z)
  }
  else
    r <- numeric(0)
  switch(as.character(length(a)),
         "0" =,
         "1" = r,
         "2" = sort(c(r,  - a[1]/a[2])),
         {
           a <- rev(a)
           a <- (a/a[1])[-1]
           M <- rbind( - a, cbind(diag(length(a) - 1), 0))
           sort(c(r, eigen(M, symmetric = FALSE,
                           only.values = TRUE)$values))
         })
}

#' @rdname solve.polynom
#' @export
solve.polylist <- function(a, b, ...) {
  if(!missing(b)) lapply(a, solve.polynom, b) else
    lapply(a, solve.polynom)
}

#' Polynomial summary
#'
#' Provide a succinct summary of the critical points of a polynomial, or list thereof
#'
#' @param object,x A polynomial or polylist object
#' @param ...  Currently unused
#'
#' @return A list giving the zeros, stationary points and points of inflexion of the polynomial(s)
#' @export
#'
#' @examples
#' p <- poly.calc(0:5)
#' summary(p)
summary.polynom <- function(object, ...) {
  dp <- deriv(object)
  structure(list(zeros = solve(object),
                 stationaryPoints = solve(dp),
                 inflexionPoints = solve(deriv(dp))),
            class = "summary.polynom",
            originalPolynomial = object)
}

#' @rdname summary.polynom
#' @export
summary.polylist <- function(object, ...) {
  lapply(object, summary.polynom)
}

#' @rdname summary.polynom
#' @export
print.summary.polynom <- function(x, ...) {
  cat("\n Summary information for:\n")
  print(attr(x, "originalPolynomial"))
  cat("\n Zeros:\n")
  print(x$zeros)
  cat("\n Stationary points:\n")
  print(x$stationaryPoints)
  cat("\n Points of inflexion:\n")
  print(x$inflexionPoints)
  invisible(x)
}

.monic <- function(p) {
  a <- coef(p)
  polynom(a/a[length(a)])
}

.degree <- function(x) {
  length(coef(x)) - 1
}


.effectively_zero <- function (p, tolerance = .Machine$double.eps^0.5) {
  all(abs(coef(p)) < tolerance)
}


.GCD2 <- function(x, y) {
  if(.effectively_zero(y)) x
  else if(.degree(y) == 0) polynom(1)
  else Recall(y, x %% y)
}

.LCM2 <- function(x, y) {
  if(.effectively_zero(x) || .effectively_zero(y))
    return(polynom(0))
  (x / .GCD2(x, y)) * y
}

#' Greatest common divisor
#'
#' Find a monic polynomial of maximal degree that divides each
#' of a set of polynomials exactly
#'
#' @param ... A list of polynomials or polylist objects
#'
#' @return A polynomial giving the greatest common divisor, as defined above
#' @export
#'
#' @examples
#' p <- poly.calc(0:5)
#' r <- poly.calc(1:6)
#' greatest_common_divisor(p, r)
#' solve(greatest_common_divisor(p, r))
#' lowest_common_multiple(p, r)
#' solve(lowest_common_multiple(p, r))
GCD <- function(...) {
  UseMethod("GCD")
}

#' @rdname GCD
#' @export
greatest_common_divisor <- function(...) {
  UseMethod("GCD")
}


#' @rdname GCD
#' @export
GCD.polynom <- function(...) {
  args <- c.polylist(...)
  if(length(args) < 2)
    stop("Need at least two polynoms.")
  .monic(.accumulate(.GCD2, args[[1]], args[-1], FALSE))
}

#' @rdname GCD
#' @export
GCD.polylist <- GCD.polynom

#' Lowest Common Multiple
#'
#' For a list of polynomials, find the lowest degree monic
#' polynomial into which each divides exactly
#'
#' @param ... A list of polynomials or polylist objects
#'
#' @return A polynomial giving the lowest common multiple
#' @export
#'
#' @examples
#' p <- poly.calc(0:5)
#' r <- poly.calc(1:6)
#' greatest_common_divisor(p, r)
#' solve(greatest_common_divisor(p, r))
#' lowest_common_multiple(p, r)
#' solve(lowest_common_multiple(p, r))
LCM <- function(...) {
  UseMethod("LCM")
}

#' @rdname LCM
#' @export
lowest_common_multiple <- function(...) {
  UseMethod("LCM")
}

#' @rdname LCM
#' @export
LCM.polynom <- function(...) {
  args <- c.polylist(...)
  if(length(args) < 2)
    stop("Need at least two polynoms.")
  .monic(.accumulate(.LCM2, args[[1]], args[-1], FALSE))
}

#' @rdname LCM
#' @export
LCM.polylist <- LCM.polynom
