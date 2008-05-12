\name{PolynomF-package}
\alias{PolynomF-package}
\alias{PolynomF}
\docType{package}
\title{Univariate polynomial objects in R}
\description{
  A package to implement a class of objects that behave like univariate
  polynomails.  Arithmetic operations (addition, subtraction,
  multiplication, division, remainder, raising to a non-negative integer
  power) are supported in a natural way.  The objects also act as R
  functions.  This package is a successor to the `polynom' package,
  but has a simpler and more convenient representation for the objects.
  Like `polynom' it uses S3 classes and methods.}
\details{
  \tabular{ll}{
    Package: \tab PolynomF\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2008-05-05\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes\cr
  }

  The constructor function \code{polynom} is used to create polynomial
  objects from their coefficient vector, in power series order.  Once
  polynomials are constructed they may used as objects in arithmetic
  operations, integration and differentiation, and as R functions that
  evaluate the polynomial either at a numeric or complex vector, or at
  another polynomial, i.e. substituting one polynomial into
  another. Facilities are also provided for graphical presentation and
  calculation of complex zeros.

  The constructor function \code{polylist} may be used to create a
  \emph{list} of polynomial objects.  Operations on \code{polylist}
  objects include simultaneous graphical display of all components and
  coercion to function.  The function may then be used to evaluate all
  all polynomials on the list simultaneously at the same argument.

}
\author{
  Bill Venables, with some code inherited from the original \code{package}
  by Bill Venables and Kurt Hornik.

  Maintainer:  <Bill.Venables@gmail.com>
}
\references{
  None.
}
\keyword{symbolmath}

\examples{
x <- polynom()
p <- (x-1)^2 + 1
p
plot(p)

pv <- p(-3:4); pv

p1 <- p(p-1); p1;
plot(polylist(p, p1))

## Hermite polynomials to degree 10
H <- polylist(1, x)
for(n in 2:10)
    H[[n+1]] <- x*H[[n]] - (n-1)*H[[n-1]]
H
### normalisation to unit length
for(n in 1:11)
    H[[n]] <- H[[n]]*exp(-lgamma(n)/2)

plot(H, xlim = c(-3,3))

## orthogonality relationship check:
f <- function(i,j) stats::integrate(function(z)
      dnorm(z)*H[[i+1]](z)*H[[j+1]](z), -Inf, Inf)

f(2,3)
f(4,4)
}
