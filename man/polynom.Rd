\name{polynom}
\alias{polynom}
\alias{as.polynom}
\alias{is.polynom}
\alias{polylist}
\alias{as.polylist}
\alias{is.polylist}
\alias{print.polynom}
\alias{print.polylist}

\title{ Constructor, coercion, predicate and print functions for polynom
  and polylist objects}
\description{
  The function \code{polynom} is the constructor function for objects of
  the eponymous S3 class.  \code{as.polynom} is the standad coercion
  function to the same class and \code{is.polynom} tests for inheritance
  from the class.  \code{polylist} objects are lists of \code{polynom}
  objects, again with the eponymous class.
}
\usage{
polynom(a = c(0, 1), \dots, eps = 0)
as.polynom(a)
is.polynom(a)
polylist(\dots)
as.polylist(x)
is.polylist(x)
\method{print}{polynom}(x, variable = "x",
   digits = getOption("digits"), decreasing = FALSE, \dots)
\method{print}{polylist}(x, \dots) 
}

\arguments{
  \item{a}{
    Numeric coefficient vector for the constructor, specifying
    the coefficients of the powers from 0 to the maximum, without gaps,
    in that order, that is, in power series order.  For the coercion
    and predicate functions this may be a coefficient vector, or an
    existing \code{polynom} object.}
  \item{\dots}{As yet, unused for \code{polynom} but may be used in
    future releases. With \code{polylist}, objects which are, or may be
    coerced to, \code{polynom} objects.  With \code{print} methods,
    additional agruments either for the generic or for the
    \code{print.polynom} method.}
  \item{x}{Object which may be coerced to, or tested for, class
    \code{polylist}.} 
  \item{eps}{Absolute tolerance below which components are considered
    zero.}
  \item{variable}{Character string giving name to be used for the
    independent variable.}
  \item{digits}{Integer giving the number of significant digits to use
    in the printed representation of the coefficients.}
  \item{decreasing}{Logical value.  Should the powers appear in
    decreasing order?}
}
\details{
  \code{polynom} objects are represented as \code{R} functions with the
  coefficient vector held in an enclosing environment, or closure.  They
  may be used as \code{R} functions, or as mathematical objects which behave
  under arithmetic and calculus operatons like polynomials.
  
  Note that polynomials with complex coefficients are not (yet)
  supported.  Nevertheless the root-finding methods described elsewhere
  will usually give
  complex zeros.

  \code{polylist} objects are lists of \code{polynom} objects.  Methods
  are available for dealing with all polynomials on the list simultaneously.}

\value{ For \code{polynom} and \code{as.polynom}, an \code{R} function
  with class \code{polynom}.  For \code{is.polynom}, a logical value}
\references{None }
\author{Bill Venables }

\examples{
x <- polynom()  # polynomial 'x'
p <- (x-1)^2 + 10*x^3 + 5*x^4
p
# 1 - 2*x + x^2 + 10*x^3 + 5*x^4 

plot(polylist(p, deriv(p), integral(p)),
    xlim = c(-2, 1), ylim = c(-10, 10))
abline(h=0, lty = "dashed", col = "grey")

print(p, variable = "z", decreasing = TRUE)

}

\keyword{symbolmath}

