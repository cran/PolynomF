\name{integral}
\alias{integral}
\alias{integral.polynom}
\alias{integral.polylist}
\alias{deriv.polynom}
\alias{deriv.polylist}

\title{Integral and differential calculus on polynomials or lists of
  polynomials.}
\description{
  Performs calculus operations on polynomial objects: differentiation
  and indefinite or definite integration.
}
\usage{
\method{integral}{polynom}(expr, limits = NULL, \dots)
\method{integral}{polylist}(expr, \dots)
\method{deriv}{polynom}(expr, \dots)
\method{deriv}{polylist}(expr, \dots)
}

\arguments{
  \item{expr}{A \code{polynom} or \code{polylist} object. }
  \item{limits}{Either \code{NULL}, implying an indefinite integral, or
    a numeric (or complex) vector of length two, implying a definite
    integral between two fixed real of complex limits.}
  \item{\dots}{ Not used, except in the case of \code{integral.polylist}
    where it may pass on common \code{limits} of integration to all
    integrals in the list.}
}
\details{
  \code{integral} is a local generic function.  The generic function
  \code{deriv} is already defined in the \code{stats} package.
}
\value{
  A \code{polynom} or numeric object giving the result of the calculus
  operation, or list of such results.  If the result is \code{polynom}
  the list is a \code{polylist}.
}
\references{None}





\examples{
p <- poly.from.zeros(-2:5)
ip <- integral(p)
ipv <- integral(p, limits = c(-2, 5))

plot(polylist(p, deriv(p)))

x <- polynom()
H <- polylist(1, x)
for(n in 2:10)
    H[[n+1]] <- x * H[[n]] - (n-1)*H[[n-1]]

solve(deriv(H))

}

\keyword{symbolmath }

