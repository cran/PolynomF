\name{solve.polynom}
\alias{solve.polynom}
\alias{solve.polylist}
\title{Find roots of polynomial equations}
\description{
  Method functions for the \code{solve} generic function to find the
  roots of a single polynomial equation or the roots of a list of
  polynomial equations
}
\usage{
\method{solve}{polynom}(a, b, \dots)
\method{solve}{polylist}(a, b, \dots)
}

\arguments{
  \item{a}{A \code{polynom} or \code{polylist} object }
  \item{b}{A \code{polynom} or scalar numeric object.  The polynomial
    equation(s) to be solved is \code{a - b = 0}.}
  \item{\dots}{Not used.}
}
\details{
  The roots of the equation will in general be complex.  In the case of
  \code{polylist} the same equation, \code{a - b = 0}, is solved for all
  polynomials in the list.  In this case the result is a list of root
  vectors.
}
\value{
  A vector of roots, or a list of such vectors.
}
\references{None}
\author{Bill Venables }

\examples{
x <- polynom()
H <- polylist(x^0, x)
for(n in 2:10)
    H[[n+1]] <- x * H[[n]] - (n-1)*H[[n-1]]
solve(H)
solve(deriv(H))

}


\keyword{symbolmath }

