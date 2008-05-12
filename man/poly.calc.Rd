\name{poly.calc}
\alias{poly.calc}
\alias{poly.from.values}
\alias{poly.from.zeros}
\alias{poly.from.roots}
\alias{poly.orth}

\title{Functions to generate polynomials in several standard ways }
\description{
  \code{poly.calc} (alias \code{poly.from.values}) computes the Lagrange
  interpolating polynomial.  \code{poly.from.zeros} (alias
  \code{poly.from.roots}) computes the monic polynomial with specified
  zeros. \code{poly.orth} calculates polynomials orthogonal over a
  discrete set of $x-$values, as done numerically by the standard \code{R}
  function \code{poly}. 
}
\usage{
poly.calc(x, y, tol = sqrt(.Machine$double.eps), lab = dimnames(y)[[2]])
poly.from.values(x, y, tol = sqrt(.Machine$double.eps), lab = dimnames(y)[[2]])
poly.from.zeros(\dots)
poly.from.roots(\dots)
poly.orth(x, degree = length(unique(x)) - 1, norm = TRUE)
}

\arguments{
  \item{x}{A numeric vector of values for the polynomial variable. }
  \item{y}{A numeric vector or matrix specifying values for the polynomial.}
  \item{tol}{A numeric tolerance }
  \item{lab}{A character vector providing names for the \code{polylist}
    of polynomials.}
  \item{degree}{The maximum degree for the orthogonal polynomials
    required.}
  \item{norm}{Logical value.  Should the polynomials be normalised to be
    of length 1?)}
  \item{\dots}{Not presently used.}
}
\details{
  Given a vector of distinct values \code{x} and a vector \code{y} of
  the same length, \code{poly.calc} computes the Lagranging
  intrepolating polynomial they define.  If \code{y} is a matrix, its
  row size must match the length of \code{x} and interpolating
  polynomials are computed for all columns.  In this case the value is a
  \code{polylist} object.

  \code{poly.from.values} is a complete alias for \code{poly.calc}.

  The function \code{poly.from.zeros} computes the monic polynomial with
  zeros as given by the arguments.  The zeros may be specified either as
  separate artuments or as a single numeric vector.

  \code{poly.from.roots} is a complete alias for \code{poly.from.zeros}.

  \code{poly.orth} calculates polynomials orthogonal with respect to the
  uniform measure over a discrete set of $x-$values given by the
  artument \code{x}.  These are the polynomials for which the standard
  function \code{poly} can be used to compute numerical values.
}
\value{
  A \code{polynom} object, or, in the case of \code{poly.calc} and
  \code{poly.orth}, possibly a \code{polylist} object
}
\references{None }
\author{Bill Venables}




\seealso{\code{poly}}
\examples{
x <- polynom()
H <- polylist(x^0, x)
for(j in 2:10)
  H[[j+1]] <- x*H[[j]] - (j-1)*H[[j-1]]
Hf <- as.function(H)
x0 <- -5:5
y0 <- Hf(x0)
J <- poly.from.values(x0, y0)
all.equal(H[[11]], J[[11]])

p1 <- poly.from.zeros(-3:2)
p2 <- poly.from.zeros(0:4)
p3 <- GCD(p1, p2)
p4 <- LCM(p1, p2)

solve(polylist(p1, p2, p3, p4))

po <- poly.orth(-4:4, degree = 4)
plot(po)

round(crossprod(as.function(po)(-4:4)), 10)

}


\keyword{symbolmath }

