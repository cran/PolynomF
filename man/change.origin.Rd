\name{change.origin}
\alias{change.origin}
\alias{change.origin.default}
\alias{change.origin.polynom}
\alias{change.origin.polylist}

\title{Shift oritin}
\description{
  These small convenience functions shift polynomials by relocating the
  origin to a new position.  The primary function is generic and active
  methods are provided for \code{polynom} and \code{polylist} objects.
}
\usage{
\method{change.origin}{default}(p, o, \dots)
\method{change.origin}{polynom}(p, o, \dots)
\method{change.origin}{polylist}(p, o, \dots)

}
\arguments{
  \item{p}{A \code{polynim} or \code{polylist} object. }
  \item{o}{The value of the original variable to become the new origin,
    that is the zero value in the shifted variable. }
  \item{\dots}{Not used }
}
\details{
  The function \code{change.origin.default} is a trap for unimplemeted
  methods.
}
\value{
  A \code{polynom} or \code{polylist} object with the origin shifted to \code{0}.
}
\references{ None }
\author{Bill Venables }

\examples{
x <- polynom()
p <- 1 - 2*x + x^2
change.origin(p, 1)
## x^2
}

\keyword{symbolmath }

