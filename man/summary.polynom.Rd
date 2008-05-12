\name{summary.polynom}
\alias{summary.polynom}
\alias{summary.polylist}
\alias{print.summary.polynom}
\alias{coef.polynom}
\alias{coef.polylist}
\alias{predict.polynom}
\alias{predict.polylist}

\title{Summary, coefficient and prediction methos for polynomial objects.}

\description{
  These function implement methods for \code{summary}, \code{coef}, and
  \code{predict} generic funcitons for \code{polynom} and \code{polylist}
  objects. 
}

\usage{
\method{summary}{polynom}(object, \dots)
\method{summary}{polynom}(object, \dots)
\method{print}{summary.polynom}(x, \dots)
\method{coef}{polynom}(object, \dots)
\method{coef}{polylist}(object, \dots)
\method{predict}{polynom}(object, newdata, \dots)
\method{predict}{polylist}(object, newdata, \dots)
}

\arguments{
  \item{object}{An object of class \code{polynom} or \code{polylist}. }
  \item{x}{A \code{summary.polynom} object to be printed.}
  \item{newdata}{A \code{numeric} or \code{polynom} object.}
  \item{\dots}{Not currently used. }
}

\details{
  The summary method for \code{polynom} objects provides information on
  the zeros, stationary points and points of inflexion for the object.
  For \code{polylist} objects this information is provided for each
  polynomial in the list.  The result is a list.
}
\value{
  For \code{summary.polynom} and \code{summary.polylist} a list of
  numeric vectors.

  For \code{predict.polynom} and \code{predict.polylist} an object of
  the same class as the input argument: either numeric or
  \code{polynom}.

  For \code{coef.polynom} and \code{coef.polylist} a numeric vector or
  matrix, or a list of numeric vectors.  
}
\references{None. }
\author{Bill Venables}

\examples{
x <- polynom()
L <- polylist(1, 1-x)
for(j in 2:10)
  L[[j+1]] <- (2*j - 1 - x)*L[[j]] - (j-1)^2*L[[j-1]]

summary(L[[5]])
predict(L[[5]], x-1)
L[[5]](x-1)

coef(L)

}
\keyword{symbolmath}

