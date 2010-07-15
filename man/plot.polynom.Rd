\name{plot.polynom}
\alias{plot.polynom}
\alias{plot.polylist}
\alias{lines.polylist}
\alias{lines.polynom}
\alias{points.polylist}
\alias{points.polynom}

\title{Graphical display of polynomial objects}
\description{
  Standard displsy methods for polynomial or lists of polynomials.
}
\usage{
\method{plot}{polynom}(x, xlim = 0:1, ylim = range(Px), type = "l",
   xlab = "x", ylab = "p(x)", \dots, len = 1000)
\method{plot}{polylist}(x, xlim = 0:1, ylim = range(Px), type = "l",
   xlab = "x", ylab = "P(x)", \dots, len = 1000)
\method{lines}{polylist}(x, ..., len = 1000)
\method{lines}{polynom}(x, ..., len = 1000)
\method{points}{polylist}(x, ..., len = 100)
\method{points}{polynom}(x, ..., at = seq(pu[1], pu[2], len = len), len = 100) 
}

\arguments{
  \item{x}{A \code{polynom} or \code{polylist} object. }
  \item{xlim, ylim, xlab, ylab, type}{As for \code{plot}. }
  \item{\dots}{Additional arguments sent to \code{plot}, \code{points}
or \code{lines} }
  \item{len}{The number of linear line segments to use to present
    the polynomial curve. }
  \item{at}{the $x-$values where the points are to appear.}
}
\details{
  \code{plot.polynom} will by default choose \code{x} limits to cover
  the (real parts of) the zeros, stationary points and points of
  inflexion of the ploynomial being plotted.  \code{plot.polylist}
  chooses by default an x region to accommodate all polynomials on the
  list in this way.  The current palette of colours is used for
  different components.  \code{lines.polynom} may be used to add
  individual polynomials to the plot.  The argument \code{len} may be
  used to increase or decrease the number of straight line segments used
  to represent the curves.
}
\value{
  Nothing of use.
}
\references{Nont }
\author{Bill Venables, with contributions by Kurt Hornik.}




\seealso{\code{curve} }
\examples{
x <- polynom()
L <- polylist(1, 1-x)
for(j in 2:10) L[[j+1]] <- (2*j - 1 - x)*L[[j]] - (j-1)^2*L[[j-1]]
plot(L[1:5], xlim = c(0,5), xaxs = "r", ylab = expression(L[j](z)),
xlab = "z", main = "Laguerre polynomials to degree 4")

lines(L[[6]], col = "grey", lwd = 2)
}


\keyword{symbolmath}

