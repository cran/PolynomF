## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "",
                      fig.height = 4.5,
                      fig.width = 6)
if(knitr::is_latex_output()) {
  knitr::opts_chunk$set(fig.align = "center",
                        out.height = "30%")
}
library(PolynomF)
library(knitr)
setHook("plot.new",
        list(las = function() par(las = 1),
             pch = function() par(pch = 20)),
        "append")
curve <- function(..., add = FALSE, n = 1001, type = "l") {
  if(add) {
    graphics::curve(..., add = TRUE, n = n)
  } else {
    graphics::curve(..., add = FALSE, type = "n")
    grid(lty = "dashed")
    graphics::curve(..., add = TRUE, n = 1001, type = type)
  }
}
# setHook("plot.new",
#         list(las = function() par(las = 1),
#              pch = function() par(pch = 20)),
#         "replace")


## ----fig.width = 8------------------------------------------------------------
Discrete <- function(p, q = p, x, w = function(x, ...) 1, ...) sum(w(x, ...)*p(x)*q(x))
PC <- poly_orth_general(inner_product = Discrete, degree = 4, 
                        x = 0:100, w = dpois, lambda = 1)
plot(PC, lwd = 2, legend = TRUE)
title(main = "Poisson-Charlier(1) polynomials to degree 4", font.main = 3)

## -----------------------------------------------------------------------------
(p1 <- poly_calc(1:6))       ## a monic polynomial with given zeros
solve(p1)                    ## check that the zeros are as specified
polyroot(coef(p1))           ## check using the Traub-Jenkins algorithm
p2 <- -change_origin(p1, 3)  ## a negative shifted version of p1
plot(p1, xlim = c(-2.5, 6.5), ylim = c(-20, 20), lwd = 2, col = "steel blue", bty = "n")
lines(p2, col = "brown", lwd = 2)
abline(h = 0)
points(cbind(solve(p1), 0), pch = 1,  cex = 1,   col = "brown")
points(cbind(solve(p2), 0), pch = 19, cex = 0.5, col = "steel blue")
## stationary points
stat <- solve(deriv(p1))
lines(tangent(p1, stat), limits = cbind(stat-0.5, stat + 0.5),
      lty = "solid", col = "black")
points(stat, p1(stat), pch = 19, cex = 0.5, col = "red")
stat <- solve(deriv(p2))
lines(tangent(p2, stat), limits = cbind(stat-0.5, stat + 0.5),
      lty = "solid", col = "cadet blue")
points(stat, p2(stat), pch = 19, cex = 0.5, col = "red")
## various checks
z <- (-2):6
setNames(p1(z), paste0("z=", z))
setNames(p2(z), paste0("z=", z))
setNames((p1*p2)(z), paste0("z=", z))
p3 <- (p1 - 2 * p2)^2                         ## moderately complicated expression.
setNames(p3(0:4), paste0("z=", 0:4))          ## should have zeros at 1, 2, 3

## -----------------------------------------------------------------------------
x0 <- c(0:3, 5)
op <- poly_orth(x0, norm = TRUE)
plot(op, lwd = 2, legend = TRUE)
fop <- as.function(op)        ## Explicit coercion needed for polylist
zap(crossprod(fop(x0)))       ## Verify orthonormality

## -----------------------------------------------------------------------------
x <- polynomial()
Tr <- polylist(1, x)
for(j in 3:15) {
  Tr[[j]] <- 2*x*Tr[[j-1]] - Tr[[j-2]]
}
Tr <- setNames(Tr, paste0("T", sub(" ", "_", format(seq_along(Tr)-1))))

## -----------------------------------------------------------------------------
ChebyT <- function(p, q = p) {
  integrate(function(x) 1/sqrt(1-x^2)*p(x)*q(x), lower = -1, upper = 1,
            subdivisions = 500, rel.tol = .Machine$double.eps^0.5)$value
}
zap(outer(Tr, Tr, Vectorize(ChebyT))*2/pi) ## check of orthogonality

## -----------------------------------------------------------------------------
fx <- function(x) dnorm(x, sd = 0.25)
b <- sapply(Tr, ChebyT, q = fx)/sapply(Tr, ChebyT)
fx_approx <- sum(b * Tr)
curve(fx, xlim = c(-1,1))
lines(fx_approx, col = "red", limits = c(-1, 1))
curve(1e4*(fx(x) - fx_approx(x)), xlim = c(-1,1)) ## error pattern x 10000


## -----------------------------------------------------------------------------
Fx <- function(x) pnorm(x, sd = 0.25) 
Fx_approx <- integral(fx_approx) + 0.5
curve(Fx, xlim = c(-1, 1))
lines(Fx_approx, col = "red", limits = c(-1,1))
curve(1e4*(Fx(x) - Fx_approx(x)), xlim = c(-1,1)) ## error pattern x 10000


## -----------------------------------------------------------------------------
x <- polynomial()
Ur <- polylist(1, 2*x)

for(j in 3:15) {
  Ur[[j]] <- 2*x*Ur[[j-1]] - Ur[[j-2]]
}
Ur <- setNames(Ur, paste0("U", sub(" ", "_", format(seq_along(Ur)-1))))
plot(Ur, ylim = c(-2,2))
ChebyU <- function(p, q = p) {
  integrate(function(x) sqrt(1-x^2)*p(x)*q(x), lower = -1, upper = 1,
            subdivisions = 500, rel.tol = .Machine$double.eps^0.5)$value
}
b <- sapply(Ur, ChebyU, q = asin)/(pi/2)
asin_approx <- sum(b * Ur)

curve(asin, xlim = c(-1,1))
lines(asin_approx, col = "red", limits = c(-1,1))

curve(1e4*(asin(x) - asin_approx(x)), xlim = c(-1,1), ylim = c(-50,50)) ## errors by 10000


## -----------------------------------------------------------------------------
P <- polynomial(c(1, 5, 3, 1))/10
s <- polynomial(c(0, 1))
(mean_offspring <- deriv(P)(1))
pretty_poly <- bquote(italic(P)(italic(s)) == ' '*
                      .(parse(text = gsub("x", "italic(' '*s)", as.character(P)))[[1]]))
plot(s, xlim = c(0,1), ylim = c(0,1), bty = "n", type = "n", main = pretty_poly,
     xlab = expression(italic(s)), ylab = expression(italic(P)(italic(s))))
x <- c(0,1,1,0)
y <- c(0,0,1,1)
segments(x, y, 1-y, x, lty = "solid", lwd = 0.2)
lines(s, limits = 0:1, col = "grey")
lines(P, limits = 0:1)
lines(tangent(P, 1), col = "red", limits = c(0.5, 1.0), lwd = 1.5)
ep <- solve((P - s)/(1 - s)) ## extinction; factor our the known zero at s = 1
(sg <- ep[ep > 0 & ep < 1])  ## extract the appropriate value (may be complex, in general)
plot(s, xlim = c(0,1), ylim = c(0,1), type = "n", bty = "n", main = pretty_poly,
     xlab = expression(italic(s)), ylab = expression(italic(P)(italic(s))))
segments(x, y, 1-y, x, lty = "solid", lwd = 0.2)
lines(s,     col = "grey", limits = 0:1)  ## higher generations
lines(P,          col = 1, limits = 0:1)
lines(P(P),       col = 2, limits = 0:1)
lines(P(P(P)),    col = 3, limits = 0:1)
lines(P(P(P(P))), col = 4, limits = 0:1)
arrows(sg, P(sg), sg, par("usr")[3], angle = 15, length = 0.125)
text(sg, 0.025, bquote(sigma == .(round(sg, 4))), pos = 4, adj = 0, cex = 0.7)

## -----------------------------------------------------------------------------
x0 <- 80:89
y0 <- c(487, 370, 361, 313, 246, 234, 173, 128, 88, 83)
p <- poly_calc(x0, y0)          ## leads to catastropic numerical failure!
range(p(x0) - y0)               ## these should be "close to zero"!
##
## Try another algorithm for finding the Lagrange interpolating polynomial
##
pn <- neville(x0, y0)                             ## Neville's algorithm
range(pn(x0) - y0)                 
pl <- lagrange(x0, y0)                            ## lagrange's formula
range(pl(x0) - y0)
##
## Try by relocating the x-values for computational purposes.
##
pr <- local({
  .p <- poly_calc(x0 - 84, y0)  ## changing origin fixes the problem
  function(x) {
    .p(x - 84)
  }
})
range(pr(x0) - y0)       ## these are 'close to zero'.
plot(p, xlim = c(80, 89), ylim = c(-100, 800), bty = "n")    ## this looks ugly!
lines(pn, col = 2, xpd = NA)                                 ## Just as bad!
lines(pl, col = 3, xpd = NA)                                 ## but different.
curve(pr, add = TRUE, lwd = 2, col = 4)                      ## Success at last!
points(x0, y0, col = 1)
legend("top", legend = c("poly_calc", "neville", "lagrange", "relocation"),
       lty = "solid", col = 1:4, cex = 0.7, title = "Method", ncol = 2)

