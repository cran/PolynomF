## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "",
                      fig.height = 5.5,
                      fig.width = 7,
                      fig.align = "center",
                      out.height = "0.3\\textheight")
library(PolynomF)
setHook("plot.new",
        list(las = function() par(las = 1),
             pch = function() par(pch = 20)),
        "append")

## ------------------------------------------------------------------------
He <- polylist(polynom(1), polynom(0:1))
x <- polynom()
for (n in 3:6) {
  He[[n]] <- x * He[[n-1]] - (n-2) * He[[n-2]] ## R indices start from 1, not 0
}
plot(He)
plot(deriv(He))
plot(integral(He))

## ------------------------------------------------------------------------
x <- c(0, 1, 2, 3, 5)
(op <- poly.orth(x))
(fop <- as.function(op))
zapsmall(crossprod(fop(x)))  ## Verify orthonormality

## ------------------------------------------------------------------------
(p1 <- poly.calc(1:6))
(p2 <- change.origin(p1, 3))
p1(0:7)
p2(0:7)
p2(0:7 - 3)
(p3 <- (p1 - 2 * p2)^2)         # moderate arithmetic expression.
p3(0:4)                         # should have 1, 2, 3 as zeros

## ------------------------------------------------------------------------
P <- polynom(c(1, 5, 2, 2)/10)
s <- polynom(c(0, 1))
(mean_offspring <- deriv(P)(1))
plot(P, xlim = c(0,1), ylim = c(0,1), xlab = "s", ylab = "P(s)")
abline(h=0:1, v=0:1, lty = "dotted")
lines(tangent(P, 1), lty = "dashed", col = "red")  
lines(s, lty = "dotted", col = "grey")
## higher generations
plot(P, xlim = c(0,1), ylim = c(0,1), xlab = "s", ylab = "P(s)")
lines(P(P), col = 2)
lines(P(P(P)), col = 3)
lines(P(P(P(P))), col = 4)
solve(P, s)                ## for the extinction probability
(ep <- solve((P-s)/(s-1))) ## factor our the known zero at s = 1
ex <- ep[2]                ## extract the appropriate value
segments(ex, 0, ex, P(ex), col = "red", lwd = 2)
abline(h=0:1, v=0:1, lty = "dotted")

## ------------------------------------------------------------------------
x <- 80:89
y <- c(487, 370, 361, 313, 246, 234, 173, 128, 88, 83)

p <- poly.calc(x, y)        ## leads to catastropic numerical failure!
p(x) - y                    ## these should be "close to zero"!

p1 <- poly.calc(x - 84, y)  ## changing origin fixes the problem
p1(x - 84) - y              ## these are 'close to zero'.

plot(p1, xlim = c(80, 89) - 84, xlab = "x - 84")
points(x - 84, y, col = "red", cex = 2)

#### Can we now write the polynomial in "raw" form?
z <- polynom()
p0 <- p1(z - 84) ## attempting to change the origin back to zero
                 ## leads to severe numerical problems again
plot(p0, xlim = c(80, 89))
points(x, y, col = "red", cex = 2) ## major numerical errors due to finite precision

