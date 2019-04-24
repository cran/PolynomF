## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "",
                      fig.height = 5.5,
                      fig.width = 7,
                      fig.align = "center",
                      out.height = "0.3\\textheight")
library(PolynomF)
library(knitr)
setHook("plot.new",
        list(las = function() par(las = 1),
             pch = function() par(pch = 16)),
        "append")

## ------------------------------------------------------------------------
Discrete <- function(p, q = p, x, w = function(x) rep(1, length(x))) {
  sum(w(x)*p(x)*q(x))
}
poly_orth_general(inner_product = Discrete, degree = 5, 
                  x = 0:20, w = function(x) dpois(x, 1))

## ---- echo=FALSE---------------------------------------------------------
old <- ls("package:PolynomF", pattern = "\\.")
new <- gsub("\\.", "_", old)
tab <- data.frame(`Old name` = old, `New name` = new, check.names = FALSE)
kable(tab, caption = "Function name changes in version 2.0.0")

## ------------------------------------------------------------------------
(p1 <- poly_calc(1:6))
(p2 <- change_origin(p1, 3))
p1(0:7)
p2(0:7)
p2(0:7 - 3)
(p3 <- (p1 - 2 * p2)^2)         # moderate arithmetic expression.
p3(0:4)                         # should have 1, 2, 3 as zeros

## ------------------------------------------------------------------------
x0 <- c(0, 1, 2, 3, 5)
(op <- poly_orth(x0, norm = TRUE))
fop <- as.function(op)       # Explicit coercion needed for polylist
zapsmall(crossprod(fop(x0)))  # Verify orthonormality

## ---- out.height="0.26\\textheight"--------------------------------------
x <- polynomial()
He <- polylist(polynomial(1), x)
for (n in 3:6) {
  He[[n]] <- x * He[[n-1]] - (n-2) * He[[n-2]] # R indices start from 1, not 0
}
(He <- setNames(He, paste0("He", seq_along(He) - 1)))
He0 <- poly_orth_general(Hermite, 5)  # using brute computation
sapply(He-He0, function(p) max(abs(coef(p)))) # accuracy check
plot(He, lty = "solid")   # plots, with a bit of annotation
He5 <- He[["He5"]]
stat <- solve(deriv(He5))
points(stat, He5(stat), pch = 19)
lines(tangent(He5, stat), limits = cbind(stat-0.5, stat+0.5),
      lty = "dashed", col="black")
plot(deriv(He), lty = "solid")
plot(integral(He), lty = "solid")

## ------------------------------------------------------------------------
Hea <- poly_orth_general(Discrete, degree = 5, 
                         x = seq(-5, 5, length.out = 101),
                         w = function(x) exp(-x^2/2), norm = FALSE)
plot(He, lty = "dashed", col = "black") # accurate polynomials
lines(Hea, lty = "solid")               # discrete approximation
He5 <- Hea[[6]]
stat <- solve(deriv(He5))
points(stat, He5(stat), pch = 19)
lines(tangent(He5, stat), limits = cbind(stat-0.5, stat+0.5),
      lty = "dashed", col="black")

## ------------------------------------------------------------------------
P <- polynomial(c(1, 5, 2, 2)/10)
s <- polynomial(c(0, 1))
(mean_offspring <- deriv(P)(1))
plot(s, xlim = c(0,1), ylim = c(0,1), xlab = "s", ylab = "P(s)", col = "grey")
abline(h=0:1, v=0:1, lty = "solid", lwd = 0.1)
lines(P, limits = 0:1)
lines(tangent(P, 1), lty = "dashed", col = "red", limits = c(0.5, 1.5))  
# higher generations
plot(s, xlim = c(0,1), ylim = c(0,1), xlab = "s", ylab = "P(s)", col = "grey")
abline(h=0:1, v=0:1, lty = "solid", lwd = 0.1)
lines(P, limits = 0:1)
lines(P(P), col = 2, limits = 0:1)
lines(P(P(P)), col = 3, limits = 0:1)
lines(P(P(P(P))), col = 4, limits = 0:1)
solve(P, s)                # for the extinction probability
(ep <- solve((P-s)/(s-1))) # factor our the known zero at s = 1
ex <- Re(ep[length(ep)])   # extract the appropriate value (may be complex)
segments(ex, 0, ex, P(ex), col = "red", lwd = 2)
abline(h=0:1, v=0:1, lty = "dotted")

## ------------------------------------------------------------------------
x0 <- 80:89
y0 <- c(487, 370, 361, 313, 246, 234, 173, 128, 88, 83)

p <- poly_calc(x0, y0)        # leads to catastropic numerical failure!
p(x0) - y0                    # these should be "close to zero"!

p1 <- poly_calc(x0 - 84, y0)  # changing origin fixes the problem
p1(x0 - 84) - y0              # these are 'close to zero'.

plot(p1, xlim = c(80, 89) - 84, xlab = "x0 - 84")
points(x0 - 84, y0, col = "red")

# Can we now write the polynomial in "raw" form?
x <- polynomial()
p0 <- p1(x - 84) # attempting to change the origin back to zero
                 # leads to severe numerical problems again
plot(p0, xlim = c(80, 89))
points(x0, y0, col = "red") # major errors due to finite precision

