## from Eric Weese
library(bbmle)
f <- function(x=2,a=1) x^2 - a
f.g <- function(x,a) 2*x
f.g2 <- function(x,a) c(2*x,0)
options(digits=3)
m1 <- mle2(f,fixed=list(a=1))
m2 <- mle2(f,gr=f.g,fixed=list(a=1))
m3 <- mle2(f,gr=f.g2,fixed=list(a=1))
stopifnot(all.equal(coef(m1),coef(m2)))
stopifnot(all.equal(coef(m1),coef(m3)))
tt <- function(x) x@details$hessian
stopifnot(all.equal(tt(m1),tt(m2),tolerance=1e-6))
stopifnot(all.equal(tt(m1),tt(m3),tolerance=1e-6))
