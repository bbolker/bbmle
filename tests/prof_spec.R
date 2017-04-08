## test whether profiling works when custom optimizer is defined
##  inside a function (GH #7)

library(bbmle)
test <- function(t, X) {
  likfun <- function(p) {
    mu <- with(as.list(p), {
      exp(a+b*t)
    })
    -sum(dpois(X, mu, log=TRUE))
  }
  parnames(likfun) <- c("a", "b")
  
  optimfun <- function(par, fn, gr = NULL, ...,
                       method = NULL, lower = -Inf, upper = Inf,
                       control = NULL, hessian = FALSE) {
    ## cat("using custom optimfun!\n")
    optim(par, fn=fn, gr=gr, ...,
          method="BFGS", control=control, hessian=TRUE)
  }
  
  mle2(likfun, start=c(a=1,b=1), optimizer="user", optimfun=optimfun)
}

f <- test(0:5, round(exp(1:6)))
pp <- profile(f,skiperrs=FALSE)
stopifnot(inherits(pp,"profile.mle2"))


