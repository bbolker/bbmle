library(bbmle)

## Simulate data

set.seed(1)
x <- 1:5
y <- 2*x+1
noise <- rnorm(5, 0, 0.1)
mydata <- data.frame(x = x, y=y+noise)

## Model definition

model <- function(a, b) with(mydata, a*x+b)

## Negative log-likelihood

nll <- function(par) with(mydata, {
  a <- par[1]
  b <- par[2]
  sum(0.5*((y-model(a,b))/0.1)^2)
  
})

gr <- function(par) with(mydata, {
  a <- par[1]
  b <- par[2]
  dnllda <- -sum(((y-model(a,b))/0.1)*x/0.1)
  dnlldb <- -sum(((y-model(a,b))/0.1)*1/0.1)
  return(c(dnllda, dnlldb))
})

## optimization

parnames(nll) <- c("a", "b")
parnames(gr) <- c("a", "b")

fit <- mle2(nll, c(a = 1, b=2), gr=gr)

myprof <- profile(fit)
myprof_c <- profile(fit,continuation="naive")
confint(myprof)
confint(myprof_c)

fit <- mle2(nll, c(a = 1, b=2), gr=gr, skip.hessian=TRUE)
myprof2 <- profile(fit,std.err=c(0.1,0.1))

## incomplete!
model2 <- ~a+b*x+c*x^2
f0 <- deriv(model2,"x",function.arg=c("a","b","c"))
## chain rule
f1 <- function() {
## memoize
lastpar <- NULL
lastval <- NULL
}

f2 <- function(par) {
    if (par==lastpar) {
        return(c(lastval))
    }
    lastpar <<- par
    lastval <<- do.call(f0,par)
    f1(par)
}
f2.gr <- function(par) {
    if (par==lastpar) {
        return(attr(lastval,".grad"))
    }
    lastpar <<- par
    lastval <<- do.call(f0,par)
    f1.gr(par)
}
parnames(f2) <- parnames(f2.gr) <- c("a","b","c")
