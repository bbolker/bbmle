\name{predict-methods}
\docType{methods}
\alias{gfun}
\alias{predict-methods}
\alias{predict,mle2-method}
\alias{residuals,mle2-method}
\alias{simulate,mle2-method}
\title{Predicted values from an mle2 fit}
\description{
  Given an \code{mle2} fit and an optional list
  of new data, return predictions (more generally,
  summary statistics of the predicted distribution)
}
\section{Methods}{
\describe{
  \item{x = "mle2"}{an \code{mle2} fit}
}}
\usage{
    \S4method{predict}{mle2}(object, newdata=NULL,
                         location="mean", newparams=NULL, \dots)
    \S4method{simulate}{mle2}(object, nsim,
                         seed, newdata=NULL, newparams=NULL, \dots)
    \S4method{residuals}{mle2}(object,type=c("pearson","response"),
                   location="mean",\dots)
}
\arguments{
  \item{object}{an mle2 object}
  \item{newdata}{optional list of new data}
  \item{newparams}{optional vector of new parameters}
  \item{location}{name of the summary statistic to return}
  \item{nsim}{number of simulations}
  \item{seed}{random number seed}
  \item{type}{residuals type}
  \item{\dots}{additional arguments (for generic compatibility)}
}
\note{For some models (e.g. constant models), \code{predict} may
  return a single value rather than a vector of the appropriate length.}
\examples{
set.seed(1002)
lymax <- c(0,2)
lhalf <- 0
x <- runif(200)
g <- factor(rep(c("a","b"),each=100))
y <- rnbinom(200,mu=exp(lymax[g])/(1+x/exp(lhalf)),size=2)
dat <- data.frame(y,g,x)

fit3 <- mle2(y~dnbinom(mu=exp(lymax)/(1+x/exp(lhalf)),size=exp(logk)),
    parameters=list(lymax~g),
    start=list(lymax=0,lhalf=0,logk=0),
data=dat)

plot(y~x,col=g)
## true curves
curve(exp(0)/(1+x/exp(0)),add=TRUE)
curve(exp(2)/(1+x/exp(0)),col=2,add=TRUE)
## model predictions
xvec = seq(0,1,length=100)
lines(xvec,predict(fit3,newdata=list(g=factor(rep("a",100),levels=c("a","b")),
                                x = xvec)),col=1,lty=2)
lines(xvec,predict(fit3,newdata=list(g=factor(rep("b",100),levels=c("a","b")),
                                x = xvec)),col=2,lty=2)


## comparing automatic and manual predictions
p1 = predict(fit3)
p2A =
with(as.list(coef(fit3)),exp(`lymax.(Intercept)`)/(1+x[1:100]/exp(lhalf)))
p2B =
with(as.list(coef(fit3)),exp(`lymax.(Intercept)`+lymax.gb)/(1+x[101:200]/exp(lhalf)))
all(p1==c(p2A,p2B))
##
simulate(fit3)
}
\keyword{methods}
