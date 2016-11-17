\name{as.data.frame.profile.mle2}
\alias{as.data.frame.profile.mle2}
\alias{coerce,profile.mle2-method}
\alias{coerce,profile.mle2,data.frame-method}
\title{convert profile to data frame}
\description{
  converts a profile of a fitted mle2 object
  to a data frame
}
\usage{
\S3method{as.data.frame}{profile.mle2}(x, row.names=NULL,
optional=FALSE, \dots)
}
\arguments{
  \item{x}{a profile object}
  \item{row.names}{row names (unused)}
  \item{optional}{unused}
  \item{\dots}{unused}
}
\value{
  a data frame with columns
  \item{param}{name of parameter being profiled}
  \item{z}{signed square root of the deviance difference from the
    minimum}
  \item{parameter values}{named par.vals.parname}
  \item{focal}{value of focal parameter: redundant, but included for
    plotting convenience}
}
\examples{
  ## use as.data.frame and lattice to plot profiles
  x <- 0:10
  y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
  library(bbmle)
  LL <- function(ymax=15, xhalf=6)
  -sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE))
  ## uses default parameters of LL
  fit1 <- mle2(LL)
  p1 <- profile(fit1)
  d1 = as.data.frame(p1)
  library(lattice)
  xyplot(abs(z)~focal|param,data=d1,
    subset=abs(z)<3,
    type="b",
    xlab="",
    ylab=expression(paste(abs(z),
    " (square root of ",Delta," deviance)")),
    scale=list(x=list(relation="free")))
}
\author{Ben Bolker}
\keyword{misc}
