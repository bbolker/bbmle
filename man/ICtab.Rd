\name{ICtab}
\alias{ICtab}
\alias{AICtab}
\alias{BICtab}
\alias{AICctab}
\alias{print.ICtab}
\title{Compute table of information criteria and auxiliary info}
\description{
  Computes information criteria for a series of models, optionally
  giving information about weights, differences between ICs, etc.
}
\usage{
ICtab(\dots, type=c("AIC","BIC","AICc","qAIC","qAICc"),
    weights = FALSE, delta = TRUE, base = FALSE,
logLik=FALSE, sort = TRUE,
nobs=NULL, dispersion = 1, mnames, k = 2)
AICtab(\dots,mnames)
BICtab(\dots,mnames)
AICctab(\dots,mnames)
\method{print}{ICtab}(x,\dots,min.weight)
}
\arguments{
  \item{\dots}{a list of (logLik or?) mle objects; in the case of
    \code{AICtab} etc., could also include other arguments to \code{ICtab}}
  \item{type}{specify information criterion to use}
  \item{base}{(logical) include base IC (and log-likelihood) values?}
  \item{weights}{(logical) compute IC weights?}
  \item{logLik}{(logical) include log-likelihoods in the table?}
  \item{delta}{(logical) compute differences among ICs (and log-likelihoods)?}
  \item{sort}{(logical) sort ICs in increasing order?}
  \item{nobs}{(integer) number of observations: required for
    \code{type="BIC"} or \code{type="AICc"} unless objects have
    a \code{\link{nobs}} method}
  \item{dispersion}{overdispersion estimate, for computing qAIC:
    required for \code{type="qAIC"} or \code{type="qAICc"} unless
    objects have a \code{"dispersion"} attribute}
  \item{mnames}{names for table rows: defaults to names of objects passed}
  \item{k}{penalty term (largely unused: left at default of 2)}
  \item{x}{an ICtab object}
  \item{min.weight}{minimum weight for exact reporting (smaller values
    will be reported as "<[min.weight]")}
}
\value{
  A data frame containing:
  \item{IC}{information criterion}
  \item{df}{degrees of freedom/number of parameters}
  \item{dIC}{difference in IC from minimum-IC model}
  \item{weights}{exp(-dIC/2)/sum(exp(-dIC/2))}
}
\note{(1) The print method uses sensible defaults; all ICs are rounded
  to the nearest 0.1, and IC weights are printed using
  \code{\link{format.pval}} to print an inequality for
  values <0.001. (2) The computation of degrees of freedom/number of
  parameters (e.g., whether
  variance parameters are included in the total) varies enormously
  between packages.  As long as the df computations
  for a given set of models is consistent, differences
  don't matter, but one needs to be careful with log likelihoods
  and models taken from different packages.  If necessary
  one can change the degrees of freedom manually by
  saying \code{attr(obj,"df") <- df.new}, where \code{df.new}
  is the desired number of parameters.
  (3) Defaults have changed to \code{sort=TRUE}, \code{base=FALSE},
  \code{delta=TRUE}, to match my conviction that it rarely makes
  sense to report the overall values of information criteria}
\references{Burnham and Anderson 2002}
\author{Ben Bolker}
\examples{
  set.seed(101)
  d <- data.frame(x=1:20,y=rpois(20,lambda=2))
  m0 <- glm(y~1,data=d)
  m1 <- update(m0,.~x)
  m2 <- update(m0,.~poly(x,2))
  AICtab(m0,m1,m2,mnames=LETTERS[1:3])
  AICtab(m0,m1,m2,base=TRUE,logLik=TRUE)
  AICtab(m0,m1,m2,logLik=TRUE)
  AICctab(m0,m1,m2,weights=TRUE)
  print(AICctab(m0,m1,m2,weights=TRUE),min.weight=0.1)
}
\keyword{misc}
