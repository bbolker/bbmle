\name{mle2-class}
\docType{class}
\alias{mle2-class}
\alias{coef,mle2-method}
\alias{show,mle2-method}
\alias{slice,mle2-method}
\alias{summary,mle2-method}
\alias{update,mle2-method}
\alias{vcov,mle2-method}
\alias{deviance,mle2-method}
\alias{coerce,mle,mle2-method}
\alias{formula,mle2-method}
\alias{stdEr}
\alias{stdEr,mle2-method}
\title{Class "mle2". Result of Maximum Likelihood Estimation.}
\description{This class encapsulates results of a generic maximum
  likelihood procedure.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("mle2", \dots)}, but
  most often as the result of a call to \code{\link{mle2}}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{(language) The call to \code{\link{mle2}}.}
    \item{\code{call.orig}:}{(language) The call to \code{\link{mle2}},
      saved in its original form (i.e. without data arguments
      evaluated).}
    \item{\code{coef}:}{(numeric) Vector of estimated parameters.}
    \item{\code{data}:}{(data frame or list) Data with which to evaluate the negative log-likelihood function}
    \item{\code{fullcoef}:}{(numeric) Fixed and estimated parameters.}
    \item{\code{vcov}:}{(numeric matrix) Approximate variance-covariance
      matrix, based on the second derivative matrix at the MLE.}
    \item{\code{min}:}{(numeric) Minimum value of objective function =
      minimum negative log-likelihood.}
    \item{\code{details}:}{(list) Return value from \code{\link{optim}}.}
    \item{\code{minuslogl}:}{(function) The negative log-likelihood
      function.}
    \item{\code{optimizer}:}{(character) The optimizing function used.}
    \item{\code{method}:}{(character) The optimization method used.}
    \item{\code{formula}:}{(character) If a formula was specified, a
      character vector giving the formula and parameter specifications.}
}
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "mle2")}: Extract coefficients.
      If \code{exclude.fixed=TRUE} (it is \code{FALSE} by default),
      only the non-fixed parameter values are returned.}t
    \item{confint}{\code{signature(object = "mle2")}: Confidence
      intervals from likelihood profiles, or quadratic approximations,
      or root-finding.}
    \item{show}{\code{signature(object = "mle2")}: Display object
      briefly.}
    \item{show}{\code{signature(object = "summary.mle2")}: Display object briefly.}
    \item{summary}{\code{signature(object = "mle2")}: Generate object summary.}
    \item{update}{\code{signature(object = "mle2")}:  Update fit.}
    \item{vcov}{\code{signature(object = "mle2")}: Extract
      variance-covariance matrix.}
    \item{formula}{\code{signature(object="mle2")}: Extract formula}
    \item{plot}{\code{signature(object="profile.mle2,missing")}: Plot
      profile. }
  }
}
\details{
  When the parameters in the original fit are constrained using
  \code{lower} or \code{upper}, or when \code{prof.lower} or
  \code{prof.upper} are set, and the confidence intervals lie
  outside the constraint region, \code{confint} will return \code{NA}.
  This may be too conservative -- in some cases, the appropriate
  answer would be to set the confidence limit to the lower/upper
  bound as appropriate -- but it is the most general answer.

  (If you have a strong opinion about the need for a new
  option to \code{confint} that sets the bounds to the limits
  automatically, please contact the package maintainer.)
  
}
\examples{
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
lowerbound <- c(a=2,b=-0.2)
d <- data.frame(x,y)
fit1 <- mle2(y~dpois(lambda=exp(a+b*x)),start=list(a=0,b=2),data=d,
method="L-BFGS-B",lower=c(a=2,b=-0.2))
(cc <- confint(fit1,quietly=TRUE))
## to set the lower bounds to the limit
na_lower <- is.na(cc[,1])
cc[na_lower,1] <- lowerbound[na_lower]
cc
}
\keyword{classes}
