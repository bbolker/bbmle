\name{profile-methods}
\docType{methods}
\alias{profile-methods}
\alias{profile,mle2-method}
\alias{profile.mle2}
\title{Likelihood profiles }
\description{
  Compute likelihood profiles for a fitted model
}

\usage{
\S4method{profile}{mle2}(fitted, which = 1:p, maxsteps = 100,
                    alpha = 0.01, zmax = sqrt(qchisq(1 - alpha/2, p)),
                    del = zmax/5, trace = FALSE, skiperrs=TRUE,
                    std.err, 
                    tol.newmin = 0.001, debug=FALSE,
                    prof.lower, prof.upper,
skip.hessian = TRUE, try_harder=FALSE, \dots)
}
\arguments{
  \item{fitted}{A fitted maximum likelihood model of class
\dQuote{mle2}}
  \item{which}{a numeric or character vector describing which parameters
to profile (default is to profile all parameters)}
  \item{maxsteps}{maximum number of steps to take looking for an upper
value of the negative log-likelihood}
  \item{alpha}{maximum (two-sided) likelihood ratio test confidence
    level to find}
  \item{zmax}{maximum value of signed square root of deviance difference
    to find (default value corresponds to a
  2-tailed chi-squared test at level alpha)}
  \item{del}{step size for profiling}
  \item{trace}{(logical) produce tracing output?}
  \item{skiperrs}{(logical) ignore errors produced during profiling?}
  \item{std.err}{Optional numeric vector of standard errors, for cases
    when the Hessian is badly behaved.  Will be replicated if necessary,
    and NA values will be replaced by the corresponding values from
    the fit summary}
  \item{tol.newmin}{tolerance for diagnosing a new minimum below the
    minimum deviance estimated in initial fit is found}
  \item{debug}{(logical) debugging output?}
  \item{prof.lower}{optional vector of lower bounds for profiles}
  \item{prof.upper}{optional vector of upper bounds for profiles}
  \item{skip.hessian}{skip hessian (defunct?)}
  \item{try_harder}{(logical) ignore \code{NA} and flat spots in the
profile, try to continue anyway?}
  \item{\dots}{additional arguments (not used)}
}
\details{
  See the vignette (\code{vignette("mle2",package="bbmle")})
  for more technical details of how profiling is done.
 }
 \seealso{\code{\link{profile.mle-class}}}
\keyword{methods}
