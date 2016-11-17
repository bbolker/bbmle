\name{profile.mle2-class}
\docType{class}
\alias{profile.mle2-class}
\alias{confint,profile.mle2-method}
\alias{confint,mle2-method}
\alias{confint.mle2}  %% bogus but good ref link
\alias{plot,profile.mle2-method}
\alias{plot,profile.mle2,missing-method}
\alias{show,profile.mle2-method}
\alias{plot.profile.mle2}
\title{Methods for likelihood profiles}
\description{Definition of the mle2 likelihood profile class,
  and applicable methods}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("profile.mle2",
    ...)}, but most often by invoking \code{profile} on an "mle2" object.
}
\section{Slots}{
  \describe{
    \item{\code{profile}:}{Object of class \code{"list"}. List of
      profiles, one for each requested parameter. Each profile is a data
      frame with the first column called \code{z} being the signed square
      root of the deviance, and the others being the
      parameters with names prefixed by \code{par.vals.}}
    \item{\code{summary}:}{Object of class \code{"summary.mle2"}. Summary
      of object being profiled.}
  }
}
\section{Methods}{
  \describe{
    \item{confint}{\code{signature(object = "profile.mle2")}: Use profile
      to generate approximate confidence intervals for parameters.}
    \item{plot}{\code{signature(x = "profile.mle2", y = "missing")}: Plot
      profiles for each parameter.}
    \item{summary}{\code{signature(x = "profile.mle2")}: Plot
      profiles for each parameter.}
    \item{show}{\code{signature(object = "profile.mle2")}: Show object.}
  }
}
\usage{
\S4method{plot}{profile.mle2}(x,
   levels, which=1:p, conf = c(99, 95, 90, 80, 50)/100,
   plot.confstr = TRUE,
   confstr = NULL, absVal = TRUE, add = FALSE,
   col.minval="green", lty.minval=2,
   col.conf="magenta", lty.conf=2,
   col.prof="blue", lty.prof=1,
   xlabs=nm, ylab="z",
   onepage=TRUE,
   ask=((prod(par("mfcol")) < length(which)) && dev.interactive() &&
               !onepage),
   show.points=FALSE,
   main, xlim, ylim, \dots)
\S4method{confint}{mle2}(object, parm, level = 0.95, method,
          trace=FALSE,quietly=!interactive(),
          tol.newmin=0.001,\dots)
\S4method{confint}{profile.mle2}(object, parm, level = 0.95, trace=FALSE, \dots)
}
\arguments{
  \item{x}{An object of class \code{profile.mle2}}
  \item{object}{An object of class \code{mle2} or \code{profile.mle2}
    (as appropriate)}
  \item{levels}{levels at which to plot likelihood cutoffs (set by conf
    by default)}
  \item{level}{level at which to compute confidence interval}
  \item{which}{(numeric or character) which parameter profiles to plot}
  \item{parm}{(numeric or character) which parameter(s) to find
    confidence intervals for}
  \item{method}{(character) "spline", "uniroot", or "quad", for
    spline-extrapolation-based (default), root-finding, or quadratic
    confidence intervals.  By default it uses the value of
    \code{mle2.options("confint")} -- the factory setting is "spline".}
  \item{trace}{trace progress of confidence interval calculation when using
    \sQuote{uniroot} method?}
  \item{conf}{(1-alpha) levels at which to plot likelihood
    cutoffs/confidence intervals}
  \item{quietly}{(logical) suppress \dQuote{Profiling ...} message
    when computing profile to get confidence interval?}
  \item{tol.newmin}{see \code{\link{profile-methods}}}
  \item{plot.confstr}{(logical) plot labels showing confidence levels?}
  \item{confstr}{(character) labels for confidence levels (by default, constructed from conf levels)}
  \item{absVal}{(logical) plot absolute values of signed square
    root deviance difference ("V" plot rather than straight-line
    plot)?}
  \item{add}{(logical) add profile to existing graph?}
  \item{col.minval}{color for minimum line}
  \item{lty.minval}{line type for minimum line}
  \item{col.conf}{color for confidence intervals}
  \item{lty.conf}{line type for confidence intervals}
  \item{col.prof}{color for profile}
  \item{lty.prof}{line type for profile}
  \item{xlabs}{x labels}
  \item{ylab}{y label}
  \item{onepage}{(logical) plot all profiles on one page,
    adjusting par(mfcol) as necessary?}
  \item{ask}{(logical) pause for user input between plots?}
  \item{show.points}{(logical) show computed profile points
    as well as interpolated spline?}
  \item{main}{(logical) main title}
  \item{xlim}{x limits}
  \item{ylim}{y limits}
  \item{\dots}{other arguments}
}
\seealso{
  \code{\link{mle2}}, \code{\link{mle2-class}}, \code{\link{summary.mle2-class}} 
}
\details{
  The default confidence interval calculation computes a likelihood
  profile and uses the points therein, or uses the computed points in
  an existing \code{profile.mle2} object, to construct an interpolation
  spline (which by default has three times as many points as were in
  the original set of profile points).  It then uses linear
  interpolation between these interpolated points (!)
}
\examples{
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
## we have a choice here: (1) don't impose boundaries on the parameters,
##  put up with warning messages about NaN values: 
fit1 <- mle2(y~dpois(lambda=ymax/(1+x/xhalf)),
     start=list(ymax=1,xhalf=1),
     data=d)
p1 <- suppressWarnings(profile(fit1))
plot(p1,main=c("first","second"),
     xlab=c(~y[max],~x[1/2]),ylab="Signed square root deviance",
     show.points=TRUE)
suppressWarnings(confint(fit1)) ## recomputes profile
confint(p1)  ## operates on existing profile
suppressWarnings(confint(fit1,method="uniroot"))
## alternatively, we can use box constraints to keep ourselves
##  to positive parameter values ...
fit2 <- update(fit1,method="L-BFGS-B",lower=c(ymax=0.001,xhalf=0.001))
p2 <- profile(fit2)
plot(p2,show.points=TRUE)
## but the fit for ymax is just bad enough that the spline gets wonky
confint(p2)  ## now we get a warning
confint(fit2,method="uniroot")
## bobyqa is a better-behaved bounded optimizer ...
##  BUT recent (development, 2012.5.24) versions of
##    optimx no longer allow single-parameter fits!
\dontrun{
if (require(optimx)) {
  fit3 <- update(fit1,
      optimizer="optimx",
      method="bobyqa",lower=c(ymax=0.001,xhalf=0.001))
   p3 <- profile(fit3)
   plot(p3,show.points=TRUE)
  confint(p3)
}
}
}
\keyword{classes}
