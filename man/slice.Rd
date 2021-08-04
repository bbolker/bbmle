\name{slice}
\alias{slice}
\alias{sliceOld}
\alias{slicetrans}
\alias{slice1D}
\alias{slice2D}
\title{Calculate likelihood "slices"}
\description{
  Computes cross-section(s) of a multi-dimensional likelihood surface
}
\usage{
slice(x, dim=1, ...)
sliceOld(fitted, which = 1:p, maxsteps = 100,
                       alpha = 0.01, zmax = sqrt(qchisq(1 - alpha/2, p)),
                       del = zmax/5, trace = FALSE,
                       tol.newmin=0.001, \dots)
slice1D(params,fun,nt=101,lower=-Inf,
                    upper=Inf,verbose=TRUE, tranges=NULL,
                    fun_args = NULL,
                    \dots)
slice2D(params,fun,nt=31,lower=-Inf,
                    upper=Inf,
                    cutoff=10,verbose=TRUE,
                    tranges=NULL,
                    \dots)
slicetrans(params, params2, fun, extend=0.1, nt=401,
                       lower=-Inf, upper=Inf)
}
\arguments{
  \item{x}{a fitted model object of some sort}
  \item{dim}{dimensionality of slices (1 or 2)}
  \item{params}{a named vector of baseline parameter values}
  \item{params2}{a vector of parameter values}
  \item{fun}{an objective function}
  \item{fun_args}{additional arguments to pass to \code{fun}}
  \item{nt}{(integer) number of slice-steps to take}
  \item{lower}{lower bound(s) (stub?)}
  \item{upper}{upper bound(s) (stub?)}
  \item{cutoff}{maximum increase in objective function to allow
    when computing ranges}
  \item{extend}{(numeric) fraction by which to extend range beyond specified points}
  \item{verbose}{print verbose output?}
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
  \item{tol.newmin}{tolerance for diagnosing a new minimum below the
    minimum deviance estimated in initial fit is found}
  \item{tranges}{a two-column matrix giving lower and upper bounds for
    each parameter}
  \item{\dots}{additional arguments (not used)}
}
\value{
  An object of class \code{slice} with
  \describe{
    \item{slices}{a list of individual parameter (or parameter-pair)
    slices, each of which is a data frame with elements
    \describe{
      \item{var1}{name of the first variable}
      \item{var2}{(for 2D slices) name of the second variable}
      \item{x}{parameter values}
      \item{y}{(for 2D slices) parameter values}
      \item{z}{slice values}
      \item{ranges}{a list (?) of the ranges for each parameter}
      \item{params}{vector of baseline parameter values}
      \item{dim}{1 or 2}
    }
  }
  \code{sliceOld} returns instead a list with elements \code{profile}
  and \code{summary} (see \code{\link{profile.mle2}})
}
}
\details{
  Slices provide a lighter-weight way to explore likelihood surfaces
  than profiles, since they vary a single parameter rather than
  optimizing over all but one or two parameters.

\describe{
\item{slice}{is a generic method}
\item{slice1D}{creates one-dimensional slices, by default of all
parameters of a model}
\item{slice2D}{creates two-dimensional slices, by default of all pairs
  of parameters in a model. In each panel the closed point represents
  the parameters given (typically the MLEs), while the open point
  represents the observed minimum value within the 2D slice.
  If everything has gone according to plan, these points should
  coincide (at least up to grid precision).
}
\item{slicetrans}{creates a slice along a transect between two specified
points in parameter space (see \code{calcslice} in the \code{emdbook}
package)}
}
}
\author{Ben Bolker}
\seealso{\code{\link{profile}}}
\examples{
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
fit1 <- mle2(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
   start=list(lymax=0,lhalf=0),
   data=d)
s1 <- bbmle::slice(fit1,verbose=FALSE)
s2 <- bbmle::slice(fit1,dim=2,verbose=FALSE)
require(lattice)
plot(s1)
plot(s2)
## 'transect' slice, from best-fit values to another point
st <- bbmle::slice(fit1,params2=c(5,0.5))
plot(st)
}
\keyword{misc}
