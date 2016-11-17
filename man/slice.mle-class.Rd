\name{slice.mle2-class}
\docType{class}
\alias{slice.mle2-class}
\title{likelihood-surface slices}
\description{evaluations of log-likelihood along transects in parameter space}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("slice.mle2", ...)}.
The objects are similar to likelihood profiles, but don't involve
any optimization with respect to the other parameters.
}
\section{Slots}{
  \describe{
    \item{\code{profile}:}{Object of class \code{"list"}. List of
      slices, one for each requested parameter. Each slice is a data
      frame with the first column called \code{z} being the signed square
      root of the -2 log likelihood ratio, and the others being the
      parameters with names prefixed by \code{par.vals.}}
    \item{\code{summary}:}{Object of class \code{"summary.mle2"}. Summary
      of object being profiled.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "profile.mle2", y = "missing")}: Plot
      profiles for each parameter.}
  }
}
\seealso{
  \code{\link{profile.mle2-class}}
}
\keyword{classes}
