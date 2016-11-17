\name{summary.mle2-class}
\docType{class}
\alias{summary.mle2-class}
\alias{coef,summary.mle2-method}
\alias{show,summary.mle2-method}
\title{Class "summary.mle2", summary of "mle2" objects}
\description{Extract of "mle2" object}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.mle2",
  ...)}, but most often by invoking \code{summary} on an "mle2" object.
They contain values meant for printing by \code{show}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"language"} The call that
      generated the "mle2" object.}
    \item{\code{coef}:}{Object of class \code{"matrix"}. Estimated
      coefficients and standard errors }
    \item{\code{m2logL}:}{Object of class \code{"numeric"}. Minus twice
      the log likelihood.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.mle2")}: Pretty-prints
      \code{object}  }
    \item{coef}{\code{signature(object = "summary.mle2")}: Extracts the
    contents of the \code{coef} slot}
  }
}
\seealso{
  \code{\link{summary}}, \code{\link{mle2}}, \code{\link{mle2-class}} 
}
\keyword{classes}
