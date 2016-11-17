\name{mle2.options}
\alias{mle2.options}
\title{Options for maximum likelihood estimation}
\description{
  Query or set MLE parameters
}
\usage{
mle2.options(...)
}
\arguments{
  \item{\dots}{names of arguments to query,
    or a list of values to set}
}
\details{
  \itemize{
    \item{optim.method}{name of optimization method (see
      \code{\link{optim}} for choices)}
    \item{confint}{name of confidence-interval: choices
      are "spline", "uniroot", "hessian" corresponding
      to spline inversion, attempt to find best answer via uniroot,
      information-matrix approximation}
    \item{optimizer}{optimization function to use by default
(choices: "optim", "nlm", "nlminb", "constrOptim")}
}
}
\value{
  Values of queried parameters, or (invisibly) the full list of parameters
}
\seealso{
  \code{\link{mle2-class}}
}
\keyword{models}
