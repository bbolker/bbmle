\name{parnames}
\alias{parnames}
\alias{parnames<-}
\title{get and set parameter names}
\description{
  Gets and sets the "parnames" attribute on a
  negative log-likelihood function
}
\usage{
parnames(obj)
parnames(obj) <- value
}
\arguments{
  \item{obj}{a negative log-likelihood function}
  \item{value}{a character vector of parameter names}
}
\details{
  The \code{parnames} attribute is used by \code{mle2()}
  when the negative log-likelihood function takes a parameter
  vector, rather than a list of parameters; this allows
  users to use the same objective function for \code{optim()}
  and \code{mle2()}
}
\value{
  Returns the \code{parnames} attribute (a character vector of
  parameter names) or sets it.
}
\author{Ben Bolker}
\examples{
x <- 1:5
set.seed(1001)
y <- rbinom(5,prob=x/(1+x),size=10)
mfun <- function(p) {
  a <- p[1]
  b <- p[2]
  -sum(dbinom(y,prob=a*x/(b+x),size=10,log=TRUE))
}
optim(fn=mfun,par=c(1,1))
parnames(mfun) <- c("a","b")
mle2(minuslogl=mfun,start=c(a=1,b=1),method="Nelder-Mead")
}
\keyword{misc}

