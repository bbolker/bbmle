\name{dnorm_n}
\alias{dnorm_n}
\title{
  Normal distribution with profiled-out standard deviation
}
\description{
  Returns the Normal probability densities for a distribution
  with the given mean values and the standard deviation equal
  to the root mean-squared deviation between x and mu
}
\usage{
dnorm_n(x, mean, log = FALSE)
}
\arguments{
  \item{x}{numeric vector of data}
  \item{mean}{numeric vector or mean values}
  \item{log}{logical: return the log-density?}
}
\details{
  This is a convenience function, designed for the case
  where you're trying to compute a MLE for the mean but
  don't want to bother estimating the MLE for the standard
  deviation at the same time
}
\value{
  Numeric vector of probability densities
}
\examples{
set.seed(101)
x <- rnorm(5,mean=3,sd=2)
dnorm_n(x,mean=3,log=TRUE)
}
\keyword{distribution}
