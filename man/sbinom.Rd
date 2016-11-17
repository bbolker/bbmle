\name{sbinom}
\alias{sbinom}
\alias{spois}
\alias{snbinom}
\alias{snorm}
\alias{sbeta}
\alias{sbetabinom}
\title{Abstract definitions of distributions}
\description{
  Functions returning values for summary statistics
  (mean, median, etc.) of distributions
}
\usage{
sbeta(shape1, shape2)
sbetabinom(size, prob, theta)
sbinom(size, prob)
snbinom(size, prob, mu)
snorm(mean, sd)
spois(lambda)
}
\arguments{
  \item{prob}{probability as defined for \code{\link{dbinom}},
    \code{\link{dnbinom}}, or beta-binomial distribution
    (\code{dbetabinom} in the \code{emdbook} package)}
  \item{size}{size parameter as defined for
    \code{\link{dbinom}} or \code{dbetabinom}
    in the \code{emdbook} package, or size/overdispersion parameter
    as in \code{\link{dnbinom}}}
  \item{mean}{mean parameter as defined for \code{\link{dnorm}}}
  \item{mu}{mean parameter as defined for \code{\link{dnbinom}}}
  \item{sd}{standard deviation parameter as defined for \code{\link{dnorm}}}
  \item{shape1}{shape parameter for \code{\link{dbeta}}}
  \item{shape2}{shape parameter for \code{\link{dbeta}}}
  \item{lambda}{rate parameter as defined for \code{\link{dpois}}}
  \item{theta}{overdispersion parameter for beta-binomial
    (see \code{dbetabinom} in the \code{emdbook} package)}
}
\value{
  \item{title}{name of the distribution}
  \item{[parameters]}{input parameters for the distribution}
  \item{mean}{theoretical mean of the distribution}
  \item{median}{theoretical median of the distribution}
  \item{mode}{theoretical mode of the distribution}
  \item{variance}{theoretical variance of the distribution}
  \item{sd}{theoretical standard deviation of the distribution}
}
\author{Ben Bolker}
\seealso{\code{\link{dbinom}}, \code{\link{dpois}}, \code{\link{dnorm}}, \code{\link{dnbinom}}}
\examples{
  sbinom(prob=0.2,size=10)
  snbinom(mu=2,size=1.2)
}
\note{these definitions are tentative, subject to change
as I figure this out better.  Perhaps construct functions
that return functions? Strip down results? Do more
automatically?}

\keyword{misc}
