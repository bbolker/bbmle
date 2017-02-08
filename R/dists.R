snorm <- function(mean,sd) {
  list(title="Normal",
       mean=mean,sd=sd,
       median=mean,
       mode=mean,
       variance=sd^2,
       sd=sd)
}



#' Abstract definitions of distributions
#' 
#' Functions returning values for summary statistics (mean, median, etc.) of
#' distributions
#' 
#' 
#' @aliases sbinom spois snbinom snorm sbeta sbetabinom
#' @param prob probability as defined for \code{\link{dbinom}},
#' \code{\link{dnbinom}}, or beta-binomial distribution (\code{dbetabinom} in
#' the \code{emdbook} package)
#' @param size size parameter as defined for \code{\link{dbinom}} or
#' \code{dbetabinom} in the \code{emdbook} package, or size/overdispersion
#' parameter as in \code{\link{dnbinom}}
#' @param mean mean parameter as defined for \code{\link{dnorm}}
#' @param mu mean parameter as defined for \code{\link{dnbinom}}
#' @param sd standard deviation parameter as defined for \code{\link{dnorm}}
#' @param shape1 shape parameter for \code{\link{dbeta}}
#' @param shape2 shape parameter for \code{\link{dbeta}}
#' @param lambda rate parameter as defined for \code{\link{dpois}}
#' @param theta overdispersion parameter for beta-binomial (see
#' \code{dbetabinom} in the \code{emdbook} package)
#' @return \item{title}{name of the distribution} \item{[parameters]}{input
#' parameters for the distribution} \item{mean}{theoretical mean of the
#' distribution} \item{median}{theoretical median of the distribution}
#' \item{mode}{theoretical mode of the distribution}
#' \item{variance}{theoretical variance of the distribution}
#' \item{sd}{theoretical standard deviation of the distribution}
#' @note these definitions are tentative, subject to change as I figure this
#' out better.  Perhaps construct functions that return functions? Strip down
#' results? Do more automatically?
#' @author Ben Bolker
#' @seealso \code{\link{dbinom}}, \code{\link{dpois}}, \code{\link{dnorm}},
#' \code{\link{dnbinom}}
#' @keywords misc
#' @examples
#' 
#'   sbinom(prob=0.2,size=10)
#'   snbinom(mu=2,size=1.2)
#' 
sbinom <- function(size,prob) {
  list(title="Binomial",
       prob=prob,size=size,
       mean=prob*size,
       median=qbinom(0.5,size,prob),
       mode=NA,
       variance=size*prob*(1-prob),
       sd=sqrt(size*prob*(1-prob)),
       formula="x*log(prob)+(size-x)*log(1-prob)")
}

sbeta <- function(shape1,shape2) {
  list(title="Beta",
       shape1=shape1,shape2=shape2,
       mean=shape1/(shape1+shape2),
       median=qbeta(0.5,shape1,shape2),
       mode=NA,
       variance=shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)),
       sd=sqrt(shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))))
}

snbinom <- function(size,prob,mu) {
    if (missing(mu) && !missing(prob)) {
        mupar <- FALSE
        mu = NA ## FIXME
        warning("STUB in snbinom: calc. mu as a function of prob")
    }
    if (!missing(mu) && missing(prob)) {
        mupar <- TRUE
        prob = size/(size+mu)
    }
    v <- if (mupar) mu+mu^2/size else size*(1-prob)/prob^2
    list(title="Negative binomial",
         prob=prob,mu=mu,size=size,
         mean=if (mupar) mu else size*(1-prob)/prob,
         median= if (mupar) qnbinom(0.5,mu=mu,size) else qnbinom(0.5,prob=prob,size),
         mode=NA,
         variance=v,
         sd=sqrt(v))
}

spois <- function(lambda) {
  list(title="Poisson",
       lambda=lambda,
       mean=lambda,
       median=qpois(0.5,lambda),
       mode=NA,
       variance=lambda,
       sd=sqrt(lambda))      
}

sbetabinom <- function(size,prob,theta) {
  list(title="Beta-binomial",
       prob=prob,size=size,theta=theta,
       mean=prob*size,
       median=NA, ## qbetabinom(0.5,size,prob),
       mode=NA,
       variance=size*prob*(1-prob)/theta,
       sd=sqrt(size*prob*(1-prob)))
}

sgamma <- function(shape,rate=1,scale=1/rate) {
    if (missing(rate)) rate <- 1/scale
    list(title="Gamma",
         mean=shape/rate,sd=sqrt(shape)/rate,
         median=NA,
         mode=NA,
         variance=shape/rate^2)
}
