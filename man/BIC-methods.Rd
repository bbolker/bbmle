\name{BIC-methods}
\docType{methods}
%\alias{BIC}
\alias{BIC-methods}
\alias{AIC-methods}
\alias{AICc-methods}
\alias{logLik-methods}
\alias{AICc}
\alias{AIC,mle2-method}
\alias{AICc,mle2-method}
\alias{AICc,logLik-method}
\alias{AICc,ANY-method}
\alias{AICc,ANY,mle2,logLik-method}
\alias{qAICc}
\alias{qAICc-methods}
\alias{qAICc,ANY-method}
\alias{qAICc,mle2-method}
\alias{qAICc,logLik-method}
\alias{qAIC}
\alias{qAIC-methods}
\alias{qAIC,ANY-method}
\alias{qAIC,mle2-method}
\alias{qAIC,logLik-method}
%\alias{BIC,logLik-method}
%\alias{BIC,ANY-method}
%\alias{BIC,mle2-method}
%\alias{BIC,ANY,mle2,logLik-method}
\alias{qAIC,ANY,mle2,logLik-method}
\alias{qAICc,ANY,mle2,logLik-method}
\alias{logLik,mle2-method}
\alias{anova,mle2-method}
\title{Log likelihoods and model selection for mle2 objects}
\description{
  Various functions for likelihood-based and information-theoretic
  model selection of likelihood models
}
\section{Methods}{
\describe{
    \item{logLik}{\code{signature(object = "mle2")}: Extract maximized
      log-likelihood.}
    \item{AIC}{\code{signature(object = "mle2")}: Calculate
      Akaike Information Criterion}
    \item{AICc}{\code{signature(object = "mle2")}: Calculate
      small-sample corrected Akaike Information Criterion}
    %\item{BIC}{\code{signature(object = "mle2")}: Calculate
    %Bayesian (Schwarz) Information Criterion}
    %\item{BIC}{\code{signature(object = "logLik")}: Calculate
    %Bayesian (Schwarz) Information Criterion}
    %\item{BIC}{\code{signature(object = "ANY")}: Calculate
    %Bayesian (Schwarz) Information Criterion}
    \item{anova}{\code{signature(object="mle2")}: Likelihood Ratio Test
      comparision of different models}
  }
}
\usage{
%\S4method{BIC}{ANY,mle2,logLik}(object,...)
\S4method{AICc}{ANY,mle2,logLik}(object,...,nobs,k=2)
\S4method{qAIC}{ANY,mle2,logLik}(object,...,k=2)
\S4method{qAICc}{ANY,mle2,logLik}(object,...,nobs,k=2)
}
\arguments{
  \item{object}{A \code{logLik} or \code{mle2} object}
  \item{...}{An optional list of additional \code{logLik}
    or \code{mle2} objects (fitted to the same data set).}
  \item{nobs}{Number of observations (sometimes 
    obtainable as an attribute of
    the fit or of the log-likelihood)}
  \item{k}{penalty parameter (nearly always left at its default value of 2)}
}
\details{
  Further arguments to \code{BIC} can be specified
  in the \code{...} list: \code{delta} (logical)
  specifies whether to include a column for delta-BIC
  in the output.
}
\value{
  A table of the BIC values, degrees of freedom,
  and possibly delta-BIC values relative to the
  minimum-BIC model
}
\note{This is implemented in an ugly way and could
  probably be improved!}
\examples{
  d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
  (fit <- mle2(y~dpois(lambda=ymax/(1+x/xhalf)),
      start=list(ymax=25,xhalf=3),data=d))
  (fit2 <- mle2(y~dpois(lambda=(x+1)*slope),
      start=list(slope=1),data=d))
  BIC(fit)
  BIC(fit,fit2)
  }
\keyword{methods}
