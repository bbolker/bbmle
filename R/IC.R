#' Compute table of information criteria and auxiliary info
#' 
#' Computes information criteria for a series of models, optionally giving
#' information about weights, differences between ICs, etc.
#' 
#' 
#' @aliases ICtab AICtab BICtab AICctab print.ICtab
#' @param \dots a list of (logLik or?) mle objects; in the case of
#' \code{AICtab} etc., could also include other arguments to \code{ICtab}
#' @param type specify information criterion to use
#' @param base (logical) include base IC (and log-likelihood) values?
#' @param weights (logical) compute IC weights?
#' @param logLik (logical) include log-likelihoods in the table?
#' @param delta (logical) compute differences among ICs (and log-likelihoods)?
#' @param sort (logical) sort ICs in increasing order?
#' @param nobs (integer) number of observations: required for \code{type="BIC"}
#' or \code{type="AICc"} unless objects have a \code{\link{nobs}} method
#' @param dispersion overdispersion estimate, for computing qAIC: required for
#' \code{type="qAIC"} or \code{type="qAICc"} unless objects have a
#' \code{"dispersion"} attribute
#' @param mnames names for table rows: defaults to names of objects passed
#' @param k penalty term (largely unused: left at default of 2)
#' @param x an ICtab object
#' @param min.weight minimum weight for exact reporting (smaller values will be
#' reported as "<[min.weight]")
#' @return A data frame containing: \item{IC}{information criterion}
#' \item{df}{degrees of freedom/number of parameters} \item{dIC}{difference in
#' IC from minimum-IC model} \item{weights}{exp(-dIC/2)/sum(exp(-dIC/2))}
#' @note (1) The print method uses sensible defaults; all ICs are rounded to
#' the nearest 0.1, and IC weights are printed using \code{\link{format.pval}}
#' to print an inequality for values <0.001. (2) The computation of degrees of
#' freedom/number of parameters (e.g., whether variance parameters are included
#' in the total) varies enormously between packages.  As long as the df
#' computations for a given set of models is consistent, differences don't
#' matter, but one needs to be careful with log likelihoods and models taken
#' from different packages.  If necessary one can change the degrees of freedom
#' manually by saying \code{attr(obj,"df") <- df.new}, where \code{df.new} is
#' the desired number of parameters.  (3) Defaults have changed to
#' \code{sort=TRUE}, \code{base=FALSE}, \code{delta=TRUE}, to match my
#' conviction that it rarely makes sense to report the overall values of
#' information criteria
#' @author Ben Bolker
#' @references Burnham and Anderson 2002
#' @keywords misc
#' @examples
#' 
#'   set.seed(101)
#'   d <- data.frame(x=1:20,y=rpois(20,lambda=2))
#'   m0 <- glm(y~1,data=d)
#'   m1 <- update(m0,.~x)
#'   m2 <- update(m0,.~poly(x,2))
#'   AICtab(m0,m1,m2,mnames=LETTERS[1:3])
#'   AICtab(m0,m1,m2,base=TRUE,logLik=TRUE)
#'   AICtab(m0,m1,m2,logLik=TRUE)
#'   AICctab(m0,m1,m2,weights=TRUE)
#'   print(AICctab(m0,m1,m2,weights=TRUE),min.weight=0.1)
#' @export
ICtab <- function(...,type=c("AIC","BIC","AICc","qAIC","qAICc"),
                  weights=FALSE,delta=TRUE,base=FALSE,
                  logLik=FALSE,
                  sort=TRUE,nobs=NULL,dispersion=1,mnames,k=2) {
    ## TO DO: allow inclusion of log-likelihood (or negative log-likelihood?)
    ## base or delta? or both?  Should deltas include delta-df as well?
    L <- list(...)
    if (is.list(L[[1]]) && length(L)==1) L <- L[[1]]
    type <- match.arg(type)
    if (dispersion !=1) {
        if (type=="BIC") stop("cannot specify dispersion with BIC")
        if (substr(type,1,1)!="q") {
            type = paste("q",type,sep="")
            warning("dispersion!=1, type changed to ",type)
        }
    }
    if (type=="AICc" || type=="BIC" || type=="qAICc") {
        if (is.null(nobs)) {
            ## if(is.null(attr(L[[1]],"nobs")))
            ## stop("must specify number of observations if corr=TRUE")
            ## nobs <- sapply(L,attr,"nobs")
            nobs <- sapply(L,nobs)
            if (length(unique(nobs))>1)
                stop("nobs different: must have identical data for all objects")
            nobs <- nobs[1]
        }
    }
    ICs <- switch(type,
                  AIC=sapply(L,AIC),
                  BIC=sapply(L,BIC),
                  AICc=sapply(L,AICc,nobs=nobs),
                  qAIC=sapply(L,qAIC,dispersion=dispersion),
                  qAICc=sapply(L,qAICc,nobs=nobs,dispersion=dispersion))
    logLiks <- sapply(L,function(x) c(logLik(x)))
    ## hack: protect against aod method
    if (is.matrix(ICs)) ICs <- ICs["AIC",]  
    getdf <- function(x) {
        if (!is.null(df <- attr(x,"df"))) return(df)
        else if (!is.null(df <- attr(logLik(x),"df"))) return(df)
    }
    dIC <- ICs-min(ICs,na.rm=TRUE)
    dlogLiks <- logLiks-min(logLiks,na.rm=TRUE)
    df <- sapply(L,getdf)
    tab <- data.frame(df=df)
    if (delta) {
        dName <- paste0("d",type)
        tab <- cbind(setNames(data.frame(dIC),dName),tab)
        if (logLik) {
            tab <- cbind(data.frame(dLogLik=dlogLiks),tab)
        }
    }
    if (base) {
        tab <- cbind(setNames(data.frame(ICs),type),tab)
        if (logLik) {
            tab <- cbind(data.frame(logLik=logLiks),tab)
        }
    }
    if (!delta && !base) stop("either 'base' or 'delta' must be TRUE")
    if (weights) {
        dIC_noNA <- na.exclude(dIC)
        wts <- napredict(attr(dIC_noNA,"na.action"),
                                    exp(-dIC_noNA/2)/sum(exp(-dIC_noNA/2)))
        tab <- data.frame(tab,weight=wts)
    }
    if (missing(mnames)) {
        Call <- match.call()
        if (!is.null(names(Call))) {
            xargs <- which(names(Call) %in% names(formals())[-1])
        } else xargs <- numeric(0)
        mnames <- as.character(Call)[c(-1,-xargs)]
    }
    row.names(tab) <- mnames
    if (sort) {
        tab <- tab[order(ICs),]
    }
    class(tab) <- "ICtab"
    tab
}

print.ICtab <- function(x,...,min.weight=0.001) {
    chtab <- format(do.call("cbind",lapply(x,round,1)))
    rownames(chtab) <- attr(x,"row.names")
    chtab[,"df"] <- as.character(x$df)
    if (!is.null(x$weight))
        chtab[,"weight"] <- format.pval(x$weight,eps=min.weight,
                                        digits=2)
    print(chtab,quote=FALSE)
}

#' @export
AICtab <- function(...,mnames) {
    ## fancy footwork to preserve model names
    if (missing(mnames)) mnames <- get.mnames(match.call())
    ICtab(...,mnames=mnames,type="AIC")
}

#' @export
BICtab <- function(...,mnames) {
    if (missing(mnames)) mnames <- get.mnames(match.call())
    ICtab(...,mnames=mnames,type="BIC")
}

#' @export
AICctab <- function(...,mnames) {
    if (missing(mnames)) mnames <- get.mnames(match.call())
    ICtab(...,mnames=mnames,type="AICc")
}

#' @export
setGeneric("AICc", function(object, ..., nobs=NULL, k=2) standardGeneric("AICc"))

#' @export
setMethod("AICc", "mle2",
          function (object, ..., nobs, k)  {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (is.null(nobs)) {
                      nobs <- sapply(L,nobs)
                  }
                  if (length(unique(nobs))>1)
                      stop("nobs different: must have identical data for all objects")
                  logLiks <- sapply(L, logLik)
                  df <- sapply(L,attr,"df")
                  val <- -2*logLiks+k*df*(df+1)/(nobs-df-1)
                  data.frame(AICc=val,df=df)
              } else {
                  df <- attr(object,"df")
                  c(-2*logLik(object)+k*df+k*df*(df+1)/(nobs-df-1))
              }
          })

#' @export
setMethod("AICc", signature(object="logLik"),
          function(object, ..., nobs=NULL, k){
              if (missing(nobs)) {
                  if (is.null(attr(object,"nobs")))
                      stop("number of observations not specified")
                  nobs <- attr(object,"nobs")
              }
              df <- attr(object,"df")
              ## FIXME: should second "2" also be k?
              -2 * c(object) + k*df+2*df*(df+1)/(nobs-df-1)
          })

#' @rdname IC-class
#' @export
setMethod("AICc", signature(object="ANY"),
          function(object, ..., nobs=NULL, k){
              AICc(object=logLik(object, ...), nobs=nobs, k=k)
          })

#' @rdname IC-class
#' @export
setMethod("AIC", "mle2",
          function (object, ..., k = 2) {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (!all(sapply(L,class)=="mle2")) stop("all objects in list must be class mle2")
                  logLiks <- lapply(L, logLik)
                  AICs <- sapply(logLiks,AIC,k=k)
                  df <- sapply(L,attr,"df")
                  data.frame(AIC=AICs,df=df)
              } else AIC(logLik(object), k = k)
          })

### quasi- methods

#' @export
setGeneric("qAICc", function(object, ..., nobs=NULL, dispersion, k=2)
           standardGeneric("qAICc"))

#' @export
setMethod("qAICc", signature(object="ANY"),
          function(object, ..., nobs=NULL, dispersion, k=2){
              qAICc(object=logLik(object), nobs=nobs, dispersion=dispersion, k=k)
          })

#' @export
setMethod("qAICc", "mle2",
          function (object, ..., nobs, dispersion, k)  {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (missing(nobs)) {
                      nobs <- sapply(L,nobs)
                  }
                  if (missing(dispersion) && is.null(attr(object,"dispersion")))
                      stop("must specify (over)dispersion coefficient")
                  if (length(unique(nobs))>1)
                      stop("nobs different: must have identical data for all objects")
                  nobs <- nobs[1]
                  logLiks <- sapply(L, logLik)/dispersion
                  df <- sapply(L,attr,"df")+1 ## add one for scale parameter
                  val <- logLiks+k*df*(df+1)/(nobs-df-1)
                  data.frame(AICc=val,df=df)
              } else {
                  df <- attr(object,"df")
                  c(-2*logLik(object)/dispersion+2*df+2*df*(df+1)/(nobs-df-1))
              }
          })

#' @rdname IC-class
#' @export
setMethod("qAICc", signature(object="logLik"),
          function(object, ..., nobs, dispersion, k){
              if (missing(nobs)) {
                  if (is.null(attr(object,"nobs")))
                      stop("number of observations not specified")
                  nobs <- attr(object,"nobs")
              }
              if (missing(dispersion)) {
                  if (is.null(attr(object,"dispersion")))
                      stop("dispersion not specified")
                  dispersion <- attr(object,"dispersion")
              }
              df <- attr(object,"df")+1 ## add one for scale parameter
              -2 * c(object)/dispersion + k*df+2*df*(df+1)/(nobs-df-1)
          })

#' @rdname IC-class
#' @export
setGeneric("qAIC", function(object, ..., dispersion, k=2)
           standardGeneric("qAIC"))

#' @rdname IC-class
#' @export
setMethod("qAIC", signature(object="ANY"),
          function(object, ..., dispersion, k=2){
              qAIC(object=logLik(object), dispersion=dispersion, k)
          })

#' @rdname IC-class
#' @title IC methods
#' 
#' @export
setMethod("qAIC", signature(object="logLik"),
          function(object, ..., dispersion, k){
              if (missing(dispersion)) {
                  if (is.null(attr(object,"dispersion")))
                      stop("dispersion not specified")
                  dispersion <- attr(object,"dispersion")
              }
              df <- attr(object,"df")
              -2 * c(object)/dispersion + k*df
          })

#' @rdname IC-class
#' @export
setMethod("qAIC", "mle2",
          function (object, ..., dispersion, k=2) {
              L <- list(...)
              if (length(L)) {
                  L <- c(list(object),L)
                  if (!all(sapply(L,class)=="mle2"))
                      stop("all objects in list must be class mle2")
                  logLiks <- lapply(L, logLik)
                  AICs <- sapply(logLiks,qAIC, k=k, dispersion=dispersion)
                  df <- sapply(L,attr,"df")
                  data.frame(AIC=AICs,df=df)
              } else {
                  qAIC(logLik(object), k=k, dispersion=dispersion)
              }
          })

