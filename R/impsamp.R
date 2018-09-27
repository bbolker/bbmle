#' generate population prediction sample from parameters
#'
#' @param object a fitted \code{mle2} object
#' 

pop_pred <- function(object,
                     n=1000,
                     n_imp=n*10,
                     return_wts=FALSE,
                     impsamp=FALSE,
                     ginv=FALSE,
                     tol = 1e-6) {
    vv <- vcov(object)
    ## FIXME: not sure how this will interact with fixed parameters?
    Lfun <- object@minuslogl
    min_eig <- min(eigen(vv,only.values=TRUE)$values)
    if (min_eig<tol) {
        if (!use_ginv) {
            stop("non-positive definitive variance-covariance matrix ",
                 sprintf("(min eig=%f): ",min_eig),
                 "consider ginv=TRUE (and probably impsamp=TRUE)")
        }
        ## use King 1994 to 'posdefify' variance-covariance matrix
        ## (better than e.g. Matrix::nearPD)
        vv <- crossprod(as.matrix(bdsmatrix::gchol(MASS::ginv(vv))))
    }
    mv_n <- if (impsamp) n_imp else n
    mv_vals <- MASS::mvrnorm(mv_n,mu=coef(object),Sigma=vv)
    if (!(impsamp || return_wts)) return(mv_vals)
    mv_wts <- apply(mv_vals,1,Lfun)
    ## shift negative log-likelihoods (avoid underflow);
    ## find scaled likelihood
    mv_wts <- exp(-(mv_wts - min(mv_wts)))
    mv_wts <- mv_wts/sum(mv_wts)
    eff_samp <- 1/sum(mv_wts^2)  ## check ???
    attributes(mv_vals,"eff_samp") <- eff_samp
    mv_vals <- cbind(mv_vals,wts=mv_wts)
    if (return_wts) return(mv_vals)
    ## do importance sampling
    mv_vals <- mv_vals[sample(seq(nrow(mv_vals)),
                              size=n,
                              weights=mv_wts,
                              replace=TRUE),]
    return(mv_vals)
}
    
    
    
        





