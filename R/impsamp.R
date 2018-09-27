#' generate population prediction sample from parameters
#'
#' @param object a fitted \code{mle2} object
#' @param n number of samples to return
#' @param n_imp number of total samples from which to draw, if doing importance sampling
#' @param return_wts return a column giving the weights of the samples, for use in weighted summaries?
#' @param impsamp subsample values (with replacement) based on their weights?
#' @param ginv use Gill and King generalized-inverse procedure to correct non-positive-definite variance-covariance matrix if necessary?
#' @param tol tolerance for detecting small eigenvalues
#' 
#' @references Gill, Jeff, and Gary King. "What to Do When Your Hessian Is Not Invertible: Alternatives to Model Respecification in Nonlinear Estimation." Sociological Methods & Research 33, no. 1 (2004): 54-87.

pop_pred_samp <- function(object,
                     n=1000,
                     n_imp=n*10,
                     return_wts=FALSE,
                     impsamp=FALSE,
                     ginv=FALSE,
                     tol = 1e-6) {
    vv <- vcov(object)
    cc <- coef(object)
    ## FIXME: not sure how this will interact with fixed parameters?
    Lfun <- object@minuslogl
    min_eig <- min(eigen(vv,only.values=TRUE)$values)
    if (min_eig<tol) {
        if (!ginv) {
            stop("non-positive definitive variance-covariance matrix ",
                 sprintf("(min eig=%f): ",min_eig),
                 "consider ginv=TRUE (and probably impsamp=TRUE)")
        }
        ## use King 1994 to 'posdefify' variance-covariance matrix
        ## (better than e.g. Matrix::nearPD)
        vv <- crossprod(as.matrix(bdsmatrix::gchol(MASS::ginv(vv))))
    }
    mv_n <- if (impsamp) n_imp else n
    mv_vals <- MASS::mvrnorm(mv_n,mu=cc,Sigma=vv)
    if (!(impsamp || return_wts)) return(mv_vals)
    ## compute MV sampling probabilities
    mv_wts <- dmvnorm(mv_vals,mu=cc,Sigma=vv,log=TRUE)
    ## compute likelihoods of each sample point
    L_wts <- apply(mv_vals,1,Lfun)
    ## shift negative log-likelihoods (avoid underflow);
    ## find scaled likelihood
    L_wts <- L_wts - mv_wts ## subtract log samp prob
    L_wts <- exp(-(L_wts - min(L_wts)))
    L_wts <- L_wts/sum(L_wts)
    eff_samp <- 1/sum(L_wts^2)  ## check ???
    attr(mv_vals,"eff_samp") <- eff_samp
    mv_vals <- cbind(mv_vals,wts=L_wts)
    if (return_wts) return(mv_vals)
    ## do importance sampling
    mv_vals <- mv_vals[sample(seq(nrow(mv_vals)),
                              size=n,
                              prob=L_wts,
                              replace=TRUE),]
    return(mv_vals)
}
    
## copy (!!) dmvnorm from emdbook to avoid cyclic dependency problems
dmvnorm <- function (x, mu, Sigma, log = FALSE, tol = 1e-06) {
    if (is.vector(x)) 
        x = t(as.matrix(x))
    n = length(mu)
    if (is.vector(mu)) {
        p <- length(mu)
        if (is.matrix(x)) {
            mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
        }
    }
    else {
        p <- ncol(mu)
    }
    if (!all(dim(Sigma) == c(p, p)) || nrow(x) != nrow(mu)) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE) 
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1]))) 
        stop("Sigma is not positive definite")
    z = t(x - mu)
    logdetS = try(determinant(Sigma, logarithm = TRUE)$modulus)
    attributes(logdetS) <- NULL
    iS = try(solve(Sigma))
    if (class(iS) == "try-error" || class(logdetS) == "try-error") {
        warning("difficulty inverting/taking determinant of Var-Cov matrix")
        return(NA)
    }
    ssq = diag(t(z) %*% iS %*% z)
    loglik = -(n * (log(2*pi)) +  logdetS + ssq)/2
    if (log) loglik else exp(loglik)
  }
