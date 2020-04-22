## STUB: starting to play with/think about TMB back-end
form <- quote({
    A <- exp(-alpha*t)
})
recruits ~ dnorm(settlers*A/(1+(beta*settlers*(1-A)/alpha)), sd)
TMB_mle2_function <- function(formula,
                              parameters,
                              links,
                              start,
                              parnames,
                              data=NULL) {

    ## find symbols

    

}
