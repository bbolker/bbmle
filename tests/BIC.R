require(bbmle)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
fit <- mle2(y~dpois(lambda=ymax/(1+x/xhalf)), start=list(ymax=25,xhalf=3),data=d)
fit2 <- mle2(y~dpois(lambda=(x+1)*slope), start=list(slope=1),data=d)
BIC(fit)
BIC(fit,fit2)

## GH #35: nobs() method was missing for mle2 objects, so BICtab()/ICtab()
## failed outright (fell back to nobs.default -> "no 'nobs' method is
## available"), and even when it worked, an explicit 'nobs' argument to
## ICtab(type="BIC") was silently ignored (only AICc/qAICc used it).
stopifnot(identical(nobs(fit), 11L))

bictab1 <- BICtab(fit, fit2, mnames=c("fit","fit2"), base=TRUE)
stopifnot(isTRUE(all.equal(bictab1$BIC, c(BIC(fit), BIC(fit2)))))

## explicit nobs override should change the BIC value (e.g. for binomial
## data where the correct sample size is the number of trials, not rows)
bictab2 <- BICtab(fit, mnames="fit", base=TRUE, nobs=1000)
expected <- -2*c(logLik(fit)) + attr(logLik(fit),"df")*log(1000)
stopifnot(isTRUE(all.equal(bictab2$BIC, expected)))
stopifnot(!isTRUE(all.equal(bictab2$BIC, bictab1$BIC[1])))
