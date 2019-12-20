library(bbmle)
old_opt <- options(digits=3)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)


## throws warning on some CRAN system (BDR report 2019-12-20)
## /tests-noLD/bbmle.Rcheck’
## * using R Under development (unstable) (2019-12-19 r77606)
## * using platform: x86_64-pc-linux-gnu (64-bit)
fit0 <- suppressWarnings(
    mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
         method="L-BFGS-B",lower=10)
)

stopifnot(is.na(confint(fit0)[1]))

fit1 <- mle2(y~dpois(lambda=exp(a+b*x)),start=list(a=0,b=2),data=d,
             method="L-BFGS-B",lower=-0.2)

suppressWarnings(confint(fit1))

fit2 <- mle2(y~dpois(lambda=exp(a+b*x)),start=list(a=0,b=2),data=d,
             method="L-BFGS-B")

pp <- profile(fit2,prof.lower=-0.2)
stopifnot(min(subset(as.data.frame(pp),param=="b")$par.vals.b)==-0.2)
## note that b does go below -0.2 when profiling a ...
options(old_opt)
