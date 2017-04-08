library(bbmle)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)

## uses default parameters of LL
fit <- mle2(y~dpois(exp(loglam)),
            data=d,
           start=list(loglam=0),control=list(maxit=2))
pp <- profile(fit)
