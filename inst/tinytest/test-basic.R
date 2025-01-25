library(bbmle)
d <- data.frame(x=0:10,
                y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))

maxit <- 1000
fit <- mle2(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
            start=list(lymax=0,lhalf=0),
            data=d,
            control=list(maxit=maxit),
            parameters=list(lymax~1,lhalf~1))
expect_equal(coef(fit),
             c(lymax = 3.21885266087612, lhalf = 1.11703461971453))
