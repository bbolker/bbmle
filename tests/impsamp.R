library(bbmle)
set.seed(1002)
lymax <- c(0,2)
lhalf <- 0
x <- runif(200)
g <- factor(rep(c("a","b"),each=100))
y <- rnbinom(200,mu=(exp(lymax[g])/(1+x/exp(lhalf)))^2,size=2)
dd <- data.frame(x,g,y)

fit3 <- mle2(y~dnbinom(mu=(exp(lymax)/(1+x/exp(lhalf)))^d,size=exp(logk)),
    parameters=list(lymax~g),
    start=list(lymax=0,lhalf=0,logk=0,d=NA),
    data=dd,
    fixed=list(d=2))

pp <- pop_pred_samp(fit3,PDify=TRUE)
stopifnot(
    !any(is.na(pp)),
    identical(colnames(pp),
              c("lymax.(Intercept)", "lymax.gb", "lhalf", "logk", "d")))
