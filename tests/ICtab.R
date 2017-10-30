library(bbmle)

set.seed(101)
z = rpois(100,lambda=5)

m1 = mle2(z~dpois(lambda=L),start=list(L=4),data=data.frame(z))

ICtab(m1,type="qAICc",dispersion=1.2,nobs=100)

m2 = glm(z~1,family=poisson)
qAICc(m2,nobs=100,dispersion=2)

## test that dAIC ignores
m3 <- glm(z~1,family=quasipoisson)
aa <- AICtab(m1,m2,m3,weights=TRUE)
stopifnot(any(!is.na(aa$dAIC)),
          any(!is.na(aa$weight)))

set.seed(101)
x <- rnorm(100)
dd <- data.frame(y=rnorm(100,2+3*x,sd=1),x)
m4A <- lm(y~x,dd)
m4B <- mle2(y~dnorm(mean=a+b*x,sd=exp(logsd)),
           data=dd,
           start=list(a=1,b=1,logsd=0))
## cosmetic differences only
stopifnot(all.equal(AIC(m4A,m4B)[,"AIC"],
                    AIC(m4B,m4A)[,"AIC"]))

