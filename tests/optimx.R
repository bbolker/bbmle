library(bbmle)
old_opt <- options(digits=3)
if (require(optimx)) {
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)

## breaks, don't try this
## optimx(fn=Lfn,par=c(15,6),method="Rvmmin")

suppressWarnings(m1 <- mle2(minuslogl=y~dpois(lambda=ymax/(1+x/xhalf)),
     start=list(ymax=15,xhalf=6),data=d,
     optimizer="optimx",
           method=c("BFGS","Nelder-Mead","CG")))

## FIXME!! fails (although not with an error, because
##  errors are caught by profiling) due to npar now
## being restricted to >1 in optimx 2012.05.24 ...

suppressWarnings(head(as.data.frame(profile(m1))))

## GH #31: optimx() can silently return NA parameter estimates when a
## method's backing package isn't installed (e.g. bobyqa/newuoa need
## minqa, spg needs BB, ucminf needs ucminf, nmkb/hjkb need dfoptim).
## mle2(optimizer="optimx") now checks up front: an unrecognized method
## name errors immediately, and a recognized method backed by a missing
## package errors with a clear message instead of quietly propagating NAs.
zz <- rpois(100,lambda=5)
dd <- data.frame(z=zz)

## unrecognized/mistyped method name (this was the original report's
## own reproducible example: passing the *package* name "minqa" instead
## of an actual optimx method such as "bobyqa")
stopifnot(inherits(try(mle2(z~dpois(lambda=L), start=list(L=4), data=dd,
                             optimizer="optimx", method="minqa"),
                        silent=TRUE),
                    "try-error"))

## sanity-check the method -> required-package mapping bbmle relies on
## (taken from optimx's own registry, optimx::ctrldefault())
pkg_for_method <- function(m) {
    ctrl <- optimx::ctrldefault(1)
    ctrl$allpkg[match(m, ctrl$allmeth)]
}
stopifnot(pkg_for_method("bobyqa")=="minqa",
          pkg_for_method("newuoa")=="minqa",
          pkg_for_method("spg")=="BB",
          pkg_for_method("ucminf")=="ucminf",
          pkg_for_method("nmkb")=="dfoptim",
          pkg_for_method("hjkb")=="dfoptim")

## when the backing package IS installed, the method should work as
## normal (regression check); if it genuinely isn't installed here,
## requesting it should fail with an informative error rather than
## silently returning NAs
if (requireNamespace("minqa", quietly=TRUE)) {
    m2 <- mle2(z~dpois(lambda=L), start=list(L=4), data=dd,
               optimizer="optimx", method="bobyqa")
    stopifnot(is.finite(coef(m2)))
} else {
    stopifnot(inherits(try(mle2(z~dpois(lambda=L), start=list(L=4), data=dd,
                                 optimizer="optimx", method="bobyqa"),
                            silent=TRUE),
                        "try-error"))
}

detach("package:optimx")
}

## GH #31: generic guard against a silent NA-parameter result from *any*
## optimizer (not just optimx) -- doesn't require optimx/minqa at all
badopt <- function(par, fn, ...) list(par=c(L=NA_real_), value=NA_real_, convergence=1)
zz2 <- rpois(100,lambda=5)
stopifnot(inherits(try(mle2(z~dpois(lambda=L), start=list(L=4),
                             data=data.frame(z=zz2),
                             optimizer="user", optimfun=badopt),
                        silent=TRUE),
                    "try-error"))

options(old_opt)
