## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- bt, echo=TRUE------------------------------------------------------
# ref BT.RES in Nash and Walker-Smith (1987)
library(lbfgsb3c)
sessionInfo()

bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

bt.badsetup<-function(n){
   x<-rep(0,n)
   lo<-rep(0,n)
   up<-lo # to get arrays set
   bmsk<-rep(1,n)
   bmsk[(trunc(n/2)+1)]<-0
   for (i in 1:n) { 
      x[i]<-2.2*i-n
      lo[i]<-1.0*(i-1)*(n-1)/n
      up[i]<-1.0*i*(n+1)/n
   }
   result<-list(x=x, lower=lo, upper=up, bdmsk=bmsk)
}

bt.setup0<-function(n){
   x<-rep(0,n)
   lo<-rep(0,n)
   up<-lower # to get arrays set
   bmsk<-rep(1,n)
   bmsk[(trunc(n/2)+1)]<-0
   for (i in 1:n) { 
      lo[i]<-1.0*(i-1)*(n-1)/n
      up[i]<-1.0*i*(n+1)/n
   }
   x<-0.5*(lo+up)
   result<-list(x=x, lower=lo, upper=up, bdmsk=bmsk)
}
nn <- 4
baddy <- bt.badsetup(nn)
lo <- baddy$lower
up <- baddy$upper
x0 <- baddy$x
baddy
## optim()
solbad0 <- optim(x0, bt.f, bt.g, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=3))
solbad0
## lbfgsb3c
solbad1 <- lbfgsb3(x0, bt.f, bt.g, lower=lo, upper=up, control=list(trace=3))
solbad1
## Possible timings
## library(microbenchmark)
## tbad0 <- microbenchmark(optim(x0, bt.f, bt.g, lower=lo, upper=up, method="L-BFGS-B"))
## t3c <- microbenchmark(lbfgsb3(x0, bt.f, bt.g, lower=lo, upper=up))
## Via optimx package
## library(optimx)
## meths <- c("L-BFGS-B", "lbfgsb3") # Note: lbfgsb3c not yet in optimx on CRAN
## allbt0 <- opm(x0, bt.f, bt.g, lower=lo, upper=up, method=meths)
## summary(allbt0, order=value)


## ---- candlestick--------------------------------------------------------
# candlestick function
# J C Nash 2011-2-3
cstick.f<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-crossprod(x)
  f<-as.double(r2+alpha/r2)
  return(f)
}

cstick.g<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-as.numeric(crossprod(x))
  g1<-2*x
  g2 <- (-alpha)*2*x/(r2*r2)
  g<-as.double(g1+g2)
  return(g)
}
library(lbfgsb3c)
nn <- 2
x0 <- c(10,10)
lo <- c(1, 1)
up <- c(10,10)
print(x0)
c2o <- optim(x0, cstick.f, cstick.g, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=3))
c2o
c2l <- lbfgsb3(x0, cstick.f, cstick.g, lower=lo, upper=up, control=list(trace=3))
c2l

## meths <- c("L-BFGS-B", "lbfgsb3c", "Rvmmin", "Rcgmin", "Rtnmin")
## require(optimx)

## cstick2a <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2a, par.select=1:2, order=value))
lo <- c(4, 4)
c2ob <- optim(x0, cstick.f, cstick.g, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=3))
c2ob
c2lb <- lbfgsb3(x0, cstick.f, cstick.g, lower=lo, upper=up, control=list(trace=3))
c2lb




## cstick2b <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2b, par.select=1:2, order=value))

nn <- 100
x0 <- rep(10, nn)
up <- rep(10, nn)
lo <- rep(1e-4, nn)
cco <- optim(x0, cstick.f, cstick.g, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=3))
cco
ccl <- lbfgsb3(x0, cstick.f, cstick.g, lower=lo, upper=up, control=list(trace=3))
ccl
## cstickc0 <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstickc0, par.select=1:5, order=value))
## lo <- rep(1, nn)
## cstickca <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstickca, par.select=1:5, order=value))
## lo <- rep(4, nn)
## cstickcb <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstickcb, par.select=1:5, order=value))

## ----exrosen-------------------------------------------------------------
# require(funconstrain) ## not in CRAN, so explicit inclusion of this function
# exrosen <- ex_rosen()
# exrosenf <- exrosen$fn
exrosenf <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    fsum <- 0
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        f_p1 <- 10 * (par[p2] - par[p1]^2)
        f_p2 <- 1 - par[p1]
        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
    }
    fsum
}
# exroseng <- exrosen$gr
exroseng <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    grad <- rep(0, n)
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]
        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]
        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
    }
    grad
}

exrosenx0 <- function (n = 20) {
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    rep(c(-1.2, 1), n/2)
}

## meths <- c("L-BFGS-B", "lbfgsb3c", "Rvmmin", "Rcgmin", "Rtnmin")
## require(optimx)
for (n in seq(2,12, by=2)) {
  cat("ex_rosen try for n=",n,"\n")
  x0 <- exrosenx0(n)
  lo <- rep(.5, n)
  up <- rep(3, n)
  print(x0)
  eo <- optim(x0, exrosenf, exroseng, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=3))
  eo
  el <- lbfgsb3(x0, exrosenf, exroseng, lower=lo, upper=up, control=list(trace=3))
  el
##   erfg <- opm(x0, exrosenf, exroseng, method=meths, lower=lo, upper=up)
##   print(summary(erfg, par.select=1:2, order=value))
}



## ---- usingFortran, eval=FALSE-------------------------------------------
#  system("cd ~/temp")
#  system("R CMD SHLIB jrosen.f")
#  dyn.load("jrosen.so")
#  is.loaded("rosen")
#  x0 <- as.double(c(-1.2,1))
#  fv <- as.double(-999)
#  n <- as.double(2)
#  testf <- .Fortran("rosen", n=as.integer(n), x=as.double(x0), fval=as.double(fv))
#  testf
#  
#  rrosen <- function(x) {
#    fval <- 0.0
#    for (i in 1:(n-1)) {
#      dx <- x[i + 1] - x[i] * x[i]
#      fval <- fval + 100.0 * dx * dx
#      dx <- 1.0 - x[i]
#      fval <- fval + dx * dx
#    }
#    fval
#  }
#  
#  (rrosen(x0))
#  
#  frosen <- function(x){
#      nn <- length(x)
#      if (nn > 100) { stop("max number of parameters is 100")}
#      fv <- -999.0
#      val <- .Fortran("rosen", n=as.integer(nn), x=as.double(x), fval=as.double(fv))
#      val$fval # NOTE--need ONLY function value returned
#  }
#  # Test the funcion
#  tval <- frosen(x0)
#  str(tval)
#  
#  mynm <- optim(x0, rrosen, control=list(trace=1))
#  mynm
#  mynmf <- optim(x0, frosen, control=list(trace=1))
#  mynmf
#  
#  library(lbfgsb3c)
#  cat("try min")
#  myopR <- lbfgsb3c(x0, rrosen, gr=NULL, control=list(trace=3))
#  myopR
#  myop <- lbfgsb3c(x0, frosen, gr=NULL, control=list(trace=3))
#  myop

