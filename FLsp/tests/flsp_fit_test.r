# Try different optimisers
setwd("c:/Projects/Deepfishman/deepfishman")

library(FLsp)

# Try SAA
data(nzrl)
data <- nzrl

p <- 1
# Polacheck's results
# NZRL
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207

# what do we get?
sp <- eval_FLsp(data$catch, data$cpue, r, k)

# Not bad!
sp[["qhat"]]

# Just the likelihood
ll <- get_logl(c(r,k),data$catch, data$cpue)

# How do we deal with NA and Infs...
ll_objfun <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  ll <- get_logl(c(r,k),catch, cpue)
#  res <- .Call("flspCpp",catch,cpue,r,1,k)
  #if (is.nan(res[["ll"]])) res[["ll"]] <- NA
  #if (is.na(res[["ll"]])) res[["ll"]] <- Inf
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  cat("ll: ", -ll, "\n")
  return(-ll)
}

ll_gradfun <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  #res <- .Call("flspCpp",catch,cpue,r,1,k)
  res <- eval_FLsp(catch, cpue, r, k)
  cat("ll_grad_r: ", -res[["ll_grads"]][1], "\n")
  cat("ll_grad_k: ", -res[["ll_grads"]][2], "\n")
  return(c(-res[["ll_grads"]][1],-res[["ll_grads"]][2]))
}

# Give the objective functions a go
ll_objfun(log(c(r,k)),data$catch,data$cpue)
ll_gradfun(log(c(r,k)),data$catch,data$cpue)

#*******************************************************************************
# Can't pass additional arguments to it
# gsl
library(gsl)

# How do we deal with NA and Infs...
ll_objfun_gsl <- function(params)
{
  r <- exp(params[1])
  k <- exp(params[2])
  ll <- get_logl(c(r,k),catch, cpue, p=1, extinct_val=1e-9)
#  res <- .Call("flspCpp",catch,cpue,r,1,k)
  if (is.nan(ll)) ll <- NA
  #if (is.na(res[["ll"]])) res[["ll"]] <- Inf
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  cat("ll: ", -ll, "\n")
  return(-ll)
}

ll_gradfun_gsl <- function(params)
{
  r <- exp(params[1])
  k <- exp(params[2])
  #res <- .Call("flspCpp",catch,cpue,r,1,k)
  res <- eval_FLsp(catch, cpue, r, k, p=1, extinct_val=1e-9)
  cat("ll_grad_r: ", -res[["ll_grads"]][1], "\n")
  cat("ll_grad_k: ", -res[["ll_grads"]][2], "\n")
  return(c(-res[["ll_grads"]][1],-res[["ll_grads"]][2]))
}

catch <- data$catch
cpue <- data$cpue

# bfgs doesn't work
init_params <- log(c(0.1,100000))
gsl_test <- multimin(x = init_params, f=ll_objfun_gsl, df=ll_gradfun_gsl, method="bfgs")

# doesn't work conjugate-fr
init_params <- log(c(0.1,100000))
gsl_test <- multimin(x = init_params, f=ll_objfun_gsl, df=ll_gradfun_gsl, method="conjugate-fr")

# conjugate-pr, doesn't work
init_params <- log(c(0.1,100000))
gsl_test <- multimin(x = init_params, f=ll_objfun_gsl, df=ll_gradfun_gsl, method="conjugate-pr")

# steepest descent
init_params <- log(c(0.1,100000))
gsl_test <- multimin(x = init_params, f=ll_objfun_gsl, df=ll_gradfun_gsl, method="steepest-descent")



# example
x0 <- c(-1.2, 1)
f <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
df <- function(x) c(-2*(1 - x[1]) + 100 * 2 * (x[2] - x[1]^2) * (-2*x[1]),
100 * 2 * (x[2] - x[1]^2))
# The simple way to call multimin.
state <- multimin(x0, f, df)



#*******************************************************************************
#test <- optim(c(1,1000),fn=ll_objfun,gr=ll_gradfun,method="L-BFGS-B", lower=c(1e-6,1e-6), upper=c(Inf,Inf),catch=saa$catch,cpue=saa$cpue)

test <- optim(log(c(0.5,100000)),fn=ll_objfun,gr=ll_gradfun,method="BFGS", catch=data$catch,cpue=data$cpue)

test <- optim(log(c(0.5,100000)),fn=ll_objfun,gr=ll_gradfun,method="CG", catch=data$catch,cpue=data$cpue)


test <- optim(log(c(0.1,200)),fn=ll_objfun,gr=ll_gradfun,method="BFGS", catch=data$catch/1000,cpue=data$cpue)
# Don't know what happens here...

# Try with parscale - shite
test <- optim(log(c(0.5,100000)),fn=ll_objfun,gr=ll_gradfun,method="BFGS", catch=data$catch,cpue=data$cpue,control=list(parscale=c(1000,1)))
# Without gradient function - even worse
test <- optim(log(c(0.5,100000)),fn=ll_objfun,method="BFGS", catch=data$catch,cpue=data$cpue,control=list(parscale=c(1000,1)))
test <- optim(log(c(0.5,100000)),fn=ll_objfun,gr=ll_gradfun,method="CG", catch=data$catch,cpue=data$cpue,control=list(parscale=c(1000,1)))

# Forget optim - shite
# Try DEoptim
library(DEoptim)

ll_objfun <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  if (is.nan(res[["ll"]])) res[["ll"]] <- -Inf
  if (is.na(res[["ll"]])) res[["ll"]] <- -Inf
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  cat("ll: ", -res[["ll"]], "\n")
  return(-res[["ll"]])
}

#test <- DEoptim(fn=ll_objfun,lower=log(c(1e-4,10000)),upper=log(c(10,1e7)),control=DEoptim.control(trace=TRUE),catch=data$catch,cpue=data$cpue)
testde <- DEoptim(fn=ll_objfun,lower=log(c(1e-9,1e-9)),upper=log(c(1e9,1e9)),control=DEoptim.control(trace=TRUE),catch=data$catch,cpue=data$cpue)
exp(testde$optim$bestmem)
r <- exp(testde$optim$bestmem)[1]
k <- exp(testde$optim$bestmem)[2]
# is ll better than that reported?
deres <- .Call("flspCpp",data$catch,data$cpue,r,1,k)
deres[["ll"]]
sp[["ll"]]
# It's higher so better (trying to maximise it)

# Works fine
# and can include bounds

# Try nls.lm
library(minpack.lm)

ll_objfun_nls <- function(params,catch,cpue)
{
  r <- (params[["r"]])
  k <- (params[["k"]])
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  residuals <- res[["res"]]
  return(residuals)
}

ll_objjac_nls <- function(params,catch,cpue)
{
    cat("Calling jacobian\n")
  r <- (params[["r"]])
  k <- (params[["k"]])
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  jacr <- res[["res_grad_r"]]
  jack <- res[["res_grad_k"]]
  return(c(jacr,jack))
}

testnls <- nls.lm(par = list(r=(1), k= (100000)), fn=ll_objfun_nls, catch=data$catch,cpue=data$cpue)
r <- (testnls$par$r)
k <- (testnls$par$k)

# with jacobian
testnls_jac <- nls.lm(par = list(r=(1), k= (100000)), fn=ll_objfun_nls, jac=ll_objjac_nls,catch=data$catch,cpue=data$cpue)

# slightly different. Which one has better residuals and likelihood
nlsres <- .Call("flspCpp",data$catch,data$cpue,testnls$par$r,1,testnls$par$k)
nlsjacres <- .Call("flspCpp",data$catch,data$cpue,testnls_jac$par$r,1,testnls_jac$par$k)
nlsres[["ll"]]
nlsjacres[["ll"]]
sp[["ll"]]

sum(nlsres[["res"]]^2)
sum(nlsjacres[["res"]]^2)




# Conclusion
# optim doesn't work at all. Even with gradients and parscaling

# DEoptim works very well. Fast, has bounds, better logl than reported
# But cannot use gradient function
# Could use it to get the hessian iside C?

# nls.lm seems to work but logl is much worse. That is because we aren't
# fitting the likelihood, we are getting sum sq (I - Ihat)
# so imposes error structure? not lognormal? Need to be careful...
# Maybe the residuals we return are taken from the logl equation?

#detach("package:FLsp")
#dyn.unload("C:/R/library/FLsp/libs/i386/FLsp.dll")

#*****************************************************************************
# Back to NZRL
# Are reported values at min ll

p <- 1
# Polacheck's results
# NZRL
data <- read.csv("FLsp/data/NZRL.csv")
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207
s2 <- sigma^2
# what do we get?
sp <- .Call("flspCpp",data$catch,data$cpue,r,p,k)

#1D profile
rvec <- seq(from=0.5, to =1.5, length=50) * r
kvec <- seq(from=0.5, to =1.5, length=50) * k
llrvec <- rep(NA,length(rvec))
llgradrrvec <- rep(NA,length(rvec))
llgradrkvec <- rep(NA,length(rvec))
llkvec <- rep(NA,length(rvec))
llgradkrvec <- rep(NA,length(rvec))
llgradkkvec <- rep(NA,length(rvec))

for (i in 1:length(rvec))
{
    rres <- .Call("flspCpp",data$catch,data$cpue,rvec[i],p,k)
    llrvec[i] <- rres[["ll"]]
    llgradrrvec[i] <- rres[["ll_grad_r"]]
    llgradrkvec[i] <- rres[["ll_grad_k"]]


    kres <- .Call("flspCpp",data$catch,data$cpue,r,p,kvec[i])
    llkvec[i] <- kres[["ll"]]
    llgradkrvec[i] <- kres[["ll_grad_r"]]
    llgradkkvec[i] <- kres[["ll_grad_k"]]
}

par(mfrow=c(2,1))
plot(rvec,llrvec,type="l")
lines(c(r,r),c(-1e-6,1e6))
plot(kvec,llkvec,type="l")
lines(c(k,k),c(-1e-6,1e6))
# Likelihood looks OK - so why the problem


par(mfrow=c(3,2))
plot(rvec,llrvec,type="l")
lines(c(r,r),c(-1e-6,1e6))
plot(kvec,llkvec,type="l")
lines(c(k,k),c(-1e-6,1e6))
# Plot grads
plot(rvec,llgradrrvec,type="l")
lines(c(r,r),c(-1e-6,1e6))
lines(c(0,1e6),c(0,0))
plot(kvec,llgradkrvec,type="l")
lines(c(k,k),c(-1e-6,1e6))
lines(c(0,1e6),c(0,0))
plot(rvec,llgradrkvec,type="l")
lines(c(r,r),c(-1e-6,1e6))
lines(c(0,1e6),c(0,0))
plot(kvec,llgradkkvec,type="l")
lines(c(k,k),c(-1e-6,1e6))
lines(c(0,1e6),c(0,0))
# gradient is being returned correctly, but grad = 0 as r or k increase

#************************************************************************
# what happens if k too low, i.e. Biomass crashes.
setwd("~/Work/deepfishman/trunk")
#setwd("m:/Projects/Deepfishman/deepfishman")

library(FLsp)
data <- read.csv("FLsp/data/NZRL.csv")
#data <- read.csv("FLsp/data/NNH.csv")
#data <- read.csv("FLsp/data/SAAL.csv")

p <- 1
# Polacheck's results
# NZRL
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207
# what do we get?
sp <- .Call("flspCpp",data$catch,data$cpue,r,p,k)
sp[["qhat"]]
sp[["B"]]

sp <- .Call("flspCpp",data$catch,data$cpue,r,p,500)
sp[["B"]]
# goes negative.
# stop that!




#*****************************************************************************
# Testing the LL equation
p <- 1
# Polacheck's results
# NZRL
data <- read.csv("FLsp/data/NZRL.csv")
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207
s2 <- sigma^2
# what do we get?
sp <- .Call("flspCpp",data$catch,data$cpue,r,p,k)
sp[["qhat"]]
# Not bad!
# Assume that q is OK.
# and Ihat is OK

vhat <- log(data$cpue) - log(sp[["Ihat"]])
n <- sum(!is.na(data$cpue))
sigma2 <- sum(vhat^2)/n
sp[["sigma2"]]
# sigma2 is spot on

l <- prod(exp(-(vhat^2)/(2*sigma2)) / sqrt(2*pi*sigma2))
ll <- log(l)
sp[["ll"]]
# ll is spot on

# Repeat with other data sets
# NNH
data <- read.csv("FLsp/data/NNH.csv")
r <- 0.379
k <- 2772.6
q <- 4.350e-4
sigma <- 0.124
s2 <- sigma^2
# what do we get?
sp <- .Call("flspCpp",data$catch,data$cpue,r,p,k)
sp[["qhat"]]
# Not bad!
# Assume that q is OK.
# and Ihat is OK

vhat <- log(data$cpue) - log(sp[["Ihat"]])
n <- sum(!is.na(data$cpue))
sigma2 <- sum(vhat^2)/n
sp[["sigma2"]]
s2
# sigma2 is spot on

l <- prod(exp(-(vhat^2)/(2*sigma2)) / sqrt(2*pi*sigma2))
ll <- log(l)
sp[["ll"]]
# ll is spot on


# Repeat with other data sets
# SAA
data <- read.csv("FLsp/data/SAA.csv")
r <- 0.328
k <- 239.6
q <- 26.71e-4 # order looks wrong
sigma <- 0.111
s2 <- sigma^2
# what do we get?
sp <- .Call("flspCpp",data$catch,data$cpue,r,p,k)
sp[["qhat"]]
# Not bad!
# Assume that q is OK.
# and Ihat is OK

vhat <- log(data$cpue) - log(sp[["Ihat"]])
n <- sum(!is.na(data$cpue))
sigma2 <- sum(vhat^2)/n
sp[["sigma2"]]
s2
# sigma2 is spot on

l <- prod(exp(-(vhat^2)/(2*sigma2)) / sqrt(2*pi*sigma2))
ll <- log(l)
sp[["ll"]]
# ll is spot on


