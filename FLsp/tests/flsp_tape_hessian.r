# Need to add index years checker in c code

# FLsp test

#setwd("~/Work/deepfishman/trunk")
setwd("m:/Projects/Deepfishman/deepfishman")

library(FLsp)
#dyn.load("C:/R/library/FLsp/libs/i386/FLsp.dll")

# Try SAA

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
tape <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)

sp[["B"]]
tape[["B"]]

sp[["Ihat"]]
tape[["Ihat"]]

sp[["qhat"]]
tape[["qhat"]]

sp[["sigma2"]]
tape[["sigma2"]]

sp[["ll"]]
tape[["ll"]]

sp[["ll_grad_r"]]
tape[["ll_grad_r"]]

sp[["ll_grad_k"]]
tape[["ll_grad_k"]]

tape[["hessian"]]
# Only computes the lower half
# For H[i][j]
# Only j > i values unallocated

# Check gradients
tiny <- 1e-9
llorig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)[["ll"]]
rtiny <- r + tiny
llrtiny <- .Call("flspCpp_tape",data$catch,data$cpue,rtiny,p,k)[["ll"]]
dlldr <- (llrtiny - llorig) / tiny
ktiny <- k + tiny
llktiny <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,ktiny)[["ll"]]
dlldk <- (llktiny - llorig) / tiny

simple_diff <- function(r,k)
{
	llorig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}



simple_diff(r,k)

# Approximate hessian
tiny <- 1e-9
dorig <- simple_diff(r,k)
drtiny <- simple_diff(r+tiny,k)
dll2dr2 <- (drtiny[1]-dorig[1])/tiny
dktiny <- simple_diff(r,k+tiny)
dll2dk2 <- (dktiny[2]-dorig[2])/tiny
dll2drk <-(dktiny[1]-dorig[1])/tiny
dll2dkr <-(drtiny[2]-dorig[2])/tiny
# My approximate Hessian looks pretty damn good
#*******************************************************************************
# Fit with DEoptim
library(DEoptim)

ll_objfun <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  #r <- exp(params[1])
  #k <- exp(params[2])
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  if (is.nan(res[["ll"]])) res[["ll"]] <- -Inf
  if (is.na(res[["ll"]])) res[["ll"]] <- -Inf
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  cat("ll: ", -res[["ll"]], "\n")
  return(-res[["ll"]])
}

testde <- DEoptim(fn=ll_objfun,lower=log(c(1e-9,1e-9)),upper=log(c(1e9,1e9)),control=DEoptim.control(trace=TRUE),catch=data$catch,cpue=data$cpue)
#testde <- DEoptim(fn=ll_objfun,lower=(c(1e-9,1e-9)),upper=(c(1e9,1e9)),control=DEoptim.control(trace=TRUE),catch=data$catch,cpue=data$cpue)
exp(testde$optim$bestmem)
r <- exp(testde$optim$bestmem)[1]
k <- exp(testde$optim$bestmem)[2]
# Get the hessian at this point
hess <- .Call("flspCpp_tape",data$catch,data$cpue,r,1,k)[["hessian"]]

hess[1,2] <- hess[2,1]

#*******************************************************************************
# Uncertainty stuff from wiki
# Estimate of the Hessian (matrix of second-order partial derivatives) at the MLE
# which approximates the (multivariate normal) variance-covariance matrix of the MLEs,
# from which we can derive standard errors and correlation:

# What is approximate hessian from optim
#test <- optim(log(c(r,k)),fn=ll_objfun,gr=ll_gradfun,method="BFGS", hessian = TRUE, catch=data$catch,cpue=data$cpue)
#test$hessian
#tape[["hessian"]]
# very different. optim approx is awful (see approximation above)

mu <- c(r,k)
vcov.matrix <- hess
# Tape only returns bottom left matrix. top right = bottom left
## derive the correlation matrix from covariance matrix
cor.matrix <- vcov.matrix
diag(cor.matrix)[] <- 1
cor.matrix[1,2] <- cor.matrix[1,2] / sqrt(prod(diag(vcov.matrix)))
cor.matrix[2,1] <- cor.matrix[1,2]
## the cholesky decomposition of the co matrix
#Psi <- chol(cor.matrix)
Psi <- t(chol(cor.matrix))

# Reconstruct covariance matrix (why?)
vcov.matrix2 <- vcov.matrix
vcov.matrix2[1,2] <- cor.matrix[1,2] * sqrt(prod(vcov.matrix))
vcov.matrix2[2,1] <- cor.matrix[2,1] * sqrt(prod(vcov.matrix))

## Individual SD of parameters
param.sd <- sqrt(diag(vcov.matrix))
# Need to be +ve?
param.sd <- sqrt(abs(diag(vcov.matrix)))

## number of Monte Carlo samples
nsam <- 1000

## set up theta matrix of samples
theta <- matrix(nrow=nsam,ncol=length(param.sd))
theta[,1] <- rnorm(nsam)
theta[,2] <- rnorm(nsam)

## use apply to "entrain" the correlation into the samples
theta <- t(apply(theta,1,function(x,Psi){x <- Psi%*%x},Psi))

## rescale to the appropriate mean and variance
theta[,1] <- r+theta[,1]*param.sd[1]
theta[,2] <- k+theta[,2]*param.sd[2]
colnames(theta) <- c("r","k")

# check
cor(theta)
cor.matrix

cov(theta)
vcov.matrix

# Actually not bad but what is going with theta?

#*******************************************************************************
# Fitting with optim
ll_objfun_tape <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  res <- .Call("flspCpp_tape",catch,cpue,r,1,k)
  if (is.nan(res[["ll"]])) res[["ll"]] <- NA
  #if (is.na(res[["ll"]])) res[["ll"]] <- Inf
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  cat("ll: ", -res[["ll"]], "\n")
  return(-res[["ll"]])
}

ll_gradfun_tape <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  res <- .Call("flspCpp_tape",catch,cpue,r,1,k)
  cat("ll_grad_r: ", -res[["ll_grad_r"]], "\n")
  cat("ll_grad_k: ", -res[["ll_grad_k"]], "\n")
  return(c(-res[["ll_grad_r"]],-res[["ll_grad_k"]]))
}


ll_objfun <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  if (is.nan(res[["ll"]])) res[["ll"]] <- NA
  #if (is.na(res[["ll"]])) res[["ll"]] <- Inf
  cat("r: ", r, "\n")
  cat("k: ", k, "\n")
  cat("ll: ", -res[["ll"]], "\n")
  return(-res[["ll"]])
}

ll_gradfun <- function(params,catch,cpue)
{
  r <- exp(params[1])
  k <- exp(params[2])
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  cat("ll_grad_r: ", -res[["ll_grad_r"]], "\n")
  cat("ll_grad_k: ", -res[["ll_grad_k"]], "\n")
  return(c(-res[["ll_grad_r"]],-res[["ll_grad_k"]]))
}



# BFGS - shite - tape takes forever
# Not using it properly
# Run once to create tape. Then a seperate function to evaluate it.

test_tape <- optim(log(c(0.5,100000)),fn=ll_objfun_tape,gr=ll_gradfun_tape,method="BFGS", catch=data$catch,cpue=data$cpue)
test <- optim(log(c(0.5,100000)),fn=ll_objfun,gr=ll_gradfun,method="BFGS", catch=data$catch,cpue=data$cpue)
# CG - shite
test <- optim(log(c(0.5,100000)),fn=ll_objfun_tape,gr=ll_gradfun_tape,method="CG", catch=data$catch,cpue=data$cpue)
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


