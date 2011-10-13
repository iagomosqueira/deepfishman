# Test Hessian and log Hessian
# Test using NZ Rock Lobster data

library(testthat)
library(FLsp)
data(nzrl)
data <- nzrl

p <- 1
# Polacheck's results
# NZRL
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207

#data$catch <- data$catch / 1000
#data$cpue <- data$cpue / 1000
#k <- 129

tape_all <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)
tape_all_log <- .Call("flspCpp_tape_log",data$catch,data$cpue,r,p,k)

# Check gradients
tiny <- 1e-9
ll_orig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)[["ll"]]
ll_r <- .Call("flspCpp_tape",data$catch,data$cpue,r+tiny,p,k)[["ll"]]
ll_k <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k+tiny)[["ll"]]
dlldr <- (ll_r - ll_orig) / tiny
dlldk <- (ll_k - ll_orig) / tiny

dlldr
tape_all[["ll_grad_r"]]

dlldk
tape_all[["ll_grad_k"]]


# Put this into a function
# If dll / d log(r) look at small change in log r, not r
grad_approx <- function(r,k,logrk=FALSE)
{
  tiny <- 1+1e-12
  ll_orig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)[["ll"]]
  ll_r <- .Call("flspCpp_tape",data$catch,data$cpue,r*tiny,p,k)[["ll"]]
  if (!logrk)
    dlldr <- (ll_r - ll_orig) / (r * tiny - r)
  if (logrk)
    dlldr <- (ll_r - ll_orig) / (log(r * tiny) - log(r))    

  ll_k <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k*tiny)[["ll"]]
  if (!logrk)
    dlldk <- (ll_k - ll_orig) / (k * tiny - k)
  if (logrk)
    dlldk <- (ll_k - ll_orig) / (log(k * tiny) - log(k))
  return(list(dlldr=dlldr,dlldk=dlldk))
}

grad_approx(r,k)
tape_all[["ll_grad_r"]]
tape_all[["ll_grad_k"]]

grad_approx(r,k,logrk=TRUE)
tape_all_log[["ll_grad_r"]]
tape_all_log[["ll_grad_k"]]

# Can we use this to get approximate hessian
# Change in gradient to change in r
tiny <- 1+1e-12
dll_orig <- grad_approx(r,k)
dll_r <- grad_approx(r*tiny,k)
dll_k <- grad_approx(r,k*tiny)
h11 <- (dll_r[[1]] - dll_orig[[1]]) / (r*tiny - r)
h22 <- (dll_k[[2]] - dll_orig[[2]]) / (k*tiny - k)
tape_all[["hessian"]]
# It's too approximate - just wrong... also method is just wrong see below
# So do it using true gradient
true_grad <- function(r,k)
{
	llorig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}

tiny <- 1+1e-9
dll_orig <- true_grad(r,k)
dll_r <- true_grad(r*tiny,k)
dll_k <- true_grad(r,k*tiny)

h11 <- (dll_r[1] - dll_orig[1]) / (r*tiny - r)
h12 <- (dll_k[1] - dll_orig[1]) / (k*tiny - k)
h21 <- (dll_r[2] - dll_orig[2]) / (r*tiny - r)
h22 <- (dll_k[2] - dll_orig[2]) / (k*tiny - k)
tape_all[["hessian"]]

# log hessian
h11_log <- (dll_r[1] - dll_orig[1]) / (log(r*tiny) - log(r))
tape_all_log[["hessian"]]
# Not right
# Hessian in log case is how much the log_gradient changes with log r,
# not how much the gradient changes with log r
true_grad_log <- function(r,k)
{
	llorig <- .Call("flspCpp_tape_log",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}
tiny <- 1+1e-9
dll_orig_log <- true_grad_log(r,k)
dll_r_log <- true_grad_log(r*tiny,k)
dll_k_log <- true_grad_log(r,k*tiny)

h11_log <- (dll_r_log[1] - dll_orig_log[1]) / (log(r*tiny) - log(r))
h12_log <- (dll_k_log[1] - dll_orig_log[1]) / (log(k*tiny) - log(k))
h21_log <- (dll_r_log[2] - dll_orig_log[2]) / (log(r*tiny) - log(r))
h22_log <- (dll_k_log[2] - dll_orig_log[2]) / (log(k*tiny) - log(k))
tape_all_log[["hessian"]]

#*******************************************************************************
# Question is, are we actually getting the right Hessian?

# We want the hessian from fitting log(r) and log(k)
# So if log(r) changes, how does ll change?
# The gradient is right
library(FLsp)
library(numDeriv)
data(nzrl)
data <- nzrl

p <- 1
# Polacheck's results
# NZRL
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207

k <- k / 1000
data$catch <- data$catch / 1000
#data$cpue <- data$cpue / 1000


ll <- function(x, catch, cpue, p)
{
  r <- x[1]
  k <- x[2]
  ll <- .Call("flspCpp_tape",catch,cpue,r,p,k)[["ll"]]
  if (is.nan(ll)) ll <- -1e10 # weight against NaNs
  return(ll)
}

ll(c(r/10,k), catch=data$catch, cpue=data$cpue, p=p)
numDeriv::hessian(func = ll, x =c(r,k), method.args=list(r=6),catch=data$catch, cpue=data$cpue, p=p)
tape_all <- Call("flspCpp_tape",data$catch,data$cpue,r,p,k)
tape_all[["hessian"]]
# Yup!

logll <- function(x, catch, cpue, p)
{
  
  r <- exp(x[1])
  k <- exp(x[2])
  ll <- .Call("flspCpp_tape",catch,cpue,r,p,k)[["ll"]]
  if (is.nan(ll)) ll <- -1e10 # weight against NaNs
  return(ll)
}

logll(log(c(r/10,k)), catch=data$catch, cpue=data$cpue, p=p)

numDeriv::hessian(func = logll, x =log(c(r,k)), method.args=list(r=6), catch=data$catch, cpue=data$cpue, p=p)
tape_all_log <- .Call("flspCpp_tape_log",data$catch,data$cpue,r,p,k)
tape_all_log[["hessian"]]
#Yup

# Our hessian is the hessian from fitting log(r) and log(k)


