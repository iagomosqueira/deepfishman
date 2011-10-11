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

# Check log gradients
ll_r_log <- .Call("flspCpp_tape",data$catch,data$cpue,exp(log(r)+tiny),p,k)[["ll"]]
ll_k_log <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,exp(log(k)+tiny))[["ll"]]
dlldr_log <- (ll_r_log - ll_orig) / tiny
dlldk_log <- (ll_k_log - ll_orig) / tiny

dlldr_log
tape_all_log[["ll_grad_r"]]

dlldk_log
tape_all_log[["ll_grad_k"]]


# Approximate hessian of non logged
# returns the gradients
simple_diff <- function(r,k)
{
	llorig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}

dorig <- simple_diff(r,k) # original gradients
drtiny <- simple_diff(r+tiny,k)
dktiny <- simple_diff(r,k+tiny)
dll2dr2 <- (drtiny[1]-dorig[1])/tiny # change in grad_ll with change in r
dll2dk2 <- (dktiny[2]-dorig[2])/tiny
dll2drk <-(dktiny[1]-dorig[1])/tiny
dll2dkr <-(drtiny[2]-dorig[2])/tiny

# Put all together
hesshat <- matrix(c(dll2dr2,dll2drk,dll2dkr,dll2dk2),nrow=2)

tape_all$hessian
hesshat

# Hessian of logged
# We have checked that the tape_log function returns dll / d log(r)
simple_diff_log <- function(r,k)
{
	llorig <- .Call("flspCpp_tape_log",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}

tiny <- 1e-9
dorig <- simple_diff_log(r,k) # original gradients
drtiny <- simple_diff_log(exp(log(r)+tiny),k)
dktiny <- simple_diff_log(r,exp(log(k)+tiny))
dll2dr2 <- (drtiny[1]-dorig[1])/tiny
dll2dk2 <- (dktiny[2]-dorig[2])/tiny
dll2drk <-(dktiny[1]-dorig[1])/tiny
dll2dkr <-(drtiny[2]-dorig[2])/tiny
hesshat <- matrix(c(dll2dr2,dll2drk,dll2dkr,dll2dk2),nrow=2)

tape_all_log$hessian
hesshat

# looks fine
