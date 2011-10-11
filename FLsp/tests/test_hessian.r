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



# returns the gradients
simple_diff <- function(r,k)
{
	llorig <- .Call("flspCpp_tape",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}

simple_diff(r,k)

# Approximate hessian of non logged
tiny <- 1e-9
dorig <- simple_diff(r,k) # original gradients
drtiny <- simple_diff(r+tiny,k)
dktiny <- simple_diff(r,k+tiny)
dll2dr2 <- (drtiny[1]-dorig[1])/tiny
dll2dk2 <- (dktiny[2]-dorig[2])/tiny
dll2drk <-(dktiny[1]-dorig[1])/tiny
dll2dkr <-(drtiny[2]-dorig[2])/tiny

# Put all together
hesshat <- matrix(c(dll2dr2,dll2drk,dll2dkr,dll2dk2),nrow=2)

tape_all$hessian
hesshat

# Check gradient of logged function

# Hessian of logged
simple_diff_log <- function(r,k)
{
	llorig <- .Call("flspCpp_tape_log",data$catch,data$cpue,r,p,k)
	return(c(llorig[["ll_grad_r"]],llorig[["ll_grad_k"]]))
}

tiny <- 1e-6
dorig <- simple_diff_log(r,k) # original gradients
drtiny <- simple_diff_log(r+tiny,k)
dktiny <- simple_diff_log(r,k+tiny)
dll2dr2 <- (drtiny[1]-dorig[1])/tiny
dll2dk2 <- (dktiny[2]-dorig[2])/tiny
dll2drk <-(dktiny[1]-dorig[1])/tiny
dll2dkr <-(drtiny[2]-dorig[2])/tiny
hesshat <- matrix(c(dll2dr2,dll2drk,dll2dkr,dll2dk2),nrow=2)

tape_all_log$hessian
hesshat
