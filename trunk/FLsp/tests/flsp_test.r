# FLsp test

library(FLsp)

# Testing the C code
catch <- c(100,200,300)
r <- 1
p <- 1
k <- 1000
index <- c(90,100,60)
#test <- .Call("testflspCpp",catch)
test <- .Call("flspCpp",catch,index,r,p,k)

exp(sum(log(index / test[["B"]]),na.rm=T) / length(which(!is.na(index))))
test[["qhat"]]

test[["B"]] * test[["qhat"]];
test[["Ihat"]]

sum((log(index / test[["Ihat"]]))^2,na.rm=T) / length(which(!is.na(index)))
test[["sigma2"]]

# vhat log(index / test[["Ihat"]])
vhat <- log(index / test[["Ihat"]])
ll <- log(prod(exp(-(vhat^2) / (2*test[["sigma2"]])) / sqrt(2*pi*test[["sigma2"]])))

# Test gradients - look good
tiny <- 1e-9
rtiny <- r + tiny
testtiny <- .Call("flspCpp",catch,index,rtiny,p,k)
(testtiny[["ll"]] - test[["ll"]]) / tiny
test[["ll_grad_r"]]

ktiny <- k + tiny
testtiny <- .Call("flspCpp",catch,index,r,p,ktiny)
(testtiny[["ll"]] - test[["ll"]]) / tiny
test[["ll_grad_k"]]

# With missing indices
index <- c(NA,800,700)
test <- .Call("flspCpp",catch,index,r,p,k)
exp(sum(log(index / test[["B"]]),na.rm=T) / length(which(!is.na(index))))
test[["qhat"]]

detach("package:FLsp")
dyn.unload("C:/R/library/FLsp/libs/i386/FLsp.dll")







