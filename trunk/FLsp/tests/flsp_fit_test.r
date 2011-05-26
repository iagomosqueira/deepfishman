# FLsp test

setwd("m:/Projects/Deepfishman/deepfishman")

library(FLsp)
dyn.load("C:/R/library/FLsp/libs/i386/FLsp.dll")

# Try SAA

saa <- read.csv("FLsp/data/NZRL.csv")
# Polacheck's results
p <- 1
r <- 0.0659
k <- 129000
q <- 2.461e-5
sigma <- 0.207
# what do we get?
saasp <- .Call("flspCpp",saa$catch,saa$cpue,r,p,k)
saasp[["qhat"]]
# Not bad!

ll_objfun <- function(params,catch,cpue)
{
  r <- params[1]
  k <- params[2]
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  return(-res[["ll"]])
}

ll_gradfun <- function(params,catch,cpue)
{
  r <- params[1]
  k <- params[2]
  res <- .Call("flspCpp",catch,cpue,r,1,k)
  return(c(-res[["ll_grad_r"]],-res[["ll_grad_k"]]))
}

test <- optim(c(1,1000),fn=ll_objfun,gr=ll_gradfun,method="L-BFGS-B", lower=c(1e-6,1e-6), upper=c(Inf,Inf),catch=saa$catch,cpue=saa$cpue)

test <- optim(c(1,1000),fn=ll_objfun,gr=ll_gradfun,method="BFGS", catch=saa$catch,cpue=saa$cpue)



detach("package:FLsp")
dyn.unload("C:/R/library/FLsp/libs/i386/FLsp.dll")







