# FLsp constructor test
#setwd("c:/Projects/Deepfishman/deepfishman")
#setwd("m:/Projects/Deepfishman/deepfishman")
setwd("~/Work/deepfishman")
dat <- read.csv("trunk/FLsp/data/NZRL.csv")
#dat <- read.csv("FLsp/data/NZRL.csv")
library(FLsp)

catch <- FLQuant(dat$catch,dimnames=list(year=dat$year))
index <- FLQuant(dat$cpue,dimnames=list(year=dat$year))

#test <- FLsp() # should return an empty object?
## Write the model function - returns ll etc
#test <- FLsp(model=sp) # Fails as no catch or index (needed for validity check)


# Sets model as sp by default
#test <- FLsp(catch=catch, index=FLQuants(i1 = index))

# Why does no longer work?
#test <- FLsp(catch=catch,index=FLQuants(i1=index),model=sp)

#test@logl

test <- FLsp(catch=catch, index=FLQuants(i1 = index))
test <- fitsp(test)
p <- predict(test)
evalC(test)
biomass(test)
qhat(test)
sigma2(test)
ll(test)

# Test fixed - broken
test <- fitsp(test, fixed=list(k = 1.4e5))

# Iter tests
# index with multiple iters
# catch with miters
# catch and index with multiple iters
# r and k with multiple iters
iters <- 10
indexit <- propagate(index,iters)
indexit <- indexit * rlnorm(prod(dim(indexit)),0,0.2)

catchit <- propagate(catch,iters)
catchit <- catchit * rlnorm(prod(dim(catchit)),0,0.2)

test1 <- FLsp(catch=catch, index=FLQuants(indexit))
test2 <- FLsp(catch=catchit, index=FLQuants(index))
test3 <- FLsp(catch=catchit, index=FLQuants(indexit))

test1 <- fitsp(test1)
test2 <- fitsp(test2)
test3 <- fitsp(test3)

hist(c(params(test1)["r"]))
hist(c(params(test2)["r"]))
hist(c(params(test3)["r"]))

# Problem is iter and index FLQuants
evalC(test1)
biomass(test1)


