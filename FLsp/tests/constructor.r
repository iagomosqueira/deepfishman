# FLsp constructor test
setwd("c:/Projects/Deepfishman/deepfishman")
dat <- read.csv("trunk/FLsp/data/NZRL.csv")
library(FLsp)

catch <- FLQuant(dat$catch,dimnames=list(year=dat$year))
index <- FLQuant(dat$cpue,dimnames=list(year=dat$year))

test <- FLsp() # should return an empty object?
# Write the model function - returns ll etc
test <- FLsp(model=sp) # Fails as no catch or index (needed for validity check)


# This doesn't fail but should set model as sp by default
test <- FLsp(catch=catch, index=FLQuants(i1 = index))

test <- FLsp(catch=catch,index=FLQuants(i1=index),model=sp)

test@logl


