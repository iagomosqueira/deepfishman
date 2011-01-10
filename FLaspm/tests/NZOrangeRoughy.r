#****************************************************************************
# Testing the common FLaspm methods for both methods on the NZ Orange Roughy data
#****************************************************************************

# Still to come:
# Methods for other outputs etc
# Test multiple iters
# Test multiple indices

# C code
# ADOLC code

#****************************************************************************
# Preliminaries
# Should have all this in a package to be happy it works
library(FLCore)
library(FLaspm)
#setwd("~/Work/deepfishman/FLaspm")
#source("~/Sandbox/flr/pkg/FLCore/R/FLAccesors.R")
#source("code/R/FLaspm_multipleindices.r")
#source("code/R/Edwards_model.r")
#source("code/R/Francis_model.r")
#source("code/R/FLModel_overloads.r")

# Fix parscale
# and pop.dyn methods
# then C
# then Charlie methods

#****************************************************************************

# Data and parameter values
index <- FLQuant(c(NA,NA,NA,NA,NA,164835,149425,102975,80397,97108,66291,NA),dimnames=list(year=1978:1989))
catch <- FLQuant(c(15340,40430,36660,32354,20064,32263,38142,39098,39896,31478,42621,37228),dimnames=list(year=1978:1989))
# Setting up the test data
amin    <- 1
amax    <- 70
# steepness of BH recruitment
hh <- FLQuant(0.95)
# natural mortality ,not age specific
M <- FLQuant(0.05)
# age-length
Linf  <- 42.5 #cm
k     <- 0.059
t0    <- -0.346
# length-weight
alpha <- 0.0963 # grams
beta  <- 2.68
am <- 23
as <- 23

# Set up FLQuants for weights, maturity and selectivity
# Could make a constructor function for this
w <- FLQuant(age_to_weight(amin:amax,Linf,k,t0,alpha,beta), dimnames=list(age=amin:amax)) / 1e6 # convert to tonnes
# age at maturity
# maturity vector
m <- FLQuant(0, dimnames=list(age=amin:amax))
m[(am - amin + 1):length(amin:amax),] <- 1
# age at commercial selectivity (age of recruitment to fishery)
s <- FLQuant(0, dimnames=list(age=amin:amax))
s[(as - amin + 1):length(amin:amax),] <- 1

#****************************************************************************
# Create the FLaspm object
res <- new("FLaspm")
res <- FLModel(class='FLaspm')


ed <- FLaspm(catch=catch, index=FLQuants(index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
ed <- FLaspm()
ed <- FLaspm(catch=catch)
# Set the Francis model
model(ed) <- aspm.Edwards()

fr <- FLaspm(catch=catch, index=FLQuants(index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()

# Dummy parameters - not fitted
B0 <- 411000* exp(0.5*c(M))
sigma2 <- 0.1 # for Charlie

ed@params['B0',] <- B0
ed@params['sigma2',] <- sigma2
fr@params['B0',] <- B0


# Test these methods
# pop.dyn
# exp.biomass
# n
# mat.biomass
# harvest (f or h)
ed.pop.dyn <- pop.dyn(ed)
fr.pop.dyn <- pop.dyn(fr)
exp.biomass(ed)
exp.biomass(fr)
mat.biomass(ed)
mat.biomass(fr)
n(ed)
n(fr)
harvest(ed)
harvest(fr)

plot(ed)
plot(fr)


fr@logl(B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)

#********************************************************************************
# Fitting test

B0 <- 426210#400000#500000
B0 <- 2262100#400000#500000
B0 <- 456000
B0 <- 700000
B0 <- 852420
B0 <- 213105
sigma2 <- 0.1
frstart <- list(B0 = B0)
frlower <- 1
frupper <- 1e10

edstart <- list(B0 = B0, sigma2=sigma2)
edlower <- c(1e-9,1e-9)
edupper <- c(1e10,1e10)

# Fit using default values
fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=TRUE)
fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=FALSE)
#fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1,parscale=10000))
fr.res@params


ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1),autoParscale=TRUE)
ed.res@params

# Hurrah

# Use default initial and bounds
fr.res <- fmle(fr,control=list(trace=1),autoParscale=TRUE)
fr.res@params

#********************************************************************************
# What happens if B0 too low, i.e. impossible catches
# Francis
fr <- FLaspm(catch=catch, index=FLQuants(index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()
B0 <- 20000
test <- aspm.pdyn.Francis(fr@catch,B0,c(fr@hh),c(fr@M),c(fr@mat),c(fr@sel),c(fr@wght),fr@amin,fr@amax) 

fr@logl(B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)
# B0 = 20000, f by year 7 has maxed out at 100, i.e. impossible catch
# bexp = 37576, C = 38142
# f = 100, return warning at end, B0 too low
# Screws up likelihood
# could put in test, if f=100 (or set maxf var) then logl = inf

# problem in pdyn and logl
# if bmid = 0, then index / bmid = inf
# set to 1? very hacky...
# how to fail gracefully...

#********************************************************************************
# Profile
# Method exists?
B0seq <- seq(from=200000, to=1000000,length=100)
frll <- rep(NA,length(B0seq))
for (i in 1:length(B0seq))
    frll[i] <- c(fr@logl(B0seq[i],fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index))
plot(B0seq,frll,type="l",xlab = "B0 start", ylab="logl")

#********************************************************************************
# Test creator

amin    <- 1
amax    <- 70
# steepness of BH recruitment
hh <- FLQuant(0.95)
# natural mortality ,not age specific
M <- FLQuant(0.05)
# age-length
Linf  <- 42.5 #cm
k     <- 0.059
t0    <- -0.346
# length-weight
a <- 0.0963 # grams
b  <- 2.68
matage <- 23
selage <- 23

# Pass nothing
test <- FLaspm()

# Pass the model only
test <- FLaspm(model=aspm.Francis())
model(test)

# Everything but model
test <- FLaspm(amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh,M=M)
test <- FLaspm(amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh)
test@M

# Everything
test <- FLaspm(model=aspm.Francis(),amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh,M=M)

# Everything and catch and indices
fr <- FLaspm(model=aspm.Francis(),amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh,M=M, catch=catch,index=index)
#test@fpm <- 1
res <- fmle(fr)
res@params
profile(res)
profile(res,maxsteps=30)
profile(res,maxsteps=30,range=0.3)


plot(res)
# Ace
