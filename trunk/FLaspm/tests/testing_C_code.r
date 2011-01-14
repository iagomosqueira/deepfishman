#****************************************************************************
# Testing the common FLaspm methods for both methods on the NZ Orange Roughy data
# Checking that C code gives the same pop.dyn results
#****************************************************************************

#****************************************************************************
# Preliminaries
library(FLCore)
library(FLaspm)

#****************************************************************************

# Data and parameter values
data(NZOR)
catch <- NZOR[["catch"]]
index <- NZOR[["index"]]
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
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(ed) <- aspm.Edwards()

fr <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()

#****************************************************************************
# Test Edwards - R vs C
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

# Test Ed
edC <- .Call("aspm_ad", ed@catch, ed@index, B0, sigma2,
                ed@hh, ed@M, ed@mat, ed@sel,
                ed@wght, ed@amin, ed@amax, dim(ed@catch)[2], 1)

# Bexp
ed.pop.dyn[["bexp"]]
edC[["bexp"]]
all.equal(c(ed.pop.dyn[["bexp"]]),edC[["bexp"]])

#Bmat
ed.pop.dyn[["bmat"]]
edC[["bmat"]]
all.equal(c(ed.pop.dyn[["bmat"]]),edC[["bmat"]])

# h (or f)
ed.pop.dyn[["harvest"]]
edC[["harvest"]] # ?
all.equal(c(ed.pop.dyn[["harvest"]]),edC[["harvest"]])

# n
ed.pop.dyn[["n"]]
edC[["n"]]
all.equal(c(ed.pop.dyn[["n"]]),c(edC[["n"]]))

# qhat
# calced internally for R code atm
exp(mean(log(ed@index[[1]]/as.vector(ed.pop.dyn[["bexp"]])),na.rm=T))
edC[["qhat"]]

# index_hat
ed.pop.dyn[["bexp"]] * exp(mean(log(ed@index[[1]]/as.vector(ed.pop.dyn[["bexp"]])),na.rm=T))
edC[["indexhat"]]

# logl
ed@logl(B0,sigma2, ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
edC[["logl"]]

#****************************************************************************
# Test Francis - R vs C
B0 <- 411000* exp(0.5*c(M))
fr@params['B0',] <- B0


# Test these methods
# pop.dyn
# exp.biomass
# n
# mat.biomass
# harvest (f or h)
fr.pop.dyn <- pop.dyn(fr)

# Test Ed
frC <- .Call("aspm_ad", fr@catch, fr@index, B0, 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)

# Bexp
fr.pop.dyn[["bexp"]]
frC[["bexp"]]
all.equal(c(fr.pop.dyn[["bexp"]]),frC[["bexp"]])

#Bmat
fr.pop.dyn[["bmat"]]
frC[["bmat"]]
all.equal(c(fr.pop.dyn[["bmat"]]),frC[["bmat"]])

# h (or f)
fr.pop.dyn[["harvest"]]
frC[["harvest"]] # ?
all.equal(c(fr.pop.dyn[["harvest"]]),frC[["harvest"]])

# !!!! MODEL HAS CHANGED IN R !!! TIMING OF RECRUITMENT ETC.
# n
fr.pop.dyn[["n"]]
frC[["n"]]
all.equal(c(fr.pop.dyn[["n"]]),c(frC[["n"]]))


# qhat
# calced internally for R code atm
bmid <- fr.pop.dyn[["bexp"]] * exp(-0.5*(c(fr@M)+fr.pop.dyn[["harvest"]]))
qhat <- apply(fr@index[[1]]/bmid,c(1,6),sum,na.rm=T) / sum(!is.na(fr@index[[1]]))
qhat
frC[["qhat"]]

# indexhat
# Not calced for logl but used for residuals I think
sweep(bmid,1,qhat,"*")
frC[["indexhat"]]

#logl
fr@logl(B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)
frC[["logl"]]

#********************************************************************************
# Stop
#********************************************************************************
# Fitting test
# If B0 is set too low then are serious issues.
# The likelihood for a B0 lower than B0_min_for_survival is unclear. At the moment
# it is completely flat. This causes problems for the solver.
# So initial B0 has to be above this minimum
# Could try setting some minimum B0 but I don't know of a way to quickly get B0_min_for_survival

# Using default initial and bounds
fr.res <- fmle(fr)
fr.res@params

ed.res <- fmle(ed)
ed.res@params


# Or set the initial and boundary values by hand
# Try setting B0 has 50 x max catch
B0 <- 100*max(fr@catch)
frstart <- list(B0 = B0)
frlower <- 1e-9
frupper <- 1e10

fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=FALSE)
fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=TRUE)
fr.res@params

B0 <- 100*max(ed@catch)
sigma2 <- 1
edstart <- list(B0 = B0, sigma2=sigma2)
edlower <- c(1,1e-9)
edupper <- c(Inf,Inf)

# Same as default
ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1))
# This should be the same
ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1,parscale=c(2e8,0.4)))
# just shite
ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1,parscale=c(1e5,1)))
ed.res@params

#********************************************************************************
# Post fitting
# Profile
# Method exists?

# Is the fit any good?
profile(fr.res,maxsteps=30,range=0.3)
profile(ed.res,maxsteps=30,range=0.3)
profile(ed.res,maxsteps=30,range=0.2)
# ed looks a bit dodgy - take a closer look

# fitted logl
ed@logl(params(ed.res)["B0"],params(ed.res)["sigma2"],ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax,ed@catch,ed@index)
# better logl? Yes. Problem with fit
ed@logl(400000,0.024,ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax,ed@catch,ed@index)

#********************************************************************************
# What happens if B0 too low, i.e. impossible catches
# Francis
fr <- FLaspm(catch=catch, index=FLQuants(index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()
B0 <- 10#300000#20000
fr@params["B0",] <- B0
#test <- aspm.pdyn.Francis(fr@catch,B0,c(fr@hh),c(fr@M),c(fr@mat),c(fr@sel),c(fr@wght),fr@amin,fr@amax)
test <- pop.dyn(fr)
# harvest maxes out
# and abundance collapses

fr@logl(B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)

# problem in pdyn and logl
# if bmid = 0, then index / bmid = inf
# set to 1? very hacky...
# how to fail gracefully...

# Plot profile to reveal this
B0seq <- seq(from=10, to=800000,length=100)
frll <- rep(NA,length(B0seq))
for (i in 1:length(B0seq))
    frll[i] <- c(fr@logl(B0seq[i],fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index))
plot(B0seq,frll,type="l",xlab = "B0 start", ylab="logl")



#********************************************************************************

#.Call(aspm_ad, CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
#                SEXP steepnessSEXP, SEXP MSEXP, SEXP matSEXP, SEXP selSEXP,
#                SEXP wghtSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP nyrsSEXP


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
