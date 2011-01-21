#****************************************************************************
# Testing the common FLaspm methods for both methods on the NZ Orange Roughy data
# Only looking at the R versions here
#****************************************************************************

# Still to come:
# Methods for other outputs etc
# Test multiple iters
# Test multiple indices
# Test Edwards model (compare to source paper)

# Make sure the pop.dyn sets the right units for harvest

#****************************************************************************
# Preliminaries
library(FLCore)
library(FLaspm)

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

# Set up FLQuant for weights
# Could make a constructor function for this
w <- FLQuant(age_to_weight(amin:amax,Linf,k,t0,alpha,beta), dimnames=list(age=amin:amax)) / 1e6 # convert to tonnes

#****************************************************************************
# Create the FLaspm objects
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=as, mat=am, wght=w, amax=amax, amin=amin)
model(ed) <- aspm.Edwards()

fr <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=as, mat=am, wght=w, amax=amax, amin=amin)
model(fr) <- aspm.Francis()

#****************************************************************************
# Accessor method check
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
ed.pop.dyn <- calc.pop.dyn(ed)
fr.pop.dyn <- calc.pop.dyn(fr)

exp.biomass(ed)
exp.biomass(fr)

mat.biomass(ed)
mat.biomass(fr)

n(ed)
n(fr)

harvest(ed)
harvest(fr)

#********************************************************************************
# Fitting test
# If B0 is set too low then are serious issues.
# The likelihood for a B0 lower than B0_min_for_survival is unclear. At the moment
# it is completely flat. This causes problems for the solver.
# So initial B0 has to be above this minimum
# Could try setting some minimum B0 but I don't know of a way to quickly get B0_min_for_survival

# Using default initial and bounds
# Takes a little time with the initial function but the results are much better
fr.res <- fmle(fr)
fr.res@params

ed.res <- fmle(ed)
ed.res@params


# Or set the initial and boundary values by hand
# Often gives quite bad results
#B0 <- 100*max(fr@catch)
#frstart <- list(B0 = B0)
#frlower <- 1e-9
#frupper <- 1e10

#fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=FALSE)
#fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=TRUE)
#fr.res@params

#B0 <- 100*max(ed@catch)
#sigma2 <- 1
#edstart <- list(B0 = B0, sigma2=sigma2)
#edlower <- c(1,1e-9)
#edupper <- c(Inf,Inf)

# Same as default
#ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1))
## This should be the same
#ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1,parscale=c(2e8,0.4)))
# just shite
#ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1,parscale=c(1e5,1)))
#ed.res@params

#********************************************************************************
# Post fitting
# Profile

# Is the fit any good?
# This is really slow for pure R version
profile(fr.res,maxsteps=30,range=0.3)
profile(ed.res,maxsteps=10,range=0.3)
profile(ed.res,maxsteps=30,range=0.3)
profile(ed.res,maxsteps=30,range=0.2)

#********************************************************************************
# Scratch below here

# What happens if B0 too low, i.e. impossible catches
# Francis
fr <- FLaspm(catch=catch, index=FLQuants(index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()
B0 <- 10#300000#20000
fr@params["B0",] <- B0
#test <- aspm.pdyn.Francis(fr@catch,B0,c(fr@hh),c(fr@M),c(fr@mat),c(fr@sel),c(fr@wght),fr@amin,fr@amax)
test <- calc.pop.dyn(fr)
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
