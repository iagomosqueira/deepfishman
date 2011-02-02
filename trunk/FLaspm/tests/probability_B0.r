# Estimation of a probability distribution of B0
# using Monte Carlo simulation
# "If this were the true value of B0, how likely is it that we would get
# an estimate of virgin biomass as large as or greater than B0?"

# Procedure
# 1) Choose a trial value of B0
# 2) Calculate biomass history for that B0
# 3) Generate simulated indices where we assume indices are normally distributed with
# mean value equal to the biomass history (i.e. q = 1), and CV (c) = 15%.
# 4) Calculate max likelihood B0sim for the simulated survey indices
# (i.e. estimate B0sim for the simulated indices)
# 5) Repeat steps 3 and 4 100 times and calc the proportion of B0sim that are
# greater than 'true' B0.
# 6) Repeat 1 - 5 for a range of trial B0 values

library(FLaspm)
data(NZOR)
catch <- NZOR[["catch"]]
index <- NZOR[["index"]]

# Parameter values
amin    <- 1
amax    <- 70
hh <- 0.95
M <- 0.05
mat_age <- 23
sel_age <- 23
Linf  <- 42.5 #cm
k     <- 0.059
t0    <- -0.346
alpha <- 0.0963 # grams
beta  <- 2.68
mass_at_age <- age_to_weight(amin:amax,Linf,k,t0,alpha,beta)
w <- FLQuant(mass_at_age, dimnames=list(age=amin:amax)) / 1e6

index.cv <- 0.15
niters <- 500

# Deterministic projection
or <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or) <- aspm.Francis.C()
or <- fmle(or)

B0det <- params(or)['B0',]
profile(or)

# Set up simulated object
orsim <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(orsim) <- aspm.Francis.C()
#model(orsim) <- aspm.Francis()
# Choose a trial B0
B0trial <- 600000
B0trial <- 500000
B0trial <- 450000
B0trial <- 400000
B0trial <- 350000 # bmid = 0? or 1e-9? 0.
# what the hell is going on with h? it's like 2e31?
B0trial <- 300000 # qhat limits at 1e10, good. But logl still goes to NaN. Why?
# Because chat goes to inf

params(orsim)['B0',] <- B0trial
# Get the population trajectory
orsim.traj <- calc.pop.dyn(orsim)[["bexp"]]

# cv = sd(x) / mean (x)
# Generate spoof index - do 100 iterations
#index.sim <- rnorm(dim(orsim.traj)[2],orsim.traj,index.cv*orsim.traj)
set.seed(0)
index.sim <- rnorm(niters,orsim.traj,index.cv*orsim.traj)
#index.2sim <- iter(index.sim,2:100)

# check this is right
#yr <- 6
#mean(c(index.sim[,yr]))
#sd(c(index.sim[,yr])) / mean(c(index.sim[,yr]))

# ormonte - fit with the spoof index - multiple iterations
or.monte <- FLaspm(catch=catch, index=index.sim/100, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, amax=amax, amin=amin)
model(or.monte) <- aspm.Francis.C()
#model(or.monte) <- aspm.Francis() # works with R version
or.monte <- fmle(or.monte)
params(or.monte)

# now try vs R version - same? for 350000
# That fmax thing might have been it. No need for INFINITY check at all

# Do the methods work for multiple iters?
bexp <- exp.biomass(or.monte)
bmat <- mat.biomass(or.monte)
h <- harvest(or.monte)
n <- n(or.monte)
qhat <- calc.qhat(or.monte)
ll <- calc.logl(or.monte)
# Don't know why there is a warning - results look OK
ih <- indexhat(or.monte)
sweep(bexp,c(1,3:6),qhat[[1]],"*")

# Proportion of B0 estimates from or.monte that are greater than B0hat
# If B0trial was the 'true' value of B0, how likely is it that the estimated
# B0 (from or.monte, where q ~ 1), is greater than B0hat (from the deterministic fit)
# B0hat uses an index that has a CV. i.e. the index has a random component
length(which(c(params(or.monte)) > c(B0det))) / niters
# B0trial = 600000, p = 1
# B0trial = 500000, p = 1
# B0trial = 450000, p = 0.98
# B0trial = 400000, p = 0.5
# B0trial = 350000, p = 0


testiter <- 48 # 350000
testor <- FLaspm(catch=catch, index=iter(index.sim,testiter)/100, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(testor) <- aspm.Francis.C()
#model(testor) <- aspm.Francis()
testor <- fmle(testor)
profile(testor, maxsteps=30,range=0.1)
params(testor)

# R version works... eh...

# Set B0 to be first value of index
params(testor)["B0",] <- iter(index.sim[,1],testiter)
profile(testor,range=0.3)







B0 <- 448682.9
testor@logl(B0, testor@hh, testor@M, testor@mat, testor@sel, testor@wght,
            testor@amin, testor@amax, testor@catch, testor@index)
# -125.7144
or.monte@logl(B0, or.monte@hh, or.monte@M, or.monte@mat, or.monte@sel, or.monte@wght,
            or.monte@amin, or.monte@amax, or.monte@catch, iter(or.monte@index,2))
# -125.7144 too
# This is not what I get inside fmle for second iter though...
# I get -125.3268
# close - but why different?

params(or.monte)

c(params(or.monte)) > c(B0true)
# Proportion of B0monte that are greater than B0true

#*******************************************************************************
set.seed(0)
B0trial <- seq(from=300000,to=500000,length=40)
prop <- rep(NA,length(B0trial))
B0iters <- array(NA,dim=c(length(B0trial),niters))
for (i in 1:length(B0trial))
{
  cat("i: ", i, "\n")
  params(orsim)['B0',] <- B0trial[i]
  orsim.traj <- calc.pop.dyn(orsim)[["bexp"]]
  index.sim <- rnorm(niters,orsim.traj,index.cv*orsim.traj)/100

  # ormonte - fit with the spoof index - multiple iterations
  or.monte <- FLaspm(catch=catch, index=index.sim, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
  model(or.monte) <- aspm.Francis.C()
  #model(or.monte) <- aspm.Francis()
  or.monte <- fmle(or.monte,control=list(trace=0))
  B0iters[i,] <-c(params(or.monte))
}

# Proportion of B0monte that are greater than B0true
prop_B0det <- apply(B0iters>=c(B0det),1,sum) / niters
plot(B0trial,prop_B0det,type="l")
# Looks a bit off still
# Ahhhh, but Francis plot is of Bmid, so everything is slightly shunted back.
# Do this!




# in C code qhat returns Inf or NaN - need to return some maximum like 1e10
# How to test for infinite

# i = 1
# Test R version
testiter = 51
orr <- FLaspm(catch=catch, index=iter(index.sim,testiter), M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(orr) <- aspm.Francis()
orr <- fmle(orr,control=list(trace=0))
params(orr)
profile(orr,maxsteps=30)

orc <- FLaspm(catch=catch, index=iter(index.sim,testiter)/10, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(orc) <- aspm.Francis.C()
orc <- fmle(orc,control=list(trace=0))

params(orc)["B0",] <- 336177.1
calc.logl(orc) # but it works here...
calc.qhat(orc) # huge... could be this getting too big

params(orr)["B0",] <- 336177.1
calc.logl(orr) 
calc.qhat(orr) # Inf...


params(orc)["B0",] <- params(orr)["B0",]
profile(orc,maxsteps=30)
calc.logl(orc)

# why?
params(testor)["B0",] <- 487202.2
profile(testor,maxsteps=40)
# horrible!
# With R, even worse. What the F? Singularity?
# Is it the weird lumpy LL the otherside of the peak? Try the R version with the smooth ll
# Or try more steps in the initial
# WIth R version fails differently
# i=1, iter 16