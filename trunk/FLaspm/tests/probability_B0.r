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

# Deterministic projection
or <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or) <- aspm.Francis.C()
or <- fmle(or)

B0hat <- params(or)['B0',]

# CV = 15%
index.cv <- 0.15
niters <- 100

# Set up simulated object
orsim <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(orsim) <- aspm.Francis.C()
# Choose a trial B0
B0trial <- 500000
params(orsim)['B0',] <- B0trial
# Get the population trajectory
orsim.traj <- calc.pop.dyn(orsim)[["bexp"]]

# Figure out CV - cv = sd(x) / mean (x)
# Generate spoof index - do 100 iterations
#index.sim <- rnorm(dim(orsim.traj)[2],orsim.traj,index.cv*orsim.traj)
index.sim <- rnorm(niters,orsim.traj,index.cv*orsim.traj)
# check this is right
#yr <- 6
#mean(c(index.sim[,yr]))
#orsim.traj[,yr]
#sd(c(index.sim[,yr])) / mean(c(index.sim[,yr]))

index81 <- iter(index.sim,81)
or81 <- FLaspm(catch=catch, index=index81, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or81) <- aspm.Francis.C()
or81 <- fmle(or81)

# ormonte - fit with the spoof index - multiple iterations
or.monte <- FLaspm(catch=catch, index=index.sim, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or.monte) <- aspm.Francis.C()
# step 4 - fitting fails
# Problem with start
# only calls intitial once - not every time
or.monte <- fmle(or.monte)
params(or.monte)

or.monte@index # 100 iters
dims(or.monte)$iter # only picked up 1 - check with Iago how this works
# Might only just look at FLQuant slots, not FLQuants
# might need to overload dims for FLaspm
# Not just for picking up iters in fmle(), but in predict() too
# does iter work for FLQuants? Think so - browser through fmle to check


