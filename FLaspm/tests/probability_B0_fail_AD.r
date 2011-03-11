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

index.cv <- 0.15 # upper estimate 0.19
#index.cv <- 0.19
niters <- 250

or_true <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, amax=amax, amin=amin)
model(or_true) <- aspm.Francis.CAD()
or_true <- fmle(or_true,method="BFGS")
B0_true <- params(or_true)['B0',]

or <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, amax=amax, amin=amin)
model(or) <- aspm.Francis.CAD()

Bmidtrial <- seq(from=350000,to=600000,by=25000)
B0trial <- Bmidtrial*exp(0.5*M)
B0store <- array(NA,dim=c(length(B0trial),niters))

set.seed(0)
na.years <- dimnames(index)$year[which(is.na(index))]
## one more year
#na.years <- as.character(1978:1982)
## two more
#na.years <- as.character(1978:1981)
## 3 more
#na.years <- as.character(1978:1980)
## 4 more
#na.years <- as.character(1978:1979)
## 5 more
#na.years <- as.character(1978)

for (i in 1:length(B0trial))
{
  cat("i: ", i, "\n")
  params(or)['B0',] <- B0trial[i]
  or.traj.mid <- exp.biomass.mid(or,yrfrac=0.5,virgin=F)

  # Cannot fit with a trajectory that goes extinct so see if
  # F has maxed out over index years
  if(all(harvest(or)<100))
  {
    index.sim <- rnorm(niters,or.traj.mid,index.cv*or.traj.mid)
    # index.sim should have same years as original index
    index.sim[,na.years] <- NA
    #apply(index.sim,2,mean)
    #    cv(index.sim)

    or.sim <- FLaspm(catch=catch, index=index.sim, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, amax=amax, amin=amin)
    model(or.sim) <- aspm.Francis.CAD()
    #model(or.sim) <- aspm.Francis()
    #or.sim <- fmle(or.sim,always_eval_initial=TRUE,SANN_maxit=0,control=list(trace=0))
    or.sim <- fmle(or.sim,always_eval_initial=TRUE,method="BFGS",control=list(trace=0))
    B0store[i,] <-(params(or.sim))
  }
}



#save(B0store,B0trial,file="C:/Projects/Deepfishman/deepfishman/trunk/B0storeAD.Rdata")
#save(B0store,na.years,B0trial,file="C:/Projects/Deepfishman/deepfishman/trunk/B0storeAD_plus6.Rdata")
save(B0store,na.years,B0trial,file="~/Work/deepfishman/trunk/B0storeAD.Rdata")
#load("C:/Projects/Deepfishman/deepfishman/trunk/B0store.Rdata")

propB0 <- apply(B0store>=c(B0_true),1,sum) / niters
plot(B0trial,propB0,type="l")
lines(x=c(B0_true,B0_true),y=c(0,1),lty=2)

Bmid_true <- B0_true * exp(-0.5*M)
Bmidstore <- B0store * exp(-0.5*M)
propBmid <- apply(Bmidstore>=c(Bmid_true),1,sum) / niters
plot(Bmidtrial,propBmid,type="l")
lines(x=c(Bmid_true,Bmid_true),y=c(0,1),lty=2)

# Cut off NA
Bmidtrial <- Bmidtrial[!is.na(propB0)]
propB0 <- propB0[!is.na(propB0)]
propBmid <- propBmid[!is.na(propBmid)]
# Fit a Johnson's Su distribution to get the cumulative distribution function
# Try fitting Johnson distribution with real data
initialparms <- c(g=-1,delta=1,xi=0,lambda=1)
jcpars <- optim(par=initialparms,fn=Johnsonll,b=Bmidtrial/1000,m=niters,p=propBmid)$par
Bmidplot <- seq(from=300,to=600,length=100)
jp <- JohnsonPDF(jcpars,Bmidplot)
plot(Bmidplot,jp,type="l")
lines(x=c(Bmid_true,Bmid_true)/1000,y=c(0,1),lty=2)

# Reading from graph Francis has
# Bmid   ~p    me c = 0.15
# 350    0      0       
# 380    0.08   0.76      
# 400    0.35   0.35
# 450    0.8    0.77
# 500    0.92   0.92

rbind(Bmidtrial,propBmid)



#*******************************************************************
# Deterministic projection
or1 <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or1) <- aspm.Francis.C()
or1 <- fmle(or1)

# Index is the relative index of abundance
# Bmid ~ I, Bmid = QI
# So if set up the index to be same as the biomass, do I fit the same B0
for1 <- calc.pop.dyn(or)[["harvest"]]
bexpor1 <- calc.pop.dyn(or)[["bexp"]]
bmidor1 <- bexpor1 * exp(-0.5*(for1+c(or@M)))


or2 <- FLaspm(catch=catch, index=bmidor1, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or2) <- aspm.Francis.C()
or2 <- fmle(or2)
profile(or2,maxsteps=50)
# Add parameter of first fit
lines(x=c(params(or1)["B0",],params(or1)["B0",]),y=c(-1000,1000),lty=2)
params(or2)
params(or1)
# The same
calc.qhat(or2)
# 1
bexpor2 <- calc.pop.dyn(or2)[["bexp"]]
for2 <- calc.pop.dyn(or2)[["harvest"]]
bmidor2 <- bexpor2 * exp(-0.5*(for2+c(or2@M)))
bmidor2
bmidor1


# Add noise to index, taking the best estimate of bmid as the index
set.seed(0)
# normally distributed, cv 15%, sd / mean
niters <- 500
index.sim <- rnorm(niters,bmidor1,index.cv*bmidor1)
# check this is right
apply(index.sim,2,mean)
bmidor1
# check sd and cv is right
cv(index.sim)
apply(index.sim,2,function(x)sd(x)/mean(x))

or3 <- FLaspm(catch=catch, index=index.sim, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
model(or3) <- aspm.Francis.C()
or3 <- fmle(or3)
vbmid <- c(params(or3)) * exp(-0.5*c(or@M))
truehist(vbmid)
mean(c(params(or3)))
params(or1)
length(which(c(params(or3)) >= c(params(or1)))) / niters
# not 0.5
# it won't be

#*******************************************************************
#model(orsim) <- aspm.Francis()
# Choose a trial B0
B0trial <- 600000
B0trial <- 500000
B0trial <- 450000
B0trial <- 400000
B0trial <- 350000
B0trial <- 300000
B0trial <- 250000

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
