#********************************************************************************
# Orange roughy assessment using FLaspm
#********************************************************************************

# Preliminaries
library(FLCore)
setwd("~/Work/deepfishman/FLaspm")

source("~/Sandbox/flr/pkg/FLCore/R/FLAccesors.R")
source("code/R/FLaspm_multipleindices.r")
source("code/R/FLModel_overloads.r")

# Input data
or_dat <- read.csv("Orange_roughy_assessment/FLASPM/or_dat.csv")

# Make these into FLQuant objects
landings <- FLQuant(or_dat$landings, dimnames=list(year=or_dat$year),units="t")
index <- FLQuant(or_dat$cpue, dimnames=list(year=or_dat$year),units="kg/hr")

# DEFINE LH PARAMETERS
# ages for model
amin    <- 0
amax    <- 100
ages <- amin:amax
nag     <- length(ages)

# steepness of BH recruitment
# Why is this an FLQuant? Mutiple iters?
hh <- FLQuant(0.75)

# natural mortality vector
# Not age specific
M <- FLQuant(0.45)

# age-length
Linf  <- 60 #cm
k     <- 0.65
t0    <- 0

# length-weight
alpha <- 0.022 # grams
beta  <- 2.95

age_to_length <- function(age,Linf,k,t0)
  return(Linf * (1 - exp(-k * (age - t0))))

length_to_weight <- function(l,a,b)
  return(a * l^b)
  
age_to_weight <- function(age,Linf,k,t0,a,b)
  return(length_to_weight(age_to_length(age,Linf,k,t0),a,b))

#l <- age_to_weight(ages,Linf,k,t0,alpha,beta)
l <- age_to_length(ages,Linf,k,t0) # in cm
w <- FLQuant(length_to_weight(l,alpha,beta),dimnames=list(age=ages)) # in g
w <- w/1e6 # tonnes, same as catch

# age at maturity
am <- 35

# maturity vector
m <- FLQuant(0, dimnames=list(age=ages))
m[(am - amin + 1):nag,] <- 1

# age at commercial selectivity (age of recruitment to fishery)
as <- 32

# selectivity vector
s <- FLQuant(0, dimnames=list(age=ages))
s[(as - amin + 1):nag,] <- 1

#*******************************************************************************
# Create the FLaspm objects

# pure R
or <- FLaspm(catch=landings, index=FLQuants(index),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(or) <- aspm()

# C version
or_c <- FLaspm(catch=landings, index=FLQuants(index),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(or_c) <- aspm_ad()

# AD version
or_ad <- FLaspm(catch=landings, index=FLQuants(index),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(or_ad) <- aspm_ad()
or_ad@gr <- get_gradient

# Chop off years
# some start years horrible
# Easy forward search
# But profile still horrible
# Check equib is OK
# Doesn't look right atm, recruits is off, but b is OK (?)
# Check timing of fishing too

# initial guess for B0 and sigma2
B0 <- 4000#500#1000 #5600
#B0 <- 1000
sigma2 <- 0.01#800
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

dyn.load("Orange_roughy_assessment/FLASPM/flaspm_ad.so")
# Fit
or.res <- fmle(or,start=start, lower=lower,upper=upper, control=list(trace=1))
or_c.res <- fmle(or_c,start=start, lower=lower,upper=upper, control=list(trace=1))
or_ad.res <- fmle(or_ad,start=start, lower=lower,upper=upper, control=list(trace=1))
dyn.unload("Orange_roughy_assessment/FLASPM/flaspm_ad.so")

# Parameter estimates the same
params(or.res)
params(or_c.res)
params(or_ad.res)

# Hessian slightly different - useful for bootstrapping
or.res@hessian
or_ad.res@hessian

# PMOD
# B0 = 5500
# lambda = 27
# estq = 1.03e-4
# cv = 0.24
# index = q * biomass

#******************************************************************************
# 2D Profiling
B0range <- params(or_ad.res)["B0"] * seq(from=0.2, to=1.8, length=100)
sigma2range <- params(or_ad.res)["sigma2"]* seq(from=0.2, to=1.8, length=100)

param_grid <- expand.grid(B0=B0range,sigma2=sigma2range,logl=NA)
for (i in 1:nrow(param_grid))
    param_grid[i,"logl"] <- logl(or_ad.res)(param_grid[i,"B0"],param_grid[i,"sigma2"],or_ad.res@hh, or_ad.res@M, or_ad.res@mat, or_ad.res@sel, or_ad.res@wght, or_ad.res@amin, or_ad.res@amax, or_ad.res@catch, or_ad.res@index)

logl_grid <- array(param_grid[,"logl"],dim=c(length(B0range),length(sigma2range)),dimnames=list(B0=B0range,sigma2=sigma2range))

# x and y the right way round?
image(x=B0range,y=sigma2range,z=logl_grid)
contour(x=B0range,y=sigma2range,z=logl_grid,add=T)
# Add out fitted points
points(x=params(or_ad.res)["B0"],y=params(or_ad.res)["sigma2"])


persp(x=B0range,y=sigma2range,z=logl_grid)

#???? B0 = 400?

#******************************************************************************
# 1D profiling
# Fixed B0, estimate sigma2
# Not use AD at the moment (returns two gradients - not one)
#B0range <- params(or_ad.res)["B0"] * seq(from=0.2, to=1.8, length=100)
B0range <- seq(from=1,to=3000,length=50)
B0range <- seq(from=1,to=10000,length=200)
prof1D <- cbind(B0=B0range,sigma2=NA,logl=NA)
dyn.load("Orange_roughy_assessment/FLASPM/flaspm_ad.so")
for (i in 1:nrow(prof1D))
{
    out <- fmle(or_c,start=list(sigma2=1), lower=1e-9,upper=1e10, fixed=list(B0=prof1D[i,"B0"]),control=list(trace=0))
    prof1D[i,"sigma2"] <- params(out)["sigma2"]
    prof1D[i,"logl"] <- c(logLik(out))
}
dyn.unload("Orange_roughy_assessment/FLASPM/flaspm_ad.so")

plot(x=prof1D[,"B0"],y=prof1D[,"logl"],type="l",xlab="B0",ylab="Log likelihood")
# plot fitted max
lines(x=c(params(or_ad.res)["B0"],params(or_ad.res)["B0"]),y=c(-1e6,1e6),lty=2)
maxcrude <- prof1D[which.max(prof1D[,"logl"]),"B0"]
lines(x=c(maxcrude,maxcrude),y=c(-1e6,1e6),lty=2)
pmodB0 <- 5500
lines(x=c(pmodB0,pmodB0),y=c(-1e6,1e6),lty=3)

# Local minima...
# so if you start at B0 = 1000, wrong side of the minima
# Really weird...

#******************************************************************************
# Plots
# B = N * mat * wght i.e. ssb
# Bexp = N * sel * wght i.e. biomass available to the fishery

# Are we returning the correct B and Bexp?
# index / q, not indexhat / q?
# and for PMOD B0 is estimated for 1991 (first index year) not 1989 (first catch year)
# eh?
# Are they basing their B0 on the estimated value of q and index?
# Not the fitted B0?

# Set up object to be as official assessment
or_ass <- or_ad.res
params(or_ass)["B0"] <- 5625

# Best, B
# index_hat
dyn.load("Orange_roughy_assessment/FLASPM/flaspm_ad.so")
out_fit <- get_aspm(or_ad.res)
out_ass <- get_aspm(or_ass)
dyn.unload("Orange_roughy_assessment/FLASPM/flaspm_ad.so")

index_data <- rbind(
	cbind(as.data.frame(or_ad.res@index[[1]]),name="index"),
	cbind(as.data.frame(out_fit$indexhat[[1]]),name="indexhat_fit"),
	cbind(as.data.frame(out_ass$indexhat[[1]]),name="indexhat_ass"),
	cbind(as.data.frame(or_ad.res@catch),name="catch")
	)

xyplot(data ~ year, groups=name, data=index_data, type="l", auto.key=T)
# Looks really bad!
# indexhat is pretty far from index

yrs <- as.numeric(dimnames(or_ad.res@catch)$year)

plot(x=yrs,y=out_ass$B,xlab="Yrs",ylab="SSB",type="l",ylim=c(0,6000))
lines(x=yrs,y=out_fit$B, lty=2)

# According to assessment
B_ass <- FLQuant(c(5625, 2974, 1668, 1386, 1361, 1402, 1442, 1451),dimnames=list(year=1991:1998))


# Plot fit at pmodB0 and 400

#******************************************************************************
test <- aspm.pdyn(or.res@catch,or.res@index[[1]],5600,c(or.res@hh),c(or.res@M),c(or.res@mat),c(or.res@sel),c(or.res@wght),or.res@amin,or.res@amax)

# Project forward with 0 catch
# Is equib population right?
C0 <- or.res@catch
C0[] <- 0
test <- aspm.pdyn(C0,or.res@index[[1]],500,c(or.res@hh),c(or.res@M),c(or.res@mat),c(or.res@sel),c(or.res@wght),or.res@amin,or.res@amax)
# b = 500 for all years, rec = same for all years
# so equb pop at start is fine
# B0 sets b in model = mature biomass = n * mat * wght

# how is q being calced?
# index = q * Bexp
# q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))

# why is index hat so low at start?

# catch in 1991 = 3502
# bexp in 1991 = 3506
# but biomass seems OK in 1992
# big recruitment coming through?

# y = 3
# bexp[3] = 1917
# C[3] = 3502
# h = C / bexp = 1.8
# but then h corrected using 
#h[y] <- min(h[y],0.999)
#    bexp[y] <- C[y] / h[y]
# b is not corrected
# next n is calced using new h (which == 1)
# It is this correction that is the problem. C > bexp should kill the system
# Instead h is set to 1 and bexp set to C
# but next n is not hammered enough?
# b is not right? - based on current b? timing

#exp(-M) * (1 * sel * h) = 0.000637
#so why is n not 10000 x smaller?
# it is, but only for age 34+

# y1, total mass of 1 - 5 year olds is very high
# n*w*mat, peak mat at 36
# n*w*sel, peak at 33,
# so even though total B will be high, reproductive B will be low
# When C > B, h <- 0.999
# n at sel ages + <- 0
# reproduction will be close to 0

# y4, n is *1e-3, of y3 (for sel range)
# b in y4 is *1e-3 of y3
# so rec should be low in y5, it is
# but then recovers as C decreases
# So everything looks sensible.
# what about the 999, make it 999999999?
# Makes no difference. Then what is it...?
# B0 year.

dyn.load("Orange_roughy_assessment/FLASPM/flaspm_ad.so")
out_fit <- get_aspm(or_ad.res)
dyn.unload("Orange_roughy_assessment/FLASPM/flaspm_ad.so")

# In Francis 1992
# B0 is the virgin exploitable biomass
# Has slightly different set up

#*******************************************************************************
# Assessment has B0 in year 1991 despite first catch year being 1989
# Chop off first two years of catch and final year

landings91 <- window(landings,start=1991,end=1998)
index91 <- window(index,start=1991,end=1998)

# Create the FLaspm objects

# pure R
or <- FLaspm(catch=landings91, index=FLQuants(index91),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(or) <- aspm()

# C version
or_c <- FLaspm(catch=landings91, index=FLQuants(index91),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(or_c) <- aspm_ad()

# AD version
or_ad <- FLaspm(catch=landings91, index=FLQuants(index91),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(or_ad) <- aspm_ad()
or_ad@gr <- get_gradient


# initial guess for B0 and sigma2
#B0 <- 6000#500#1000 #5600
B0 <- 1000
sigma2 <- 0.01#800
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

dyn.load("Orange_roughy_assessment/FLASPM/flaspm_ad.so")
# Fit
or.res <- fmle(or,start=start, lower=lower,upper=upper, control=list(trace=1))
or_c.res <- fmle(or_c,start=start, lower=lower,upper=upper, control=list(trace=1))
or_ad.res <- fmle(or_ad,start=start, lower=lower,upper=upper, control=list(trace=1))
dyn.unload("Orange_roughy_assessment/FLASPM/flaspm_ad.so")

# Parameter estimates the same
params(or.res)
params(or_c.res)
params(or_ad.res)

# inital B0 = 500, final = 313 (mature) > catch in 91 3502
# initial B0 = 

#******************************************************************************

# Doesn't go anywhere
# Why?
# Check output of logl and biomass estimates? scaling of w?

or@logl(B0,sigma2,hh,M,m,s,w,amin,amax,or@catch,or@index)
test_index <- aspm.index(or@catch,or@index,B0,hh,M,m,s,w,amin,amax) # FLQS
aspm.pdyn(or@catch,or@index,B0,c(hh),c(M),c(m),c(s),w,amin,amax)

B0range <- 5600 * seq(from=0.5, to =1.5, length=20)
loglout <- rep(NA,length(B0range))
for (i in 1:length(B0range))
  loglout[i] <- or@logl(B0range[i],sigma2,hh,M,m,s,w,amin,amax,or@catch,or@index)
# not going anywhere. sigma2?

# plot landings and index
# any info?
# maybe surface is just so flat. need adolc?
plot(x=1989:1999,y=c(landings/mean(landings)),type="l")
lines(x=1989:1999,y=c(index)/mean(c(index),na.rm=T),lty=2)


# PMOD
# B0 = 5500
# lambda = 27
# estq = 1.03e-4
# cv = 0.24
# index = q * biomass

pmod <- aspm.pdyn(or@catch,or@index,5500,c(hh),c(M),c(m),c(s),w,amin,amax)




#*******************************************************************************
# RUN FLASPM

# initial guess
#B0 <- 5600
#sigma2 <- 800

#start <- list(B0 = B0, sigma2 = sigma2)
#lower <- rep(1e-9,2)
#upper <- rep(1e10,2)

# call fmle
#sra.res <- fmle(sra.object,start=start,lower=lower,upper=upper,seq.iter=FALSE)
#params(sra.res)

# repeat with different initial values to check convergence

# PLOT OUTPUTS

# spawning biomass
xl <- 'Year'
yl <- 'Biomass (tonnes)'
ml <- 'Spawning Biomass'
plot(biomass(sra.res,'B'),type='l',lwd=2,xlab=xl,ylab=yl,main=ml)

# exploitable biomass
xl <- 'Year'
yl <- 'Biomass (tonnes)'
ml <- 'Exploitable Biomass'
plot(biomass(sra.res,'Bexp'),type='l',lwd=2,xlab=xl,ylab=yl,main=ml)

# harvest rate
xl <- 'Year'
yl <- 'Harvest Rate (Catch/Bexp)'
ml <- 'Harvest Rate'
plot(harvest.rate(sra.res),type='l',lwd=2,xlab=xl,ylab=yl,main=ml)

# fit to index
plot.fit(sra.res)

# END
