# Example with multiple indices
setwd("m:/Projects/Deepfishman/Orange_roughy_assessment/FLASPM")
#setwd("~/Work/deepfishman/FLaspm/examples/R/alb")
library(FLCore)
#load("IO_ctot.RData")

# Need to source FLAccesors from FLCore to allow us to write method for fmle and predict
# If we had a package we wouldn't need to do this
source("m:/Sandbox/flr/pkg/FLCore/R/FLAccesors.R")
source("m:/Projects/Deepfishman/deepfishman/FLaspm/code/R/FLaspm_multipleindices.r")
source("m:/Projects/Deepfishman/deepfishman/FLaspm/code/R/FLmodel_overloads.r")

# INPUT DATA
or_dat <- read.csv("or_dat.csv")

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
w <- length_to_weight(l,alpha,beta) # in g
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

# Create the FLaspm object
or <- FLaspm(catch=landings, index=FLQuants(index),
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)

model(or) <- aspm()

# initial guess for B0 and sigma2

B0 <- 1000 #5600
sigma2 <- 1#800
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)
# Fit
or.res <- fmle(or,start=start, lower=lower,upper=upper, control=list(trace=1))
#params(or.res)

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