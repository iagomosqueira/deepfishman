
library(FLCore)
source("../../../code/adolc_tapeless/FLaspm_adolc.r")
dyn.load("../../../code/adolc_tapeless/fit_adolc_tapeless.dll")

# INPUT DATA

# catch (tonnes)
cc <- scan('landings.dat')
cc.units <- 'tonnes'
ymin   <- 1989
ymax   <- 1999
yr     <- as.character(ymin:ymax)
sra.catch <- FLQuant(as.vector(cc),dimnames=list(year=yr),units=cc.units)

# abundance index
cpue <- scan('cpue.dat')
cpue.units <- 'kg/hour'
ymin <- 1991
ymax <- 1998
yr <- as.character(ymin:ymax)
sra.index <- FLQuant(as.vector(cpue),dimnames=list(year=yr),units=cpue.units)

# align dimensions
sra.index <- mcf(FLQuants(index=sra.index,catch=sra.catch))[['index']]
units(sra.index) <- cpue.units

# DEFINE LH PARAMETERS

# ages for model
amin    <- 0
amax    <- 100
age.rng <- amin:amax
nag     <- length(age.rng)

# steepness
hh <- 0.75

# natural mortality vector
M <- rep(0.045,nag)

# age-length
Linf  <- 60 #cm
k     <- 0.65
t0    <- 0

# length-weight
alpha <- 0.022 # grams
beta  <- 2.95

# age-weight
l <- vector('numeric',nag)
w <- vector('numeric',nag)
for(a in 1:nag) {
  age <- age.rng[a]
  l[a] <- Linf * (1 - exp(-k * (age - t0)))
  w[a] <- alpha * l[a]^beta
}

w <- w/1e6 # tonnes

# age at maturity
am <- 35

# maturity vector
m <- vector("numeric",length=nag)
names(m) <- age.rng
m[1:(am - amin)] <- 0
m[(am - amin + 1):nag] <- 1

# age at commercial selectivity
as <- 32

# selectivity vector
s <- vector("numeric",length=nag)
names(s) <- age.rng
s[1:(as - amin)] <- 0
s[(as - amin + 1):nag] <- 1


# PRELIMINARIES

# create the FLaspm object
sra.object <- FLaspm(catch=sra.catch,
  index=sra.index,
  M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)

# load functions
model(sra.object) <- aspm()

# load gradient function
sra.object@grad <- function(x) {
    
  B0 <- x[1]
  sigma2 <- x[2]
  -pdyn(B0,sigma2,sra.object@catch,sra.object@index,sra.object@hh,sra.object@M,sra.object@mat,sra.object@sel,sra.object@wght,sra.object@amin,sra.object@amax)[['grad']]
}

# RUN FLASPM

# initial guess
B0 <- 5600
sigma2 <- 800

start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

# call fmle
sra.res <- fmle(sra.object,start=start,lower=lower,upper=upper,seq.iter=FALSE)
params(sra.res)

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