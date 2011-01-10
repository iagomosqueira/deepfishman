# Examples with iterations
#setwd("c:/Projects/Deepfishman/deepfishman/FLaspm/examples/R/alb")
setwd("~/Work/deepfishman/FLaspm/examples/R/alb")
library(FLCore)
load("IO_ctot.RData")

source("../../../code/R/FLaspm_FS.R")

# input required life history parameters

# age for model
amax <- 9
amin <- 1
nag <- amax - amin + 1

# steepness
hh <- FLQuant(0.75)

# natural mortality
M <- FLQuant(0.34)

# mean weight
mw <- 0.014/1e3

# age at commercial selectivity
as <- 5

# selectivity vector
s <- vector("numeric",length=nag)
names(s) <- amin:amax
s[1:(as - amin)] <- 0
s[as] <- 0.5
s[(as - amin + 2):nag] <- 1
#s[(as - amin + 1):nag] <- 1
# make into FLQuant
sflq <- FLQuant(s,dimnames=list(age=amin:amax))

# age at maturity
am <- 5

# maturity vector
m <- vector("numeric",length=nag)
names(m) <- amin:amax
m[1:(am - amin)] <- 0
m[am] <- 0.5
m[(am - amin + 2):nag] <- 1
#m[(am - amin + 1):nag] <- 1
# make into FLQuant
mflq <- FLQuant(m,dimnames=list(age=amin:amax))


# input data
# catch
C <- subset(io_catch,Species=='ALB')[['Catch']]/1000
yr <- subset(io_catch,Species=='ALB')[['Year']]
alb.catch <- FLQuant(as.vector(C),dimnames=list(year=yr)) 

# abundance index
cpue_tw <- scan("cpue.dat"); ys <- 1980; yf <- 2006
names(cpue_tw) <- yr.i <- ys:yf
alb.index <- FLQuant(as.vector(cpue_tw),dimnames=list(year=ys:yf))
alb.index <- mcf(FLQuants(index=alb.index,catch=alb.catch))[['index']]  # align dimensions

#*******************************************************************************
# Single iteration fit and profile

alb <- FLaspm(catch=alb.catch,
  index=alb.index,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=mw, fpm=1, amax=amax, amin=amin)

model(alb) <- aspm()

# initial guess for B0 and sigma2

B0 <- 120000/1e3
sigma2 <- 0.1

start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

# Fit
alb.res <- fmle(alb,start=start, lower=lower,upper=upper)
params(alb.res)

plot(1950:2007,as.vector(alb.res@fitted),ylim=c(0,5.5),xlab='Year',ylab='CPUE',type='l',col=2)
points(1950:2007,as.vector(alb.res@index))

# Profile both B0 and sigma2
prof <- profile(alb.res)
prof <- profile(alb.res, maxsteps=20)

prof <- profile(alb.res, which="B0")




********************************************************************************
# Multiple iters
# Add iterations to the index as an example (for example, could bootstrap from the residuals)
niters <- 100
alb.index <- propagate(alb.index,iter=niters)
summary(alb.index)
# All iters currently the same
iter(alb.index,1)
iter(alb.index,2)
# multiply by lognormal noise to generate spoof index
alb.index <- alb.index * rlnorm(prod(dim(alb.index)),meanlog=0,sdlog=0.2)
# Now iters are different
iter(alb.index,1)
iter(alb.index,2)

#   create the FLaspm object
alb <- FLaspm(catch=alb.catch,
  index=alb.index,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=mw, fpm=1, amax=amax, amin=amin)

model(alb) <- aspm()

# check index slot has lots of iters

# initial guess for B0 and sigma2

B0 <- 120000/1e3
sigma2 <- 0.1

start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

#alb.res <- fmle(alb,start=start, lower=lower,upper=upper,seq.iter=FALSE)
# Fit all at once - might take a while
alb.res <- fmle(alb,start=start, lower=lower,upper=upper)
params(alb.res)

#plot(1950:2007,as.vector(alb.res@fitted),ylim=c(0,5.5),xlab='Year',ylab='CPUE',type='l',col=2)
#points(1950:2007,as.vector(alb.res@index))

# Have a look
bwplot(data ~ year, data = alb.res@fitted)

#*******************************************************************************
# Bootstrap

# Fit it once and look at residuals
amax <- 9
amin <- 1
nag <- amax - amin + 1
# steepness
hh <- 0.75
# natural mortality
M <- 0.34
# mean weight
mw <- 0.014/1e3
# age at commercial selectivity
as <- 5
# selectivity vector
s <- vector("numeric",length=nag)
names(s) <- amin:amax
s[1:(as - amin)] <- 0
s[as] <- 0.5
s[(as - amin + 2):nag] <- 1
#s[(as - amin + 1):nag] <- 1
# make into FLQuant
sflq <- FLQuant(s,dimnames=list(age=amin:amax))
# age at maturity
am <- 5
# maturity vector
m <- vector("numeric",length=nag)
names(m) <- amin:amax
m[1:(am - amin)] <- 0
m[am] <- 0.5
m[(am - amin + 2):nag] <- 1
#m[(am - amin + 1):nag] <- 1
# make into FLQuant
mflq <- FLQuant(m,dimnames=list(age=amin:amax))
# input data
# catch
C <- subset(io_catch,Species=='ALB')[['Catch']]/1000
yr <- subset(io_catch,Species=='ALB')[['Year']]
alb.catch <- FLQuant(as.vector(C),dimnames=list(year=yr)) 
# abundance index
cpue_tw <- scan("cpue.dat"); ys <- 1980; yf <- 2006
names(cpue_tw) <- yr.i <- ys:yf
alb.index <- FLQuant(as.vector(cpue_tw),dimnames=list(year=ys:yf))
alb.index <- mcf(FLQuants(index=alb.index,catch=alb.catch))[['index']]  # align dimensions

alb <- FLaspm(catch=alb.catch,
  index=alb.index,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=mw, fpm=1, amax=amax, amin=amin)

model(alb) <- aspm()

B0 <- 120000/1e3
sigma2 <- 0.1
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

# Fit
alb.res <- fmle(alb,start=start, lower=lower,upper=upper)

# Perform a bootstrap by resampling from residuals.
 
##sample with replacement by randomly selecting years
iters <- 100
yrs <- 1980:2006
mc.yrs <- sample(yrs,length(yrs)*iters,TRUE)

 
## then create an FLQuant for the residuals with the right dimensions
dmns     <-dimnames(alb.res@index)
dmns$iter<-1:100
dev.index<-FLQuant(NA,dimnames=dmns)
# pull out residuals from mc.yrs and 
dev.index[,ac(yrs)] <- alb.res@residuals[,ac(mc.yrs)]
# bootstrap the index 
alb.boot <- alb
alb.boot@index <- alb@index * exp(dev.index)
 

## rerun the assessment 100 times and put results in stock
alb.boot <- fmle(alb.boot,start=start, lower=lower,upper=upper)

# Plot some exciting things

#****************************************************************************
# Multiple iters on hh
niters <- 3
#hhflq <- FLQuant(hh,iter=niters)
hh <- propagate(hh,niters)
hh <- hh *rlnorm(niters,0,0.1)

C <- subset(io_catch,Species=='ALB')[['Catch']]/1000
yr <- subset(io_catch,Species=='ALB')[['Year']]
alb.catch <- FLQuant(as.vector(C),dimnames=list(year=yr))

# abundance index
cpue_tw <- scan("cpue.dat"); ys <- 1980; yf <- 2006
names(cpue_tw) <- yr.i <- ys:yf
alb.index <- FLQuant(as.vector(cpue_tw),dimnames=list(year=ys:yf))
alb.index <- mcf(FLQuants(index=alb.index,catch=alb.catch))[['index']]  # align dimensions


alb <- FLaspm(catch=alb.catch,
  index=alb.index,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=mw, fpm=1, amax=amax, amin=amin)

model(alb) <- aspm()

# check index slot has lots of iters

# initial guess for B0 and sigma2

B0 <- 120000/1e3
sigma2 <- 0.1

start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

#alb.res <- fmle(alb,start=start, lower=lower,upper=upper,seq.iter=FALSE)
# Fit all at once - might take a while
alb.res <- fmle(alb,start=start, lower=lower,upper=upper)
params(alb.res)


