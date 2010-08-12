# Example with multiple indices
setwd("c:/Projects/Deepfishman/deepfishman/FLaspm/examples/R/alb")
library(FLCore)
load("IO_ctot.RData")

source("../../../code/R/FLaspm_multipleindices.r")

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

# split this into two indices that overlap
alb.index1 <- window(alb.index,start = 1950, end = 1998)
alb.index2 <- window(alb.index,start = 1995, end = 2007)
# But these need to have the same year dims
alb.indices <- mcf(FLQuants(index1 = alb.index1,index2 = alb.index2))
#alb.indices <- FLQuants(alb.index1,alb.index2)



# fmle needs to include a line similar to
#    datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))
# to include FLQuants objects
#getSlotNamesClass(alb, 'FLQuants')
# next problem is pulling out iters of a list
# e.g. data <- lapply(alldata, iter, it) if element all of alldata is itself a list
# might have to overload fmle for aspm
#dummy <- list(alb.indices)
# data <- lapply(dummy, iter, 1) # seems to work though...

#*******************************************************************************
# Single index as an FLQuants

alb <- FLaspm(catch=alb.catch,
  index=alb.indices,
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

# problem with predict at end of fmle
