# New constructor test
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
# original with mat and sel as FLQuants
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
summary(ed)

# mat as a single value
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=55, wght=w, fpm=1, amax=amax, amin=amin)
summary(ed)
ed@mat

# sel as single value
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=67, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
summary(ed)
ed@sel

# sel and mat as values
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
ed@sel
ed@mat

# hh and M as single values
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=0.2,hh=0.8,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
# good
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=c(0.2,0.3),hh=0.8,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
# fails - good
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=FLQuant(0.2),hh=0.8,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
# good

# index as a FLQuant
ed <- FLaspm(catch=catch, index=index, M=0.2,hh=0.8,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
ed@index
class(ed@index)

# validity checks if index and catch have same year range
ed <- FLaspm(catch=catch, index=index, M=0.2,hh=0.8,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
# good
badindex <- FLQuant(0,dimnames=list(year=1968:1979))
ed <- FLaspm(catch=catch, index=badindex, M=0.2,hh=0.8,sel=60, mat=50, wght=w, fpm=1, amax=amax, amin=amin)
# should fail but FLModel needs a bit of tinkering


