
library(FLCore)
library(FLaspm)
#data(IOalbacore)
load("../data/IOalbacore.RData")


# input required life history parameters

# age for model
amax <- 9
amin <- 1

# steepness
hh <- 0.75

# natural mortality
M <- 0.34

# mean weight (tonnes)
mw <- 0.014/1e3
mw <- FLQuant(mw,dimnames=list(age=amin:amax))    # temp

# age at commercial selectivity
as <- 5

# age at maturity
am <- 5

#   create the FLaspm object
alb <- FLaspm(catch=IOalbacore$catch,
              index=IOalbacore$index,
              M=M,hh=hh,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(alb) <- aspm.Edwards()

# initial guess for B0 and sigma2

B0 <- 120000/1e3
sigma2 <- 0.1

start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

alb.res <- fmle(alb,start=start, lower=lower,upper=upper,seq.iter=FALSE)
params(alb.res)

plot(1950:2007,as.vector(alb.res@fitted),ylim=c(0,5.5),xlab='Year',ylab='CPUE',type='l',col=2)
points(1950:2007,as.vector(alb.res@index))

# using logl() and optimise

obj.fn <- function(x) {
  return(-logl(alb)(x,sigma2,hh,M,m,s,mw,amin,amax,alb@catch,alb@index))
}

xx <- optimise(obj.fn,interval=c(0,500),maximum=FALSE)[['minimum']]
xx

# plot profile
B0.seq <- 80:160
lk <- c()
for(i in 1:length(B0.seq)) lk <- c(lk,obj.fn(B0.seq[i]))
plot(B0.seq,lk,type='l')

# accessor functions
exp.biomass(alb.res)
harvest.rate(alb.res)
