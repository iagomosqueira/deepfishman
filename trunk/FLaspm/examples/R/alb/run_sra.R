
library(FLCore)
setwd("C:\\Projects\\Deepfishman\\deepfishman\\FLaspm\\examples\\R\\alb")
load("IO_ctot.RData")

source("../../../code/R/FLaspm.R")

# input required life history parameters

# age for model
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

# age at maturity
am <- 5

# maturity vector
m <- vector("numeric",length=nag)
names(m) <- amin:amax
m[1:(am - amin)] <- 0
m[am] <- 0.5
m[(am - amin + 2):nag] <- 1
#m[(am - amin + 1):nag] <- 1

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

#   create the FLaspm object
alb <- FLaspm(catch=alb.catch,
  index=alb.index,
  M=M,hh=hh,sel=s, mat=m, wght=mw, fpm=1, amax=amax, amin=amin)

model(alb) <- aspm()

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
