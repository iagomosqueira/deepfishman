# Example with multiple indices
#setwd("m:/Projects/Deepfishman/deepfishman/FLaspm/examples/R/alb")
setwd("~/Work/deepfishman/FLaspm/examples/R/alb")
library(FLCore)
load("IO_ctot.RData")

# Need to source FLAccesors from FLCore to allow us to write method for fmle and predict
# If we had a package we wouldn't need to do this
#source("m:/Sandbox/flr/pkg/FLCore/R/FLAccesors.R")
source("~/Sandbox/flr/pkg/FLCore/R/FLAccesors.R")
source("../../../code/R/FLaspm_multipleindices.r")
source("../../../code/R/FLModel_overloads.r")
#source("~/Work/deepfishman/FLaspm/code/R/FLModel_overloads.r")


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
wght <- FLQuant(mw,dimnames=list(age=amin:amax))

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

alb <- FLaspm(catch=alb.catch,
  index=alb.indices,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

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

# seems to work

#******************************************************************************
#Multiple iters and indices
iters <- 5
hh_iters <- propagate(alb@hh,iters)
hh_iters[] <- c(alb@hh) * rlnorm(iters,0,0.2)

alb_iters <- FLaspm(catch=alb.catch,
  index=alb.indices,
  M=M,hh=hh_iters,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

model(alb_iters) <- aspm()

# initial guess for B0 and sigma2
B0 <- 120000/1e3
sigma2 <- 0.1
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)
# Fit
alb_iters <- fmle(alb_iters,start=start, lower=lower,upper=upper)

params(alb_iters)
alb_iters@fitted_index
alb_iters@residuals_index

# seems to work too.

#*******************************************************************************
# Using AD


aspm_ad <- function() {
  # logl
  logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
  {
      #browser()
      nyrs <- length(c(catch))
    out <- .Call("aspm",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,nyrs)
    total.logl <- out[["logl"]][["logl"]]
    #cat("total.logl", total.logl, "\n")
    #cat("out", out, "\n")
    return(total.logl)
  }
  
   # initial parameter values
   
  initial <- structure(function(catch) {
    return(FLPar(B0=100*max(catch), sigma2=1))
	},
  # lower and upper limits for optim()
	lower=c(1, 1e-8),
	upper=c(Inf, Inf)
	)

    #model <- index ~ aspm.index(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
# has to return FLQuants
#model <- index ~ .Call("aspm",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,length(c(catch)))[["indexhat"]]
    model <- index ~ get_indexhat(catch,index,B0,sigma2,hh,M,mat,sel,wght,amin,amax)
  
  return(list(logl=logl,model=model,initial=initial))
} # }}}

get_indexhat <- function(catch,index,B0,sigma2,hh,M,mat,sel,wght,amin,amax)
{
    #browser()
    out <-  .Call("aspm",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,length(c(catch)))[["indexhat"]]
    indexhat <- FLQuants()
    dnms <- dimnames(index[[1]])
    for (i in 1:length(index))
	indexhat[[i]] <- FLQuant(out[i,],dimnames=dnms)
    return(indexhat)
}

# Must have same args as logl
#logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
get_gradient <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
{
    out <-  .Call("aspm",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,length(c(catch)))[["logl"]]
    cat("in grad function ", unlist(out), "\n")
    #    return(unlist(out)[-1])
    return(-1*c(out[["logl_grad_B0"]],out[["logl_grad_sigma2"]]))
}


get_gradient2 <- function(B0,sigma2)
{
    out <-  .Call("aspm",(c(object@catch)),object@index,B0,sigma2,(c(object@hh)),(c(object@M)),(c(object@mat)),(c(object@sel)),(c(object@wght)),object@amin,object@amax,length(c(object@catch)))[["logl"]]
    #    cat("in grad function ", unlist(out), "\n")
    #    return(unlist(out)[-1])
    return(c(out[["logl_grad_B0"]],out[["logl_grad_sigma2"]]))
}

#***************************************************************************
# One index for the moment

# Original R
alb <- FLaspm(catch=alb.catch,
    index=FLQuants(alb.indices[[1]]),
    #index=alb.indices,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

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

# What is logl at start values of pure R
alb@logl(B0,sigma2,alb@hh,alb@M,alb@mat,alb@sel,alb@wght,alb@amin,alb@amax,alb@catch,alb@index)

# Set up the AD version
alb_ad <- FLaspm(catch=alb.catch,
    index=FLQuants(alb.indices[[1]]),
    #index=alb.indices,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

model(alb_ad) <- aspm_ad()

# What is logl at start values of pure R
alb@logl(B0,sigma2,alb@hh,alb@M,alb@mat,alb@sel,alb@wght,alb@amin,alb@amax,alb@catch,alb@index)

dyn.load("~/Work/deepfishman/FLaspm/Orange_roughy_assessment/FLASPM/flaspm_ad.so")
alb_ad@logl(B0,sigma2,alb@hh,alb@M,alb@mat,alb@sel,alb@wght,alb@amin,alb@amax,alb@catch,alb@index)
ih <- get_indexhat(alb_ad@catch,alb_ad@index,B0,sigma2,alb_ad@hh,alb_ad@M,alb_ad@mat,alb_ad@sel,alb_ad@wght,alb_ad@amin,alb_ad@amax,58)
dyn.unload("~/Work/deepfishman/FLaspm/Orange_roughy_assessment/FLASPM/flaspm_ad.so")
# For one index is fine





alb_ad <- FLaspm(catch=alb.catch,
	#    index=FLQuants(alb.indices[[1]]),
    index=alb.indices,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

model(alb_ad) <- aspm_ad()
# initial guess for B0 and sigma2
B0 <- 120000/1e3
sigma2 <- 0.1
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)
# Fit
so_dir <- "~/Work/deepfishman/FLaspm/Orange_roughy_assessment/FLASPM"
dyn.load(paste(so_dir,"flaspm_ad.so",sep="/"))
alb_ad_res <- fmle(alb_ad,start=start, lower=lower,upper=upper)
dyn.unload(paste(so_dir,"flaspm_ad.so",sep="/"))
#dyn.unload("flaspm_ad.so")
alb_ad_res@params

#********************************************************************
alb_ad <- FLaspm(catch=alb.catch,
	#    index=FLQuants(alb.indices[[1]]),
    index=alb.indices,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

model(alb_ad) <- aspm_ad()
# initial guess for B0 and sigma2
B0 <- 120000/1e3
sigma2 <- 0.1
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)
# Fit
so_dir <- "~/Work/deepfishman/FLaspm/Orange_roughy_assessment/FLASPM"
# Need to return index_hat
# As an FLQ?
dyn.load(paste(so_dir,"flaspm_ad.so",sep="/"))
out <- .Call("aspm",(c(alb_ad@catch)),alb_ad@index,B0,sigma2,(c(alb_ad@hh)),(c(alb_ad@M)),(c(alb_ad@mat)),(c(alb_ad@sel)),(c(alb_ad@wght)),amin,amax,58)
dyn.unload(paste(so_dir,"flaspm_ad.so",sep="/"))

# index hat
#**************************************************************************
# Can we return gradient too?
# yes but is shite!
alb_ad_gr <- FLaspm(catch=alb.catch,
	#index=FLQuants(alb.indices[[1]]),
    index=alb.indices,
  M=M,hh=hh,sel=sflq, mat=mflq, wght=wght, fpm=1, amax=amax, amin=amin)

model(alb_ad_gr) <- aspm_ad()
alb_ad_gr@gr <- get_gradient
alb_ad_gr@gr


dyn.load(paste(so_dir,"flaspm_ad.so",sep="/"))
alb_ad_res <- fmle(alb_ad,start=start, lower=lower,upper=upper)
alb_ad_gr <- fmle(alb_ad_gr,start=start, lower=lower,upper=upper)
dyn.unload(paste(so_dir,"flaspm_ad.so",sep="/"))
#dyn.unload("flaspm_ad.so")
alb_ad_gr@params
alb_ad_nogr@params

alb_ad_gr@hessian
alb_ad_nogr@hessian

# Whoohoo!

