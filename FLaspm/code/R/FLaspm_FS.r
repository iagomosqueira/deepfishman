
library(FLCore)

# Age structured parameters like sel and mat need to be FLQuant
# Shouldn't affect the model
# (also opens the possibility of year effects like sel changing over time which will require model changes)
# class
setClass('FLaspm', representation(
  'FLModel',
  catch='FLQuant',
  index='FLQuant',
#  M='numeric',
M = 'FLQuant', # So we can use multiple iterations.
# Maybe an FLPar would be more appropriate
# Could put in a single FLPar for M and hh?
#  hh='numeric',
hh = 'FLQuant',
  sel='FLQuant',
  wght='numeric',
  mat='FLQuant',
  fpm='numeric',
  amax='numeric',
  amin='numeric'
  )
)

# creator
setGeneric("FLaspm", function(model, ...){
		standardGeneric("FLaspm")
})
setMethod('FLaspm', signature(model='ANY'),
  function(model, ...)
  {
    res <- FLModel(model, ..., class='FLaspm')
    return(res)
  }
)
setMethod('FLaspm', signature(model='missing'),
  function(...)
  {
    res <- FLModel(..., class='FLaspm')
    return(res)
  }
)

# methods
# plot

# functions

# Population dynamics

aspm.pdyn <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax) {
#browser()
  C <- as.vector(catch)
  nyr <- length(C)
  nag <- amax-amin+1
  yr <- as.numeric(dimnames(catch)[['year']])
  iyr <- as.numeric(dimnames(index)[['year']])
  dm <- dimnames(catch)
  ind <- as.vector(index)
  ys <- iyr[1]
  yf <- iyr[length(iyr)]
  y1 <- (ys-yr[1])+1
  y2 <- (yf-yr[1])+1
  n <- array(dim=c(nag,nyr))
  b <- vector("numeric",length=nyr)
  bexp <- vector("numeric",length=nyr)
  h <- vector("numeric",length=nyr)
  p <- vector("numeric",length=nag)
  
  # set up eqm population
  p[1] <- 1
  for(a in 2:nag)
    p[a] <- p[a-1]*exp(-M)
  p[nag] <- p[nag]/(1-exp(-M))
  rho <- sum(p * mat * wght)
  R0 <- B0 / rho
  n[,1] <- R0 * p
  b[1] <- sum(n[,1] * mat * wght)
  bexp[1] <- sum(n[,1] * sel * wght)
  h[1] <- C[1] / bexp[1]
  h[1] <- max(h[1],0)
  h[1] <- min(h[1],0.999)

  # set up S-R parameters

  alp <- (4*hh*R0)/(5*hh-1)
  bet <- B0*(1-hh)/(5*hh-1)

  # Loop through the years

  for(y in 2:nyr) {

    # recruitment

    n[1,y] <- alp * b[y-1]/(bet + b[y-1])

    # adult dynamics

    for(a in 2:nag)
      n[a,y] <- n[a-1,y-1]*exp(-M)*(1-sel[a-1]*h[y-1])
    n[nag,y] <- n[nag,y] + n[nag,y-1]*exp(-M)*(1-sel[nag]*h[y-1])
    bexp[y] <- sum(n[,y] * sel * wght)
    h[y] <- C[y] / bexp[y]
    h[y] <- max(h[y],0)
    h[y] <- min(h[y],0.999)
    bexp[y] <- C[y] / h[y]
    b[y] <- sum(n[,y] * mat * wght)
  }
  # return predicted index
  return(FLQuant(bexp,dimnames=dm))
}

aspm.index <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax) {

  yr <- as.numeric(dimnames(catch)[['year']])
  iyr <- as.numeric(dimnames(index)[['year']])
  dm <- dimnames(catch)
  ind <- as.vector(index)
  ys <- iyr[1]
  yf <- iyr[length(iyr)]
  y1 <- (ys-yr[1])+1
  y2 <- (yf-yr[1])+1

#browser()

  # Strip out the FLQuant to speed it up
  mat <- c(mat)
  sel <- c(sel)
  
  bexp <- aspm.pdyn(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
  
  # nuisance q
  q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
  # return predicted index
  return(FLQuant(q*bexp,dimnames=dm))
}


aspm <- function() {

  # logl
  
  logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
  {
    #bexp<-ASPM.pdyn(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
    #hr<-catch/bexp
    #penalty<-100*length(which(hr >= 0.99))
    sum(dnorm(log(index),window(log(aspm.index(catch,index,B0,hh,M,mat,sel,wght,amin,amax)),start=dims(index)$minyear,end=dims(index)$maxyear), sqrt(sigma2), TRUE), na.rm=TRUE)#-penalty
  }
  
   # initial parameter values
   
  initial <- structure(function(catch) {
    return(FLPar(B0=100*max(catch), sigma2=1))
	},
  # lower and upper limits for optim()
	lower=c(1, 1e-8),
	upper=c(Inf, Inf)
	)



  model <- index ~ aspm.index(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
  
  return(list(logl=logl,model=model,initial=initial))
} # }}}

# post-fitting accessors for biomass etc.

# exploitable biomass {{{
if (!isGeneric("exp.biomass"))
	setGeneric("exp.biomass", function(object, ...)
    	standardGeneric("exp.biomass"))
setMethod('exp.biomass', signature(object='FLaspm'),
  function(object) {
    #return(fitted(object) / as.numeric(params(object)['q',]))
    return(aspm.pdyn(object@catch,object@index,params(object)['B0'],object@hh,object@M,object@mat,object@sel,object@wght,object@amin,object@amax))
    })
# }}}

# harvest rate {{{
if (!isGeneric("harvest.rate"))
	setGeneric("harvest.rate", function(object, ...)
    	standardGeneric("harvest.rate"))
setMethod('harvest.rate', signature(object='FLaspm'),
  function(object) {
    return(object@catch / exp.biomass(object))})
# }}}


# methods
setMethod('index', signature(object='FLaspm'),
  function(object)
    return(object@index)
)

if (!isGeneric("plot.fit"))
	setGeneric("plot.fit", function(object, ...)
    	standardGeneric("plot.fit"))
setMethod('plot.fit', signature(object='FLaspm'),
  function(object, ...) {
    windows(width=18)
    par(mfrow=c(1,3))
    plot(dimnames(object@index)$year,object@index,main='Fit to index',xlab='Year',ylab=paste('Index (',units(object@index),')',sep=''), ...)
    lines(dimnames(object@index)$year,object@fitted,lty=2, ...)
    plot(object@fitted,object@residuals,main='Residuals',ylab='Residuals',xlab=paste('Fitted Values (',units(object@index),')',sep=''), ...)
    abline(h=0,lty=2)
    qqnorm(as.vector(object@residuals),main='QQ plot of residuals', ...)
    qqline(as.vector(object@residuals),lty=2, ...)
  }
)
