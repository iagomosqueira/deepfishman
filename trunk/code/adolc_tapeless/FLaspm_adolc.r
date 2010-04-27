
library(FLCore)

# class
setClass('FLaspm', representation(
  'FLModel',
  catch='FLQuant',
  index='FLQuant',
  M='numeric',
  hh='numeric',
  sel='numeric',
  wght='numeric',
  mat='numeric',
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
# dyn.load('fit_adolc_tapeless.dll')

pdyn <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {

  ymin     <- dims(catch)$minyear
  ymax     <- dims(catch)$maxyear
  ymin_obj <- dims(index)$minyear
  ymax_obj <- dims(index)$maxyear
  
  B0     <- as.numeric(B0)
  sigma2 <- as.numeric(sigma2)
  catch  <- as.vector(catch)
  index  <- as.vector(index)
  hh     <- as.numeric(hh)
  M      <- as.vector(M)
  mat    <- as.vector(mat)
  sel    <- as.vector(sel)
  wght   <- as.vector(wght)
  amin   <- as.numeric(amin)
  amax   <- as.numeric(amax)
  
  out <- .Call('fit',B0,sigma2,catch,index,hh,M,m,sel,wght,amin,amax,ymin,ymax,ymin_obj,ymax_obj)
}

pdyn.index <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {

  FLQuant(pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['Ipred']],dimnames=dimnames(index))
}

aspm <- function() {

  # logl
  logl <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
    pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['nLogLk']]
  }
  
  # initial parameter values
  initial <- structure(function(catch) return(100*max(catch)),lower=c(1,1e-8), upper=c(Inf, Inf))

  model <- index ~ pdyn.index(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)
  
  return(list(logl=logl,model=model,initial=initial))
}

# post-fitting accessors for biomass etc.

# biomass {{{
if (!isGeneric("biomass"))
	setGeneric("biomass", function(object,type,...)
    	standardGeneric("biomass"))
setMethod('biomass', signature(object='FLaspm'),
  function(object,type='B') {
    return(FLQuant(pdyn(params(object)['B0'],params(object)['sigma2'],object@catch,object@index,object@hh,object@M,object@mat,object@sel,object@wght,object@amin,object@amax)[[type]],dimnames=dimnames(object@catch),units=units(object@catch)))
    })
# }}}

# harvest rate {{{
if (!isGeneric("harvest.rate"))
	setGeneric("harvest.rate", function(object, ...)
    	standardGeneric("harvest.rate"))
setMethod('harvest.rate', signature(object='FLaspm'),
  function(object) {
    return(FLQuant(pdyn(params(object)['B0'],params(object)['sigma2'],object@catch,object@index,object@hh,object@M,object@mat,object@sel,object@wght,object@amin,object@amax)[['F']],dimnames=dimnames(object@catch)))
    })
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

