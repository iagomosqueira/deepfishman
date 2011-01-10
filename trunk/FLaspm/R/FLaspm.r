# Fix auropar
# Update plot
# write document with tests

# Then C code
# package it

#*******************************************************************************
#library(FLCore)
# Age structured parameters like sel and mat need to be FLQuant
# Shouldn't affect the model
# (also opens the possibility of year effects like sel changing over time which will require model changes)
# class
setClass('FLaspm', representation(
  'FLModel',
  catch='FLQuant',
  index='FLQuants',
  M = 'FLQuant',
  hh = 'FLQuant',
  sel='FLQuant',
  wght='FLQuant',
  mat='FLQuant',
  fpm='numeric',
  amax='numeric',
  amin='numeric',
  fitted_index = 'FLQuants',
  residuals_index = 'FLQuants'
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
    #browser()
    res <- FLModel(..., class='FLaspm')
    return(res)
  }
)

#setMethod('FLaspm', signature(model='ANY'),
#  function(model, ...)
#  {
#      # Pull out args
#      args <- list(...)
#      #res <- FLModel(model, ..., class='FLaspm')
#      res <- FLaspm(...)
#      model(res) <- model
#    return(res)
#  }
#)
# 
#
## More useful creator
## This is called if model is missing. Can we specify model too?
## Pass in selage, matage, M, hh, amin, amax, Linf, K, t0
## and catch and indices
#
# Calling FLModel with missing model, calls as if with formula
# All FLArray slots then set to same dims
# Including fitted and residual which is wrong.
#setMethod('FLaspm', signature(model='missing'),
#  function(...)
#  {
#      #browser()
#	args <- list(...)
#	# Make age based quants based on amin and amax
#	if (all(c('amin','amax') %in% names(args)))
#	{
#	    amin <- args[['amin']]
#	    amax <- args[['amax']]
#	    age.quant <- FLQuant(NA,dimnames=list(age=amin:amax))
#	    # Age based quants based on this
#	    sel <- age.quant
#	    wght <- age.quant
#	    mat <- age.quant
#	    # Fill these up if you have them
#	    if ('selage' %in% names(args))
#	    {
#		sel[] <- 0
#		sel[ac(args[['selage']]:amax),] <- 1  
#	    }
#	    if ('matage' %in% names(args))
#	    {
#		mat[] <- 0
#		mat[ac(args[['matage']]:amax),] <- 1  
#	    }
#	    if (all(c('Linf','k','t0','a','b') %in% names(args)))
#		wght[] <- args[['a']] * (args[['Linf']] * (1 - exp(-args[['k']] * ((amin:amax) - args[['t0']]))))^args[['b']]
#	    res <- FLModel(amin=amin, amax=amax, sel=sel, wght=wght, mat=mat,class='FLaspm')
#	}
#	else
#	    res <- FLModel(..., class='FLaspm')
#
#	#browser()
#	# Non-age structured Quants
#	# This is horrible
#	if('M' %in% names(args)) res@M <- M
#	else res@M <- FLQuant()
#	if('hh' %in% names(args)) res@hh <-hh 
#	else res@hh <- FLQuant()
#	if('catch' %in% names(args)) res@catch <- catch
#	else res@catch <- FLQuant()
#	if('index' %in% names(args)) 
#	{
#	    if (!is.FLQuants(index)) index <- FLQuants(index)
#	    res@index <- index   
#	}
#	else res@index <- FLQuants()
#
#	    # More checks, make catch and index have same year range
#
#	return(res)
#  }
#)
#


#********************************************************************************
# post-fitting accessors for biomass and fishing mortality etc.
# Need to fix this so that it works with Charlie and Francis models

if (!isGeneric("pop.dyn"))
    setGeneric("pop.dyn", function(object, ...)
    standardGeneric("pop.dyn"))

# Also what are all those cs about?
# This is pretty crappy, doing one iter at a time
setMethod('pop.dyn', signature(object='FLaspm'),
    function(object) {
	#browser()
	iters <- dims(object)$iter
	bexp <- FLQuant(NA,dimnames=dimnames(object@catch))
	bexp <- propagate(bexp,iters)
	bmat <- bexp
	harvest <- bexp
	n <- FLQuant(NA,dimnames=list(age=object@amin:object@amax,year = dimnames(object@catch)$year),iter=iters)
	# Need to make sure it calls the right function
	# This is pretty ugly and slow.
	if(grepl("Francis",as.character(model(object))[3]))
	    pdyn.func <- "aspm.pdyn.Francis"
	if(grepl("Edwards",as.character(model(object))[3]))
	    pdyn.func <- "aspm.pdyn.Edwards"
	for (i in 1:iters)
	{
	    # Need to make sure it calls the right function
	    op <- eval(call(pdyn.func,iter(object@catch,i),iter(params(object)['B0'],i),c(iter(object@hh,i)),c(iter(object@M,i)),c(iter(object@mat,i)),c(iter(object@sel,i)),c(iter(object@wght,i)),object@amin,object@amax))
	    #op <- aspm.pdyn(iter(object@catch,i),iter(params(object)['B0'],i),c(iter(object@hh,i)),c(iter(object@M,i)),c(iter(object@mat,i)),c(iter(object@sel,i)),c(iter(object@wght,i)),object@amin,object@amax)
	    iter(bexp,i) <- op[["bexp"]]
	    iter(bmat,i) <- op[["bmat"]]
	    iter(n,i) <- op[["n"]]
	    iter(harvest,i) <- op[["harvest"]]
	}
	return(op)
    })

# exploitable biomass {{{
if (!isGeneric("exp.biomass"))
    setGeneric("exp.biomass", function(object, ...)
    standardGeneric("exp.biomass"))

# Is this going to work with multiple iters?
# Probably not
# Also what are all those cs about?
setMethod('exp.biomass', signature(object='FLaspm'),
    function(object) {
	#iters <- dims(object)$iter
	#bexp <- FLQuant(NA,dimnames=dimnames(object@catch))
	#bexp <- propagate(bexp,iters)
	#for (i in 1:iters)
	#    iter(bexp,i) <- aspm.pdyn(iter(object@catch,i),iter(params(object)['B0'],i),c(iter(object@hh,i)),c(iter(object@M,i)),c(iter(object@mat,i)),c(iter(object@sel,i)),c(iter(object@wght,i)),object@amin,object@amax)[["bexp"]]
	#return(bexp)
	return(pop.dyn(object)[["bexp"]])
    })

# Numbers at age
if (!isGeneric("n"))
    setGeneric("n", function(object, ...)
    standardGeneric("n"))

setMethod('n', signature(object='FLaspm'),
    function(object) {
	return(pop.dyn(object)[["n"]])
    })

# Mature biomass
if (!isGeneric("mat.biomass"))
    setGeneric("mat.biomass", function(object, ...)
    standardGeneric("mat.biomass"))

setMethod('mat.biomass', signature(object='FLaspm'),
    function(object) {
	return(pop.dyn(object)[["bmat"]])
    })

# f
# Again multiple iters
if (!isGeneric('harvest'))
	setGeneric('harvest', function(object, ...)
    	standardGeneric('harvest'))

setMethod('harvest', signature(object='FLaspm'),
	function(object)
	{
	    #bexp <- exp.biomass(object)
	    #iter <- dims(object)$iter
	    #f <- FLQuant(NA,dimnames=dimnames(bexp))
	    #for (i in 1:iter)
	    #    for (y in 1:dim(object@catch)[2])
	    #        f[,y,,,,i] <- ratner_search(func=fobj,x=c(0,1,100),m=c(iter(object@M,i)),catch=c(iter(object@catch,i))[y],biomass=c(iter(bexp,i))[y])
	    #return(f)
	    all <- pop.dyn(object)
	    return(all[["harvest"]])
	}
    )

# harvest rate {{{
#if (!isGeneric("harvest.rate"))
#	setGeneric("harvest.rate", function(object, ...)
#    	standardGeneric("harvest.rate"))
#setMethod('harvest.rate', signature(object='FLaspm'),
#  function(object) {
#    return(object@catch / exp.biomass(object))})
## }}}


# methods
setMethod('index', signature(object='FLaspm'),
  function(object)
    return(object@index)
)

#********************************************************************************
# VB growth functions
#********************************************************************************
age_to_length <- function(age,Linf,k,t0)
  return(Linf * (1 - exp(-k * (age - t0))))

length_to_weight <- function(l,a,b)
  return(a * l^b)
  
age_to_weight <- function(age,Linf,k,t0,a,b)
  return(length_to_weight(age_to_length(age,Linf,k,t0),a,b))

#********************************************************************************
# Plot
#********************************************************************************
# Multiple iters?
# Currently only for the first index
setMethod('plot', signature(x='FLaspm'),
  function(x, ...) {

      # If not actually fitted, need to calc the residuals etc here


    par(mfrow=c(3,1))
    plot(dimnames(x@index[[1]])$year,x@index[[1]],main='Fit to index',xlab='Year',ylab=paste('Index (',units(x@index[[1]]),')',sep=''), ...)
    lines(dimnames(x@index[[1]])$year,x@fitted_index[[1]],lty=2, ...)
    plot(x@fitted_index[[1]],x@residuals_index[[1]],main='Residuals',ylab='Residuals',xlab=paste('Fitted Values (',units(x@index[[1]]),')',sep=''), ...)
    abline(h=0,lty=2)
    #qqplot
    qqnorm(as.vector(x@residuals_index[[1]]),main='QQ plot of residuals', ...)
    qqline(as.vector(x@residuals_index[[1]]),lty=2, ...)
  }
)

#********************************************************************************
# Functions for using the C version - no gradients
# CHECK and make worth for Francis and Charlie
#********************************************************************************
aspm_c <- function() {
    # logl
    logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
    {
	nyrs <- length(c(catch))
	#cat("Before .Call\n")
	out <- .Call("aspm",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,nyrs)
	#cat("After .Call\n")
	total.logl <- out[["logl"]][["logl"]]
	return(total.logl)
    }
  
    # initial parameter values
    initial <- structure(
	function(catch) {
	return(FLPar(B0=100*max(catch), sigma2=1))
	},
	# lower and upper limits for optim()
	lower=c(1, 1e-8),
	upper=c(Inf, Inf)
    )

    # The model function returns index_hat
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

# Runs the projection based on B0 in params slot and returns the output
get_aspm_c <- function(object)
{
    #browser()
    out_c <-  .Call("aspm",(c(object@catch)),object@index,object@params["B0"],object@params["sigma2"],(c(object@hh)),(c(object@M)),(c(object@mat)),(c(object@sel)),(c(object@wght)),object@amin,object@amax,length(c(object@catch)))
    dms <- dimnames(object@catch)
    Bexp <- FLQuant(out_c$Bexp,dimnames=dms)
    B <- FLQuant(out_c$B,dimnames=dms)
    h <- FLQuant(out_c$h,dimnames=dms)
    indexhat <- FLQuants()
    for (i in 1:length(object@index))
	indexhat[[i]] <- FLQuant(out_c$indexhat[i,],dimnames=dms)
    out = list(Bexp = Bexp, B=B, f=h, q = out_c$q, logl=out_c$logl, indexhat=indexhat)
    return(out)
}


#********************************************************************************
# Functions for using the AD version - NOT READY YET!
#********************************************************************************
aspm_ad <- function() {
    # logl
    logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
    {
	nyrs <- length(c(catch))
	#cat("Before .Call\n")
	out <- .Call("aspm_ad",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,nyrs)
	#cat("After .Call\n")
	total.logl <- out[["logl"]][["logl"]]
	return(total.logl)
    }
  
    # initial parameter values
    initial <- structure(
	function(catch) {
	return(FLPar(B0=100*max(catch), sigma2=1))
	},
	# lower and upper limits for optim()
	lower=c(1, 1e-8),
	upper=c(Inf, Inf)
    )

    # The model function returns index_hat
    model <- index ~ get_indexhat(catch,index,B0,sigma2,hh,M,mat,sel,wght,amin,amax)

    return(list(logl=logl,model=model,initial=initial))
} # }}}

get_indexhat <- function(catch,index,B0,sigma2,hh,M,mat,sel,wght,amin,amax)
{
    #browser()
    out <-  .Call("aspm_ad",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,length(c(catch)))[["indexhat"]]
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
    out <-  .Call("aspm_ad",(c(catch)),index,B0,sigma2,(c(hh)),(c(M)),(c(mat)),(c(sel)),(c(wght)),amin,amax,length(c(catch)))[["logl"]]
    # Return * -1 because 
    return(-1*c(out[["logl_grad_B0"]],out[["logl_grad_sigma2"]]))
}


# Runs the projection based on B0 in params slot and returns the output
get_aspm_ad <- function(object)
{
    #browser()
    out_ad <-  .Call("aspm_ad",(c(object@catch)),object@index,object@params["B0"],object@params["sigma2"],(c(object@hh)),(c(object@M)),(c(object@mat)),(c(object@sel)),(c(object@wght)),object@amin,object@amax,length(c(object@catch)))
    dms <- dimnames(object@catch)
    Bexp <- FLQuant(out_ad$Bexp,dimnames=dms)
    B <- FLQuant(out_ad$B,dimnames=dms)
    h <- FLQuant(out_ad$h,dimnames=dms)
    indexhat <- FLQuants()
    for (i in 1:length(object@index))
	indexhat[[i]] <- FLQuant(out_ad$indexhat[i,],dimnames=dms)
    out = list(Bexp = Bexp, B=B, f=h, q = out_ad$q, logl=out_ad$logl, indexhat=indexhat)
    return(out)
}






