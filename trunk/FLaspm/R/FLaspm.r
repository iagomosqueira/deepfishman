# Fix auropar
# Update plot
# write document with tests

# Then C code
# package it

#*******************************************************************************
# Validity
validFLaspm <- function(object)
{
#browser()

  # check that dim indices are the same as dim catch
  dim_catch <- dim(object@catch)
  for (i in 1:length(object@index))
    if(all(dim(object@index[[i]]) != dim_catch)) stop ("indices must have same dims as catch")

  # check year range of indices and catches are the same
  dnm_catch <- dimnames(object@catch)
  for (i in 1:length(object@index))
    if(!all(dimnames(object@index[[i]])$year %in% dnm_catch$year)) stop ("indices must have same year range as as catch")

# More checks about the age ranges of mat, sel, wght and amin and amax

  return(TRUE)
}

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
  sr_res = 'FLQuant',
  sel='FLQuant',
  wght='FLQuant',
  mat='FLQuant',
  fpm='numeric',
  amax='numeric',
  amin='numeric',
  fitted_index = 'FLQuants',
  fitted_flag = 'logical',
  residuals_index = 'FLQuants',
  pop.dyn = 'function',
  qhat = 'function'
  ),
  validity=validFLaspm
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

#setMethod('FLaspm', signature(model='missing'),
#  function(...)
#  {
#    #browser()
#    res <- FLModel(..., class='FLaspm')
#    return(res)
#  }
#)

setMethod('FLaspm', signature(model='missing'),
  function(...)
  {
    #browser()
    args <- list(...)
    # check through list for names (numeric and length 1)
    amin <- args[["amin"]]
    amax <- args[["amax"]]
    arg_names <- names(args)
    # Test for mat and sel
    for (slot in c("mat","sel"))
    {
      if (slot %in% arg_names & is.numeric(args[[slot]]))
      {
        if (length(args[[slot]]) == 1) 
        {
          temp <- FLQuant(0, dimnames=list(age=amin:amax))
          temp[(args[[slot]] - amin + 1):length(amin:amax),] <- 1
          args[[slot]] <- temp
        }
        if (length(args[[slot]]) > 1) 
        {
          temp <- FLQuant(args[[slot]], dimnames=list(age=amin:amax))
          args[[slot]] <- temp
        }
      }
    }
    
    # test for mean weight
    if ("wght" %in% arg_names & is.numeric(args[["wght"]]))
    {
        temp <- FLQuant(args[["wght"]], dimnames=list(age=amin:amax))
        args[["wght"]] <- temp
    }

    # Test for hh and M
    for (slot in c("hh","M"))
    {
      if (slot %in% arg_names & is.numeric(args[[slot]]) & !is.FLQuant(args[[slot]]) & length(args[[slot]]) == 1)
        args[[slot]] <- FLQuant(args[[slot]])
      if (slot %in% arg_names & is.numeric(args[[slot]]) & !is.FLQuant(args[[slot]]) & length(args[[slot]]) > 1)
        args[[slot]] <- FLQuant(args[[slot]],dim=c(1,1,1,1,1,length(args[[slot]])))
    }

    # Test for index
    if ("index" %in% arg_names & is.FLQuant(args[["index"]]))
      args[["index"]] <- FLQuants(index=args[["index"]])
      
    # use mcf to align dimensions?
    
    # default value for fpm if not specified
    if (!("fpm" %in% arg_names)) {
      args[["fpm"]] <- 1
      
    }
    
    # fill flags for fitting
    temp <- vector('logical',length(args[["index"]]))
    names(temp) <- names(args[["index"]])
    temp[] <- TRUE
    args[["fitted_flag"]] <- temp

#browser()

    # check sr_residuals
    if (!("sr_res" %in% arg_names) & "catch" %in% arg_names)
    {
      sr_res <- FLQuant(1,dimnames=dimnames(args[["catch"]]))
      args[["sr_res"]] <- sr_res
    }

    res <- do.call(FLModel,c(list(class='FLaspm'), args))



    #res <- FLModel(..., class='FLaspm')
    return(res)
  }
)


#********************************************************************************
# REWRITE pop.dyn so it evaluates the pop.dyn slot

# post-fitting accessors for biomass and fishing mortality etc.
# Need to fix this so that it works with Charlie and Francis models

# Does one iter at a time - bit crappy
# This is because pop.dyn function only handles one iter at a time
if (!isGeneric("calc.pop.dyn"))
    setGeneric("calc.pop.dyn", function(object, ...)
    standardGeneric("calc.pop.dyn"))


setMethod('calc.pop.dyn', signature(object='FLaspm'),
    function(object) {
  #browser()
	iters <- dims(object)$iter
	bexp <- FLQuant(NA,dimnames=dimnames(object@catch))
	bexp <- propagate(bexp,iters)
	bmat <- bexp
	harvest <- bexp
	n <- FLQuant(NA,dimnames=list(age=object@amin:object@amax,year = dimnames(object@catch)$year),iter=iters)

  # Sort out arguments for pop.dyn call
  parnames_not_in_params <- names(formals(object@pop.dyn)[names(formals(object@pop.dyn)) %in% slotNames(object)])
  args <- tapply(parnames_not_in_params, 1:length(parnames_not_in_params),function(x) slot(object,x), simplify=FALSE)
  names(args) <- parnames_not_in_params
  # add in params from params slot
  parnames_in_params <- names(formals(object@pop.dyn))[names(formals(object@pop.dyn)) %in% dimnames(object@params)$params]
  par_args <- tapply(parnames_in_params, 1:length(parnames_in_params),function(x) object@params[x], simplify=FALSE)
  names(par_args) <- parnames_in_params
  args <- c(args,par_args)

#browser()

	for (i in 1:iters)
	{
      iter_args <- lapply(args,function(x)iter(x,i))
      op <- do.call(object@pop.dyn,iter_args)
      iter(bexp,i)[] <- op[["bexp"]]
	    iter(bmat,i)[] <- op[["bmat"]]
	    iter(n,i)[] <- op[["n"]]
	    iter(harvest,i)[] <- op[["harvest"]]
	}
	return(list(bexp=bexp,bmat=bmat,n=n,harvest=harvest))
})

if (!isGeneric("calc.initial"))
    setGeneric("calc.initial", function(object, ...)
    standardGeneric("calc.initial"))


# Does one iter at a time - bit crappy
# This is because pop.dyn function only handles one iter at a time
setMethod('calc.initial', signature(object='FLaspm'),
  function(object) {

    iters <- dims(object)$iter
    init <- FLPar(NA,dimnames=dimnames(object@params),iter=iters)   # need iter argument?

    # Sort out arguments for call
    parnames_not_in_params <- names(formals(object@initial)[names(formals(object@initial)) %in% slotNames(object)])
    args <- tapply(parnames_not_in_params, 1:length(parnames_not_in_params),function(x) slot(object,x), simplify=FALSE)
    names(args) <- parnames_not_in_params

    for (i in 1:iters)
    {
      iter_args <- lapply(args,function(x)iter(x,i))
      iter(init,i)[] <- do.call(object@initial,iter_args)
    }
    
    return(init)
})


if (!isGeneric("calc.logl"))
    setGeneric("calc.logl", function(object, ...)
    standardGeneric("calc.logl"))


# Does one iter at a time - bit crappy
# This is because pop.dyn function only handles one iter at a time
setMethod('calc.logl', signature(object='FLaspm'),
  function(object) {

    #pop.dyn <- object@pop.dyn # not needed?
    iters <- dims(object)$iter
    logl <- FLQuant(NA,iter=iters)

    # Sort out arguments for pop.dyn call
    parnames_not_in_params <- names(formals(object@logl)[names(formals(object@logl)) %in% slotNames(object)])
  args <- tapply(parnames_not_in_params, 1:length(parnames_not_in_params),function(x) slot(object,x), simplify=FALSE)
  names(args) <- parnames_not_in_params
  # add in params from params slot
  parnames_in_params <- names(formals(object@logl))[names(formals(object@logl)) %in% dimnames(object@params)$params]
  par_args <- tapply(parnames_in_params, 1:length(parnames_in_params),function(x) object@params[x], simplify=FALSE)
  names(par_args) <- parnames_in_params
  args <- c(args,par_args)

	for (i in 1:iters)
	{
      iter_args <- lapply(args,function(x)iter(x,i))
      iter(logl,i)[] <- do.call(object@logl,iter_args)
	}
	return(logl)
})



if (!isGeneric("calc.qhat"))
    setGeneric("calc.qhat", function(object, ...)
    standardGeneric("calc.qhat"))

# Fix this for multiple iters
# Maybe the qhat function should work for multiple iters
# because here we are calling it repeatedly - pretty inefficient
setMethod('calc.qhat', signature(object='FLaspm'),
    function(object)
    {
      # call object@qhat with right arguments
      # e.g.
      #browser()
      parnames_not_in_params <- names(formals(object@qhat)[names(formals(object@qhat)) %in% slotNames(object)])
      args <- tapply(parnames_not_in_params, 1:length(parnames_not_in_params),function(x) slot(object,x), simplify=FALSE)
      names(args) <- parnames_not_in_params
      # add in params from params slot
      parnames_in_params <- names(formals(object@qhat))[names(formals(object@qhat)) %in% dimnames(object@params)$params]
      par_args <- tapply(parnames_in_params, 1:length(parnames_in_params),function(x) object@params[x], simplify=FALSE)
      names(par_args) <- parnames_in_params
      args <- c(args,par_args)
      # Set up output object
      op <- FLQuants()
      for (flq in 1:length(object@index))
        #op[[flq]] <- FLQuant(NA,iter=dims(object@index[[flq]])$iter)
        op[[flq]] <- FLQuant(NA,iter=dims(object)$iter)
      for (it in 1:dims(object)$iter)
      {
        iter.args <- lapply(args,function(x,it)iter(x,it),it=it)
        iter.op <- do.call(object@qhat,iter.args)
        for (flq in 1:length(object@index))
          iter(op[[flq]],it) <- iter.op[[flq]]
      }
      
      return(op)
  }
)

if (!isGeneric("calc.sigma2"))
    setGeneric("calc.sigma2", function(object, ...)
    standardGeneric("calc.sigma2"))
    
setMethod('calc.sigma2', signature(object='FLaspm'),
    function(object, yrfrac=0.5)
    {
      # check if sigma2 is in params slot, if so use that
      if ("sigma2" %in% dimnames(object@params)$params)
        return(object@params["sigma2",])
      # Otherwise we'll have to calculate it using Francis method
      bmid <- exp.biomass.mid(object,yrfrac=yrfrac)
      qhat <- calc.qhat(object)
      index <- object@index
      s2 <- FLPar(sigma2 = NA)
      s2 <- propagate(s2,length(index))
      for (index.count in 1:length(index))
      {
        nonnaindexyears <- !is.na(index[[index.count]])
        # number of non NA years in index
        n <- dim(index[[index.count]][nonnaindexyears])[2]
        s2["sigma2",index.count] <- apply((index[[index.count]] / sweep(bmid,1,qhat[[index.count]],"*") - 1)^2,c(1,6),sum,na.rm=T) / (n-2)
      }
      return(s2)
})

if (!isGeneric("indexhat"))
    setGeneric("indexhat", function(object, ...)
    standardGeneric("indexhat"))

setMethod('indexhat', signature(object='FLaspm'),
  function(object, yrfrac=0) {
    # yrfac should be 0 for Edwards and 0.5 for Francis
    #browser()
    bexp <- exp.biomass.mid(object,yrfrac)
    qhat <- calc.qhat(object)
    ihat <- lapply(qhat,function(x,b) sweep(b,c(1,3:6),x,"*"),b=bexp)
    return(ihat)
})

# exploitable biomass at some point through the year{{{
# Warning - this only makes sense if harvest has units of 'f'
# should put a check in but harvest needs to be done correctly elsewhere first
if (!isGeneric("exp.biomass.mid"))
    setGeneric("exp.biomass.mid", function(object, ...)
    standardGeneric("exp.biomass.mid"))

setMethod('exp.biomass.mid', signature(object='FLaspm'),
    function(object, yrfrac=0.5, virgin=T) {
      # if virgin - no fishing at harvest = 0
      pdyn <- calc.pop.dyn(object)
      bexp <- pdyn[["bexp"]]
      bmid <- bexp*exp(-yrfrac*sweep(as.numeric(!virgin)*pdyn[["harvest"]],c(1,3:6),object@M,"+"))
      return(bmid)
    })


# exploitable biomass {{{
if (!isGeneric("exp.biomass"))
    setGeneric("exp.biomass", function(object, ...)
    standardGeneric("exp.biomass"))

setMethod('exp.biomass', signature(object='FLaspm'),
    function(object) {
      return(calc.pop.dyn(object)[["bexp"]])
    })

# Numbers at age
if (!isGeneric("n"))
    setGeneric("n", function(object, ...)
    standardGeneric("n"))

setMethod('n', signature(object='FLaspm'),
    function(object) {
	return(calc.pop.dyn(object)[["n"]])
    })

# Mature biomass
if (!isGeneric("mat.biomass"))
    setGeneric("mat.biomass", function(object, ...)
    standardGeneric("mat.biomass"))

setMethod('mat.biomass', signature(object='FLaspm'),
    function(object) {
	return(calc.pop.dyn(object)[["bmat"]])
    })

# f
# Again multiple iters
if (!isGeneric('harvest'))
	setGeneric('harvest', function(object, ...)
    	standardGeneric('harvest'))

setMethod('harvest', signature(object='FLaspm'),
	function(object)
	{
    all <- calc.pop.dyn(object)
    return(all[["harvest"]])
	}
    )

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
# Johnson SU distribution and likelihood functions
# Used for estimating PDF of B0
#********************************************************************************
rJohnson <- function(n,parms)
{
  g <- parms["g"]
  delta <- parms["delta"]
  xi <- parms["xi"]
  lambda <- parms["lambda"]
  Z <- rnorm(n,mean=0,sd=1)
  B <- xi + lambda * sinh((Z-g)/delta)
  return(B)
}

JohnsonCum <- function(parms,b)
{
  g <- parms["g"]
  delta <- parms["delta"]
  xi <- parms["xi"]
  lambda <- parms["lambda"]
  q <- g + delta*asinh((b-xi)/lambda)
  pnorm(q,mean=0,sd=1)
}

JohnsonPDF <- function(parms,b)
{
  g <- parms["g"]
  delta <- parms["delta"]
  xi <- parms["xi"]
  lambda <- parms["lambda"]
  z <- (b - xi) / lambda
  return((delta / (lambda * sqrt(2*pi) * sqrt(z^2 + 1))) *
        exp(-0.5*(g + delta*asinh(z))^2))
#  q <- g + delta*asinh((b-xi)/lambda)
#  return(dnorm(q,mean=0,sd=1))
}

Johnsonll <- function(parms,b,m,p)
{
  #browser()
  g <- parms["g"]
  delta <- parms["delta"]
  xi <- parms["xi"]
  lambda <- parms["lambda"]
  jc <- JohnsonCum(parms,b)
  ll <- m * p * log(jc) + m * (1 - p) * log(1 - jc)
  return(-sum(ll))
}

#********************************************************************************
# Plot
#********************************************************************************
# Multiple iters?
# Currently only for the first index
setMethod('plot', signature(x='FLaspm'),
  function(x, ...) {
     
    ihat <- indexhat(x)

    par(mfrow=c(length(x@index),1))
    for(i in 1:length(x@index)) {
      x@fitted_index[[i]] <- ihat[[i]]
      
      y.rng1 <- range(x@index[[i]],na.rm=TRUE)
      y.rng2 <- range(x@fitted_index[[i]],na.rm=TRUE)
      y.rng  <- range(y.rng1,y.rng2)

      plot(dimnames(x@index[[i]])$year,x@index[[i]],ylim = y.rng,
              main=paste('Fit to index: ', names(x@index)[i]),xlab='Year',
              ylab=paste('Index (',units(x@index[[i]]),')',sep=''), ...)
      lines(dimnames(x@index[[i]])$year,x@fitted_index[[i]],lty=2, ...)
    
    }
    #plot(x@fitted_index[[1]],x@residuals_index[[1]],main='Residuals',ylab='Residuals',xlab=paste('Fitted Values (',units(x@index[[1]]),')',sep=''), ...)
    #abline(h=0,lty=2)
    #qqplot
    #qqnorm(as.vector(x@residuals_index[[1]]),main='QQ plot of residuals', ...)
    #qqline(as.vector(x@residuals_index[[1]]),lty=2, ...)
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






