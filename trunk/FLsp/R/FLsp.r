# To Do
# Multiple indices
# a plot
# iter check
# tests

# qhat should be FLQuants - for multiple indices

# Surplus production model class
# class FLsp
validFLsp <- function(object)
{
  # check that dim indices are the same as dim catch
  dim_catch <- dim(object@catch)
  for (i in 1:length(object@index))
    if(all(dim(object@index[[i]]) != dim_catch)) stop ("indices must have same dims as catch")

  # check year range of indices and catches are the same
  dnm_catch <- dimnames(object@catch)
  for (i in 1:length(object@index))
    if(!all(dimnames(object@index[[i]])$year %in% dnm_catch$year)) stop ("indices must have same year range as as catch")

  return(TRUE)
}

setClass('FLsp', representation(
  'FLModel',
  catch='FLQuant',
  #biomass='FLQuant',
  index='FLQuants',
  fitted_index='FLQuants',
  residuals_index='FLQuants'
),
	validity=validFLsp
)

setGeneric("FLsp", function(model, ...){
		standardGeneric("FLsp")
})

#setMethod('FLsp', signature(model='ANY'),
#  function(model, ...)
#  {
#		#browser()
#    res <- FLModel(model, ..., class='FLsp')
#    return(res)
#  }
#)
#
setMethod('FLsp', signature(model='missing'),
  function(...)
  {
    #browser()

    args <- list(...)
    catch <- args$catch
    index <- args$index
    if (class(index)=="FLQuant")
	index <- FLQuants(index=index)
    res <- FLModel(catch=catch,index=index, class='FLsp')



    #res <- FLModel(..., class='FLsp')
    model(res) <- sp
    # set up fitted_index and residuals index
    # Have to be the same dims as index slot - but empty
    # Why doesn't this work?
    #res@fitted_index <- lapply(res@index,function(x)x[] <- NA)
    res@residuals_index <- res@index
    for (i in 1:length(res@index))
	res@residuals_index[[i]][] <- NA
    res@fitted_index <- res@residuals_index


    return(res)
  }
)


#*******************************************************************************
# Only copes with one index at a time at the moment
sp <- function()
{
  logl <- function(r, k, catch, index)
  {
		res <- .Call("flspCpp",catch,index[[1]],r,1,k)[["ll"]]
		if (is.nan(res)) res <- Inf # Fix for DEoptim
  	if (is.na(res)) res <- Inf
  	return(res)
	}

  initial <- structure(function(catch)
	{
		# The function to provide initial values
		# DEoptim does not use start values
    return(FLPar(r=NA, k=NA))
	},

  # lower and upper limits for optim()
	lower=rep(1e-9, 2),
	upper=rep(1e9, 2))

	# awkward
	model  <- index ~ ihat(catch,index,r,k)

	return(list(logl=logl, model=model, initial=initial))
}

#*******************************************************************************
# Accessor things
ihat <- function(catch, index, r, k)
{
	res <- .Call("flspCpp",catch,index[[1]],r,1,k)
	ihat <- FLQuant(res[["qhat"]]*res[["B"]])
  indexhat_flqs <- FLQuants()
#  for (i in 1:length(index))
#    indexhat_flqs[[i]] <- FLQuant(indexhat_array[i,],dimnames=dimnames(catch))
	indexhat_flqs[[1]] <- ihat
  return(indexhat_flqs)
}

#*******************************************************************************
# Methods

# Fix iter to include FLQuants

if (!isGeneric("evalC"))
    setGeneric("evalC", function(object, ...)
    standardGeneric("evalC"))

setMethod('evalC', signature(object='FLsp'),
  function(object) {
      # Only works on first iter
	object <- iter(object,1)
	tape <- .Call("flspCpp",object@catch,object@index[[1]],object@params['r'],1,object@params['k'])
	return(tape)
})

#******** biomass *************
if (!isGeneric("biomass"))
    setGeneric("biomass", function(object, ...)
    standardGeneric("biomass"))

setMethod('biomass', signature(object='FLsp'),
  function(object) {
      #browser()
      # Make an FLQuant with right number of iterations
      iters <- dims(object)$iter
      dimnames <- dimnames(object@catch)
      dimnames$iter <- 1:iters
      biomass <- FLQuant(NA,dimnames=dimnames)
      for (i in 1:iters)
	  iter(biomass,i)[] <- evalC(iter(object,i))[["B"]]
      return(biomass)
})


#******** Bcurrent *************

if (!isGeneric("bcurrent"))
    setGeneric("bcurrent", function(object, ...)
    standardGeneric("bcurrent"))

setMethod('bcurrent', signature(object='FLsp'),
  function(object) {
      #biomass <- biomass(object)
      return(biomass(object)[,dim(object@catch)[2]])
})


if (!isGeneric("qhat"))
    setGeneric("qhat", function(object, ...)
    standardGeneric("qhat"))

setMethod('qhat', signature(object='FLsp'),
  function(object) {
      # Make an FLQuant with right number of iterations
      iters <- dims(object)$iter
      #dimnames <- dimnames(object@catch)
      #dimnames$iter <- iters
      qhat <- FLQuant(NA,iter=iters)
      for (i in 1:iters)
	  iter(qhat,i)[] <- evalC(iter(object,i))[["qhat"]]
      return(qhat)
})

if (!isGeneric("sigma2"))
    setGeneric("sigma2", function(object, ...)
    standardGeneric("sigma2"))

setMethod('sigma2', signature(object='FLsp'),
  function(object) {
      # Make an FLQuant with right number of iterations
      iters <- dims(object)$iter
      #dimnames <- dimnames(object@catch)
      #dimnames$iter <- iters
      sigma2 <- FLQuant(NA,iter=iters)
      for (i in 1:iters)
	  iter(sigma2,i)[] <- evalC(iter(object,i))[["sigma2"]]
      return(sigma2)
})


if (!isGeneric("ll"))
    setGeneric("ll", function(object, ...)
    standardGeneric("ll"))

setMethod('ll', signature(object='FLsp'),
  function(object) {
      # Make an FLQuant with right number of iterations
      iters <- dims(object)$iter
      #dimnames <- dimnames(object@catch)
      #dimnames$iter <- iters
      ll <- FLQuant(NA,iter=iters)
      for (i in 1:iters)
	  iter(ll,i)[] <- evalC(iter(object,i))[["ll"]]
      return(ll)
})


if (!isGeneric("indexhat"))
    setGeneric("indexhat", function(object, ...)
    standardGeneric("indexhat"))

setMethod('indexhat', signature(object='FLsp'),
  function(object) {
      # Make an FLQuant with right number of iterations
      nindex <- length(object@index)
      iters <- dims(object)$iter
      dimnames <- dimnames(object@catch)
      #dimnames$iter <- iters
      flq <- FLQuant(NA,dimnames=dimnames,iter=iters)
      indexhat <- FLQuants()
      for (i in 1:iters)
      {
	  ihat <- evalC(iter(object,i))[["Ihat"]]
	  for (j in 1:nindex)
	  {
	      iter(flq,i)[] <- ihat
	      indexhat[[j]] <- flq
	  }
      }
      return(indexhat)
})

if (!isGeneric("msy"))
    setGeneric("msy", function(object, ...)
    standardGeneric("msy"))

setMethod('msy', signature(object='FLsp'),
  function(object) {
      out <- c(t(params(object)['r',] * params(object)['k',] / 4))
      return(out)
})


#**** BMSY *****
if (!isGeneric("bmsy"))
    setGeneric("bmsy", function(object, ...)
    standardGeneric("bmsy"))

setMethod('bmsy', signature(object='FLsp'),
  function(object) {
      out <- c(params(object)['k',] / 2)
      return(out)
})

#**** Dims *****
# Need to overload this so that FLQuants slots are included
# And iters in params slot - pretty hacky...
# Why doesn't this get loaded?
setMethod("dims", signature(obj="FLsp"),
    # Returns a list with different parameters
    function(obj, ...)
	{
	    #browser()
    res <- callNextMethod()
    iters_in_index_slot <- max(unlist(lapply(obj@index,function(x)dim(x)[6])))
    iters_in_params_slot <- dim(obj@params)[2]
    res$iter <- max(res$iter,iters_in_index_slot, iters_in_params_slot)
    return(res)
	})


# New profile plot - includes gradients
setMethod("profile", signature(fitted="FLsp"),
  function(fitted, which, maxsteps=11, range=0.5, ci=c(0.25, 0.5, 0.75, 0.95),
      plot=TRUE, fixed=list(), print=FALSE, control=list(trace=0), ...)
  {
    # vars
    foo <- logl(fitted)
    params <- params(fitted)
    parnames <- dimnames(params)$params
    fixnames <- names(fixed)
    profiled <- list()
    grid <- list()
    plotfit <- TRUE

    # HACK! clean up fixed list if elements are named vectors
    fixed <- lapply(fixed, function(x){ names(x) <- NULL; x})

    # which params to profile
    if(missing(which))
      which <- parnames[!parnames %in% fixnames]
    if(length(which) > 2)
        stop("surface only works over 2 parameters")
    
    # data
    args <- list()
    data <- names(formals(foo))
    data <- data[data %in% slotNames(fitted)]
    for(i in data)
      args[i] <- list(slot(fitted, i))
      
    # use initial if model has not been estimated
    if(all(is.na(params)))
    {
      params <- do.call(initial(fitted), args)
      plotfit <- FALSE
    }

    # (1) create grid of param values for numeric range
    if(is.numeric(range) && length(range) == 1)
    {
      if(!plotfit)
        warning("model has not been fitted: initial values are used for profile range")
      for(i in which)
      {
        # steps for param[i]
        estim <- c(params[i,])
        steps <- estim * seq(1-range, 1+range, length=maxsteps)
        profiled[[i]] <- sort(steps)
      }
    # (2) and for list of ranges
    } else if (is.list(range)) 
    {
      # if missing(which), which is names in range
      if(missing(which))
        which <- names(range)
      else
        # checks all params to be profiled specified
        if(any(names(range) != which))
          stop("range not specified for parameters:", which[!which%in%names(range)])
      profiled <- lapply(range, sort)
    }

    # grid
    grid <- do.call(expand.grid, profiled)

    # col for logLik
    grid$logLik <- as.numeric(NA)

    dots <- list(...)
    # calculate logLik for grid if no fitting
    if(identical(order(c(which, fixnames)), order(parnames)))
      for(i in seq(nrow(grid)))
        grid[i, 'logLik'] <- do.call(logl(fitted), c(args, as.list(grid[i,which]), fixed))

    # or fit over grid
    else
      for(i in seq(nrow(grid)))
      {
        fixed <- as.list(grid[i,which])
        names(fixed) <- which
        grid[i, 'logLik'] <- do.call('fmle', c(list(object=fitted, fixed=fixed,
          control=control), dots))@logLik
      }
   
    surface <- tapply(grid$logLik, grid[,which], sum)

    # print
    if(print)
    {
      cat(paste("max(profile) =", format(max(grid$logLik), digits=5), " "))
      for(i in which)
        cat(paste(i, " = ", format(grid[grid$logLik==max(grid$logLik),i], digits=5), " "))
      cat("\n")
      if(plotfit)
      {
        cat(paste("logLik =", format(logLik(fitted), digits=5), " "))
        for(i in which)
          cat(paste(i, " = ", format(c(params(fitted)[i]), digits=5), " "))
        cat("\n")
      }
    }

    # CIs
    cis <- max(surface) - qchisq(ci, 2)
    
    # plot
    if(plot)
    {
      if(length(which) == 2)
      {

    profiled[["ll_grad_k"]] <- rep(NA,length(profiled[[1]]))
    profiled[["ll_grad_r"]] <- rep(NA,length(profiled[[1]]))
    dummy <- fitted
    dummy@params <- params
    # get dll/dr
    for (i in 1:length(profiled[["ll_grad_k"]]))
    {
	dummy@params["r",] <- profiled[["r"]][i]
	profiled[["ll_grad_r"]][i] <- evalC(iter(dummy,1))[["ll_grad_r"]]
    }
    dummy@params <- params
    # get dll/dk
    for (i in 1:length(profiled[["ll_grad_k"]]))
    {
	dummy@params["k",] <- profiled[["k"]][i]
	profiled[["ll_grad_k"]][i] <- evalC(iter(dummy,1))[["ll_grad_k"]]
    }
    #browser()
    # Need y and x lims for the plots
    ylim= profiled[["ll_grad_r"]][round(maxsteps/2)+1]*c(-3,3)
    xlim= profiled[["ll_grad_k"]][round(maxsteps/2)+1]*c(-3,3)

    # Now set up plot
    lay_mat <- matrix(c(1,3,0,2),nrow=2)
    lay<-layout(lay_mat,widths=c(3,1),heights=c(1,3))
    #layout.show(lay)
    # And fanny about with margins
    par(mar=c(0,5.1,1,1))
    plot(x=profiled[["r"]], y=profiled[["ll_grad_r"]], type="l", axes=FALSE, xlab="",ylab="dll/dr", ylim=ylim)
    lines(x=c(-1e9,1e9),y=c(0,0),lty=2)
    par(mar=c(5.1,0,1,1))
    plot(y=profiled[["k"]], x=profiled[["ll_grad_k"]], type="l", axes=FALSE, ylab= "", xlab="dll/dr",xlim=xlim)
    lines(y=c(-1e9,1e9),x=c(0,0),lty=2)
    par(mar=c(5.1,5.1,1,1))


        do.call('image', c(list(x=profiled[[1]], y=profiled[[2]], z=surface,
          xlab=which[1], ylab=which[2]), dots[!names(dots) %in% names(formals(optim))]))

        if(plotfit)
          points(params[which[1]], params[which[2]], pch=19)

        do.call('contour', list(x=sort(profiled[[1]]), y=sort(profiled[[2]]), z=surface,
          levels=cis, add=TRUE, labcex=0.8, labels=ci))
      }
      else if(length(which) == 1)
      {
        plot(grid[,which], grid[,'logLik'], type='l', xlab=which, ylab="logLik", axes=F)
        axis(1); box()
        points(params[which], logLik(fitted), pch=19)
      }
    }
    if(length(which) == 2)
      invisible(list(x=grid[,which[1]], y=grid[,which[2]], z=surface))
    else if(length(which) == 1)
      invisible(list(x=grid[which], y=grid['logLik']))
  }
) # }}}






