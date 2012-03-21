# Method for fitting FLsp objects using DEoptim
# Stolen wholesale from fmle()

setGeneric('fitsp', function(object, ...)
    standardGeneric('fitsp'))


# Sort out control
setMethod('fitsp',
	signature(object="FLsp"),
	function(object, fixed=list(),
    #control = DEoptim.control(trace=50),
    control = DEoptim.control(NP=50,trace=200,itermax=2000),
		lower=NULL,
		upper=NULL, start=missing, seq.iter=TRUE, ...)
	{
	    #browser()

    args <- list(...)
    call <- sys.call(1)
    logl <- object@logl

    parnm <- names(formals(logl))[names(formals(logl))%in%
			dimnames(object@params)$param]

		# Can we handle fixed parameters?
    # get fixed parameter names
    fixnm <- names(fixed)
    # fixed must match params
    if(any(!fixnm %in% parnm))
      stop("some named arguments in 'fixed' are not arguments to the
        supplied log-likelihood")
    # HACK! clean up fixed list if elements are named vectors
    fixed <- lapply(fixed, function(x){ names(x) <- NULL; x})

		# This could be simplified because we know what slots we need
#    datanm <- getSlotNamesClass(object, 'FLArray')
#    # Include FLQuants contents too
#		flqs <- getSlotNamesClass(object, 'FLQuants')
#    for (i in length(flqs))
#      datanm <- c(datanm, names(slot(object, flqs[i])))
#    datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))
#    #   get those in formals of logl
#    datanm <- datanm[datanm%in%names(formals(logl))]
		datanm <- c("catch","index")

		# limits
		if (is.null(lower))
        lower <- lower(object)[match(parnm, names(fixed), nomatch=0)==0]
    if(is.null(upper))
        upper <- upper(object)[match(parnm, names(fixed), nomatch=0)==0]

		# gr function
		if(!is.null(body(object@gr)))
    {
      gr <- function(par)
      {
      #browser()
        pars <- as.list(par)
        names(pars) <- names(start)
    	# We are estimating log(params) so we need to pass in exp(params) to log function
	    	pars <- lapply(pars,exp)
        pars[fixnm] <- lapply(fixed, iter, it)
				#cat("Calling gradient function\n")
				grs <- -1*(do.call(object@gr, args=c(pars, data)))
				#cat("Gradients ", grs, "\n")
        return(grs)
      }
    }
    else
		gr <- NULL
		
	  loglfoo <- function(par) {
	    #browser()
    	pars <- as.list(par)
    	names(pars) <- names(start)
    	# We are estimating log(params) so we need to pass in exp(params) to log function
    	pars <- lapply(pars,exp)
    	pars[fixnm] <- lapply(fixed, iter, it)
    	return(-1*(do.call(logl, args=c(pars, data))))
    }

    # input data
    alldata <- list()
    for (i in datanm)
      alldata[[i]] <- slot(object, i)

    dimna <- dimnames(slot(object, datanm[1]))[names(slot(object, datanm[1]))%in%
      all.vars(object@model)]
    if(length(dimna) > 0)
    {
      # get them in the right shape
      dimdat <- lapply(dimna, function(x)
        {
          out <- slot(object, datanm[1])
          out[] <- as.numeric(x)
          return(out)
        })
      alldata <- c(alldata, dimdat)
    }

	if(seq.iter)
    {
      iter <- dims(object)$iter
      # Problem in that dims doesn't include the dims in index

      # iters in fixed
      if(length(fixnm) >= 1)
      {
        fiter <- unlist(lapply(fixed, length))
        if(!all(fiter == 1))
        {
          fiter <- fiter[fiter > 1]
          # all ietrs in fixed are equal?
          if(any(fiter/fiter[1] != 1))
            stop("objects in fixed have different number of iters")
          # are iter in object 1 and fixiter > 1? use fixiter
          if(iter == 1 & fiter > 1)
            iter <- fiter
          # are they different and > 1? STOP
          else if(fiter > 1 & fiter != iter)
            stop("different iters in fixed and object")
        }
      }
    }
    else
      iter <- 1

	  logLik <- rep(NA, iter)
    class(logLik) <- 'logLik'
    attr(logLik, 'df') <- length(parnm) - length(fixed)
    object@logLik <- logLik

    # Correct FLPar, fitted and residuals
    if(iter > dim(object@params)[length(dim(object@params))])
    {
      params(object) <- FLPar(iter=iter, params=dimnames(object@params)$params)
    }

    fitted(object) <- propagate(fitted(object), iter)
    residuals(object) <- propagate(residuals(object), iter)

    # vcov
    #object@vcov <- array(NA, dim=c(rep(length(parnm)-length(fixed),2), iter),
    #  dimnames=list(parnm[!parnm%in%names(fixed)],parnm[!parnm%in%names(fixed)],
    #  iter=1:iter))
    #object@hessian <- object@vcov
    object@hessian <- array(NA,dim=c(2,2,iter),dimnames=list(c("r","k"),c("r","k"),iter=1:iter))
    object@vcov <- object@hessian

    object@hessian_log <- array(NA,dim=c(2,2,iter),dimnames=list(c("r","k"),c("r","k"),iter=1:iter))
    object@vcov_log <- object@hessian_log


    #browser()
    # We're solving on the log scale so upper and lower need to be logged
	    lower <- log(lower)
	    upper <- log(upper)

# Not picked up iters
# Set up fitted_index and residual_index slots
for (index.count in 1:length(object@index))
{
    index_dmns <- dimnames(object@index[[index.count]])
    object@fitted_index[[index.count]] <- propagate(FLQuant(NA,dimnames=index_dmns),iter)
    object@residuals_index[[index.count]] <- propagate(FLQuant(NA,dimnames=index_dmns),iter)
}

#browser()
    for (it in 1:iter)
    {

#if (it == 4) browser()

    cat("iter: ", it, "\n")

      # data
      if(seq.iter)
        data <- lapply(alldata, iter, it)
      else
        data <- alldata

#print(data)

	# We don't have start values but we need to set some because they are used by loglfoo
      # start values
      if(missing(start)) {
#        # add call to @initial
        if(is.function(object@initial))
         start <- as(do.call(object@initial, args=data[names(formals(object@initial))]),
           'list')
        else
          start <- formals(logl)[names(formals(logl))%in%parnm]
      }
      else
#        # HACK! clean up fixed list if elements are named vectors
        start <- lapply(start, function(x){ names(x) <- NULL; x})

     if(!is.null(fixnm))
        start[fixnm] <- NULL
      if(any(!names(start) %in% parnm))
        stop("some named arguments in 'start' are not arguments to the
          supplied log-likelihood")
      start <- start[order(match(names(start), parnm))]
      # add small number to start if 0
      start <- lapply(start, function(x) if(x == 0) x/100000 else x)
      if(is.null(start))
        stop("No starting values provided and no initial function available")

      # TODO protect environment
#      out <- do.call('optim', c(list(par=unlist(start), fn=loglfoo, method=method,
#        hessian=TRUE, control=control, lower=lower, upper=upper, gr=gr)))

    #browser()
			# Using DEoptim
        out <- do.call('DEoptim', c(list(fn=loglfoo, lower=lower, upper=upper, control=control)))
        names(out$optim$bestmem) <- names(start)
        iter(object@params[names(out$optim$bestmem),], it) <- exp(out$optim$bestmem)
        object@logLik[it] <- -out$optim$bestval
        cat("DEoptim best: ", exp(out$optim$bestmem), "\n")

        # Finish off with optim? - Doesn't like infs
        #			cat("Trying with optim\n")
        #      out <- do.call('optim', c(list(par=exp(out$optim$bestmem), fn=loglfoo, method="BFGS",
        #        hessian=TRUE, gr=gr)))


# Use ucminf to finish off
#browser()
#cat("Trying ucminf\n")
#require(ucminf)
#out <- do.call('ucminf', c(list(par=out$optim$bestmem, fn=loglfoo, gr=gr)))
#cat("ucminf out\n")
#print(exp(out$par))
#print(out$value)


			# Use Rgenoud
			#browser()
#			out <- do.call('genoud', c(list(fn=loglfoo, nvars= (2-length(fixed)),
#											control=control, gr=gr,
#											print.level=1,
#											max.generations=200,
#											#hard.generation.limit=FALSE,
#											starting.values = log(unlist(start)),
#											#starting.values = c(-1.23,6.3),
#											Domains = matrix(c(log(1e-9), log(10), log(1e-9), log(5000)),nrow=2,byrow=TRUE)))
#											)
#	    names(out$par) <- names(start)
#      iter(object@params[names(out$par),], it) <- exp(out$par)
#      object@logLik[it] <- -out$value
#			#browser()
#

			#out <- do.call('genoud', c(list(fn=loglfoo, nvars= (2-length(fixed)), control=control, gr=gr)))
	    #names(out$par) <- names(start)
      #iter(object@params[names(out$par),], it) <- exp(out$par)
      #object@logLik[it] <- -out$value



			# fixed
      if(length(fixed) > 0)
        iter(object@params, it)[fixnm,] <- unlist(lapply(fixed, iter, it))
      # TODO make details list of lists if iter > 1?
      #object@details <- list(call=call, value=out$value, count=out$counts,
      #  convergence=out$convergence, message=out$message)

      # logLik
      attr(object@logLik, 'nobs') <- length(data[[1]])




    # fitted & residuals
    # No iter <- methods for FLQuants so a bit hacky
#browser()
	for (i in 1:length(object@index))
	{
	    object@fitted_index[[i]][,,,,,it] <- predict(iter(object, it))[[i]]
	    object@residuals_index[[i]][,,,,,it] <- iter(object@index[[i]],it) - iter(object@fitted_index[[i]], it)
	}

	# Load up the hessian slots
	# leave out for the moment
	#browser()
	tape_res <- .Call("flspCpp_tape",iter(object@catch,it),iter(object@index[[1]],it),iter(object@params["r"],it),1,iter(object@params["k"],it))
	# fix the upper right part of hessian
	tape_res$hessian[1,2] <- tape_res$hessian[2,1]
	object@hessian[,,it] <- tape_res$hessian

 	tape_res_log <- .Call("flspCpp_tape_log",iter(object@catch,it),iter(object@index[[1]],it),iter(object@params["r"],it),1,iter(object@params["k"],it))
 	# fix the upper right part of hessian
	tape_res_log$hessian[1,2] <- tape_res_log$hessian[2,1]
	object@hessian_log[,,it] <- tape_res_log$hessian

	# Sort out variance-covariance matrix
	tempvcov <- try(solve(-1 * object@hessian[,,it]),silent=TRUE)
	if (class(tempvcov) == 'try-error')
	    object@vcov[,,it] <- NA
	else
	    object@vcov[,,it] <- tempvcov

	tempvcov_log <- try(solve(-1 * object@hessian_log[,,it]),silent=TRUE)
	if (class(tempvcov_log) == 'try-error')
	    object@vcov_log[,,it] <- NA
	else
	    object@vcov_log[,,it] <- tempvcov_log


    }
    # force dimnames[1:5] in 'fitted' and 'residuals' to match
    #dimnames(fitted(object))[1:5] <- dimnames(do.call(as.character(
    #  as.list(object@model)[2]), list(object)))[1:5]
    #dimnames(residuals(object)) <- dimnames(fitted(object))

    # return object
    return(object)

	}
)

#*******************************************************************************
# Also overloading predict
#setMethod('predict', signature(object='FLsp'),
#  function(object, ...)
#  {
#      #browser()
##stop("in predict for FLaspm")
#    args <- list(...)
#    if(length(args) > 0 && is.null(names(args)))
#      stop('FLQuant or FLCohort inputs must be named to apply formula')
#    # call
#    call <- as.list(object@model)[[3]]
#    fittedSlot <- as.list(object@model)[[2]]
#
#
#    # check vars in call match input in args
#    if(length(args) > 0 & !any(names(args)%in%all.vars(call)))
#      warning(paste("Input names do not match those in model formula: '",
#        paste(names(args)[!names(args)%in%all.vars(call)], collapse=','), "'", sep=""))
#
#    # create list of input data
#    #   get FLQuant/FLCohort slots' names
#    datanm <- getSlotNamesClass(object, 'FLArray')
#    datanm <- c(datanm, getSlotNamesClass(object, 'FLQuants'))
#    datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))
#
#    # add dimnames if used
#    dimna <- dimnames(slot(object, datanm[1]))[names(slot(object, datanm[1]))%in%
#      all.vars(object@model)]
#    # get them in the right shape
#    dimdat <- lapply(dimna, function(x)
#      {
#        out <- slot(object, datanm[1])
#        out[] <- as.numeric(x)
#        return(out)
#      })
#
#    # iterations
#    #   from object
#    iter <- max(unlist(qapply(object, function(x) dims(x)$iter)))
#    #   from extra input
#    if(length(args) > 0)
#    {
#      iterarg <- lapply(args, function(x) {
#        itera <- try(dims(x)$iter)
#        if(class(iter) =='try-error')
#          return(1)
#        else
#          return(itera)
#      })
#      iterarg <- max(unlist(iterarg))
#    }
#    else
#      iterarg <- 1
#    #   decision
#    if (iter == iterarg)
#      iters <- iter
#    else if(iter > iterarg && iterarg == 1)
#      iters <- iter
#    else if(iterarg > iter && iter == 1)
#      iters <- iterarg
#    else
#      stop("Iter for object and input arguments do not match")
#
#    fitted_index <- FLQuants()
#    for (index.count in 1:length(object@index))
#    {
#
#      for (it in 1:iters)
#      {
#      obj <- iter(object, it)
#
#      #   input data
#        data <- list()
#        for (i in datanm)
#          data[[i]] <- slot(obj, i)
#
#        # add covar if defined and available
#        if('covar' %in% slotNames(obj))
#        {
#          covarnm <- names(obj@covar)
#          if(length(covarnm))
#            data <- c(data, covar(obj)[covarnm])
#        }
#
#        # add newdata
#        data[names(args)] <- lapply(args, iter, it)
#
#        params <- as.vector(obj@params@.Data)
#        names(params) <- dimnames(obj@params)[['params']]
#
#        # get right dimnames
#        if(length(args) > 0)
#          dimnames <- dimnames(args[[1]])
#        else
#          dimnames <- dimnames(slot(obj, fittedSlot)[[index.count]])
#
#                #browser()
#      # check inputs
#        if(it == 1)
#        {
#          res <- propagate(eval(call,envir=c(params, data, dimdat))[[index.count]], iters, fill.iter=FALSE)
#          dimnames(res)[1:5] <- dimnames[1:5]
#        }
#        else
#        {
#          iter(res, it) <- eval(call,envir=c(params, data, dimdat))[[index.count]]
#        }
#      } # end iter count
#      fitted_index[[index.count]] <- res
#    } # end index.count
#  return(fitted_index)
#  }
#)   # }}}
#
#
#


#predict from aspm
setMethod('predict', signature(object='FLsp'),
  function(object, ...)
  {
      #browser()
#stop("in predict for FLaspm")
    args <- list(...)
    if(length(args) > 0 && is.null(names(args)))
      stop('FLQuant or FLCohort inputs must be named to apply formula')
    # call
    call <- as.list(object@model)[[3]]
    fittedSlot <- as.list(object@model)[[2]]


    # check vars in call match input in args
    if(length(args) > 0 & !any(names(args)%in%all.vars(call)))
      warning(paste("Input names do not match those in model formula: '",
        paste(names(args)[!names(args)%in%all.vars(call)], collapse=','), "'", sep=""))

    # create list of input data
    #   get FLQuant/FLCohort slots' names
    datanm <- getSlotNamesClass(object, 'FLArray')
    datanm <- c(datanm, getSlotNamesClass(object, 'FLQuants'))
    datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))

    # add dimnames if used
    dimna <- dimnames(slot(object, datanm[1]))[names(slot(object, datanm[1]))%in%
      all.vars(object@model)]
    # get them in the right shape
    dimdat <- lapply(dimna, function(x)
      {
        out <- slot(object, datanm[1])
        out[] <- as.numeric(x)
        return(out)
      })

    # iterations
    #   from object
    iter <- max(unlist(qapply(object, function(x) dims(x)$iter)))
    #   from extra input
    if(length(args) > 0)
    {
      iterarg <- lapply(args, function(x) {
        itera <- try(dims(x)$iter)
        if(class(iter) =='try-error')
          return(1)
        else
          return(itera)
      })
      iterarg <- max(unlist(iterarg))
    }
    else
      iterarg <- 1
    #   decision
    if (iter == iterarg)
      iters <- iter
    else if(iter > iterarg && iterarg == 1)
      iters <- iter
    else if(iterarg > iter && iter == 1)
      iters <- iterarg
    else
      stop("Iter for object and input arguments do not match")

    fitted_index <- FLQuants()
    for (index.count in 1:length(object@index))
    {

      for (it in 1:iters)
      {
      obj <- iter(object, it)

      #   input data
        data <- list()
        for (i in datanm)
          data[[i]] <- slot(obj, i)

        # add covar if defined and available
        if('covar' %in% slotNames(obj))
        {
          covarnm <- names(obj@covar)
          if(length(covarnm))
            data <- c(data, covar(obj)[covarnm])
        }

        # add newdata
        data[names(args)] <- lapply(args, iter, it)

        params <- as.vector(obj@params@.Data)
        names(params) <- dimnames(obj@params)[['params']]

        # get right dimnames
        if(length(args) > 0)
          dimnames <- dimnames(args[[1]])
        else
          dimnames <- dimnames(slot(obj, fittedSlot)[[index.count]])

                #browser()
      # check inputs
        if(it == 1)
        {
          res <- propagate(eval(call,envir=c(params, data, dimdat))[[index.count]], iters, fill.iter=FALSE)
          dimnames(res)[1:5] <- dimnames[1:5]
        }
        else
        {
          iter(res, it) <- eval(call,envir=c(params, data, dimdat))[[index.count]]
        }
      } # end iter count
      fitted_index[[index.count]] <- res
    } # end index.count
  return(fitted_index)
  }
)   # }}}


# Also overloading dims()
# overload dims() for FLaspm
# original dims in FLComp does not count FLQuants
# We need to here because index slot is FLQuants and can be multi iter

#setMethod("dims", signature(obj="FLsp"),
#    # Returns a list with different parameters
#    function(obj, ...)
#	{
#    res <- callNextMethod()
#    iters_in_index_slot <- max(unlist(lapply(obj@index,function(x)dim(x)[6])))
#    res$iter <- max(res$iter,iters_in_index_slot)
#    return(res)
#	})

