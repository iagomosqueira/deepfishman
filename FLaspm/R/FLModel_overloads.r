# It is necessary to overload fmle and predict
# This is because predict tries to dump the result of predict into the
# 'fitted' slot. The 'fitted' slot is of type FLQuant. However, here we
# return FLQuants (length = length(index)). The prediction is for index hat
# so there may be more than one index hat.

# FLModel structure maybe needs a rethink but here we can do a quick 'n' dirty
# hack. Make a new slot of type FLQuants and dump the result of predict in there.

# Need some validity checks
# e.g. that index FLQuants have same dims (or are forced to have same dims) as catch



setMethod('fmle',
  signature(object='FLaspm', start='ANY'),
  function(object, start, method='L-BFGS-B', fixed=list(),
    control=list(trace=1), lower=rep(-Inf, dim(params(object))[1]),
    upper=rep(Inf, dim(params(object))[1]), seq.iter=TRUE, autoParscale=TRUE,
    tiny_number=1e-6, relAutoParscale=TRUE, ...)
  {

#stop("Just checking that we are in fmle for FLasmp")

      #browser()

    # TODO Check with FL
    args <- list(...)
    call <- sys.call(1)
    logl <- object@logl



#browser()

    # get parameter names by matching elements in param slot
    parnm <- names(formals(logl))[names(formals(logl))%in%
      dimnames(object@params)$param]

    # get fixed parameter names
    fixnm <- names(fixed)
    # fixed must match params
    if(any(!fixnm %in% parnm))
      stop("some named arguments in 'fixed' are not arguments to the
        supplied log-likelihood")
    # HACK! clean up fixed list if elements are named vectors
    fixed <- lapply(fixed, function(x){ names(x) <- NULL; x})

    # create list of input data
    #   get FLQuant slots' names
    datanm <- getSlotNamesClass(object, 'FLArray')
    # Include FLQuants too
    datanm <- c(datanm, getSlotNamesClass(object, 'FLQuants'))
    datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))
    #   get those in formals of logl
    datanm <- datanm[datanm%in%names(formals(logl))]

    # limits
    if(method == 'L-BFGS-B')
    {
      if(missing(lower) && !is.null(lower(object)))
        # if is(lower, function)
        lower <- lower(object)[match(parnm, names(fixed), nomatch=0)==0]
      if(missing(upper) && !is.null(upper(object)))
        upper <- upper(object)[match(parnm, names(fixed), nomatch=0)==0]
    }
    else
    {
      lower <- -Inf
      upper <- Inf
    }

    # gr function
    if(!is.null(body(object@gr)))
      gr <- object@gr
    else
      gr <- NULL

    # create logl function
    loglfoo <- function(par) {
      pars <- as.list(par)
      names(pars) <- names(start)
      pars[fixnm] <- lapply(fixed, iter, it)
      return(-1*(do.call(logl, args=c(pars, data))))
    }

    # Hack that gradient function!
    if(is.null(gr))
	grfoo <- NULL
    else
	grfoo <- function(par) {
	  pars <- as.list(par)
	  names(pars) <- names(start)
	  pars[fixnm] <- lapply(fixed, iter, it)
	  return((do.call(gr, args=c(pars, data))))
	}

    # input data
    alldata <- list()
    for (i in datanm)
      alldata[[i]] <- slot(object, i)

    # add dimnames if used
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

    # iterations
    if(seq.iter)
    {
      iter <- dims(object)$iter
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

    # logLik
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

    #browser()

    # vcov
    object@vcov <- array(NA, dim=c(rep(length(parnm)-length(fixed),2), iter),
      dimnames=list(parnm[!parnm%in%names(fixed)],parnm[!parnm%in%names(fixed)],
      iter=1:iter))
    object@hessian <- object@vcov

# Set up fitted_index and residual_index slots
for (index.count in 1:length(object@index))
{
    index_dmns <- dimnames(object@index[[index.count]])
    object@fitted_index[[index.count]] <- propagate(FLQuant(NA,dimnames=index_dmns),iter)
    object@residuals_index[[index.count]] <- propagate(FLQuant(NA,dimnames=index_dmns),iter)
}

    for (it in 1:iter)
    {
      # data
      if(seq.iter)
        data <- lapply(alldata, iter, it)
      else
        data <- alldata

      # add covar if defined and available
      if('covar' %in% slotNames(object))
      {
        covarnm <- names(object@covar)
        covarnm <- covarnm[covarnm%in%names(formals(logl))]
        if(length(covarnm))
          data <- c(data, covar(object)[covarnm])
      }
      # start values
      if(missing(start)) {
        # add call to @initial
        if(is.function(object@initial))
         start <- as(do.call(object@initial, args=data[names(formals(object@initial))]),
           'list')
        else
          start <- formals(logl)[names(formals(logl))%in%parnm]
      }
      else
        # HACK! clean up fixed list if elements are named vectors
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

    #browser()

      # autoParscale
      if(autoParscale && !'parscale' %in% names(control))
      {
        # named vectors for logl plus/minus tiny_number and diff
        diff_logl <- logl_bump1 <- logl_bump2 <- unlist(start)

        # get logLik for start values
        logl_start <- do.call(logl, args=c(start, data, fixed))

        for(j in names(start))
        {
          # bump up & down each param by tiny_number
          bump_params <- start
          bump_params[[j]] <- bump_params[[j]] * (1 + tiny_number)
          logl_bump1[[j]] <- do.call(logl, args=c(data, bump_params, fixed))
          #
          bump_params <- start
          bump_params[[j]] <- bump_params[[j]] * (1 - tiny_number)
          logl_bump2[[j]] <- do.call(logl, args=c(data, bump_params, fixed))
        }
	#          diff_logl <- 1 / (abs(logl_bump1) + abs(logl_bump2)) / (unlist(start) *
	#    2 * tiny_number)

diff_logl <-  abs(1/(((logl_bump1 - logl_bump2) / (2 * unlist(start) * tiny_number))))

        # relative
# This fails if only one parameter = 1
#        if(relAutoParscale)
#          diff_logl <- diff_logl / max(diff_logl)

        control <- c(control, list(parscale=diff_logl))
      }

      # TODO protect environment
      out <- do.call('optim', c(list(par=unlist(start), fn=loglfoo, method=method,
			      #hessian=TRUE, control=control, lower=lower, upper=upper, gr=gr)))
        hessian=TRUE, control=control, lower=lower, upper=upper, gr=grfoo)))

#browser()


      # output
      # place out$par in right iter dim
      iter(object@params[names(out$par),], it) <- out$par
      # fixed
      if(length(fixed) > 0)
        iter(object@params, it)[fixnm,] <- unlist(lapply(fixed, iter, it))
      # TODO make details list of lists if iter > 1?
      object@details <- list(call=call, value=out$value, count=out$counts,
        convergence=out$convergence, message=out$message)
      # vcov & hessian
      coef <- out$par
      object@vcov[,,it] <-
        if (length(coef))
        {
          if(det(out$hessian) != 0)
          {
            tmphess <- try(solve(out$hessian), silent=TRUE)
            if(class(tmphess) =='try-error')
            {
              matrix(numeric(0), length(coef), length(coef), dimnames=list(names(coef),
                names(coef)))
            } else
            tmphess
          } else
            0
        } else
          0
      object@hessian[,,it] <- -out$hessian

      # logLik
      object@logLik[it] <- -out$value
      attr(object@logLik, 'nobs') <- length(data[[1]])

      #browser()
      # fitted & residuals
      for (index.count in 1:length(object@index))
      {
        iter(object@fitted_index[[index.count]],it) <- predict(iter(object, it))[[index.count]]
        iter(object@residuals_index[[index.count]], it) <-
          iter(slot(object,as.list(object@model)[[2]])[[index.count]],it) - iter(object@fitted_index[[index.count]], it)
      }
    }
    return(object)
  }
)   # }}}

# predict   {{{
setMethod('predict', signature(object='FLaspm'),
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

      #          browser()
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
