# Method for fitting FLsp objects using DEoptim
# Stolen wholesale from fmle()

setGeneric('fitsp', function(object, ...)
    standardGeneric('fitsp'))

setMethod('fitsp',
	signature(object="FLsp"),
	function(object,...)
	{

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
    datanm <- getSlotNamesClass(object, 'FLArray')
    # Include FLQuants contents too
		flqs <- getSlotNamesClass(object, 'FLQuants')
    for (i in length(flqs))
      datanm <- c(datanm, names(slot(object, flqs[i])))
    datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))
    #   get those in formals of logl
    datanm <- datanm[datanm%in%names(formals(logl))]


	# Check out what fmle() is doing and do it non generically?


	
		# Refangle the logl function (DEoptim passes variables as vector)
		loglfoo <- function(par) {
      pars <- as.list(par)
      names(pars) <- names(start)
      pars[fixnm] <- lapply(fixed, iter, it)
      return(-1*(do.call(logl, args=c(pars, data))))
    }

		deop <- DEoptim(fn =
	
	
	}
