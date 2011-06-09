# To Do
# Multiple indices
# a plot
# iter check
# rewrite dims method
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


if (!isGeneric("biomass"))
    setGeneric("biomass", function(object, ...)
    standardGeneric("biomass"))

setMethod('biomass', signature(object='FLsp'),
  function(object) {
      # Make an FLQuant with right number of iterations
      iters <- dims(object)$iter
      dimnames <- dimnames(object@catch)
      dimnames$iter <- iters
      biomass <- FLQuant(NA,dimnames=dimnames)
      for (i in 1:iters)
	  iter(biomass,i)[] <- evalC(iter(object,i))[["B"]]
      return(biomass)
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
      qhat <- FLQuant(NA,dimnames=list(iter=iters))
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
      sigma2 <- FLQuant(NA,dimnames=list(iter=iters))
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
      ll <- FLQuant(NA,dimnames=list(iter=iters))
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
      dimnames$iter <- iters
      flq <- FLQuant(NA,dimnames=dimnames)
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

# Need to overload this so that FLQuants slots are included
setMethod("dims", signature(obj="FLsp"),
    # Returns a list with different parameters
    function(obj, ...)
	{
    res <- callNextMethod()
    iters_in_index_slot <- max(unlist(lapply(obj@index,function(x)dim(x)[6])))
    res$iter <- max(res$iter,iters_in_index_slot)
    return(res)
	})



