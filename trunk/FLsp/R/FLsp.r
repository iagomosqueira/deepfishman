# To Do
# Multiple indices
# Fitting function
# methods for returning qhat, ll, sigma2 and so on
# a plot
# fix the model in the function - looks horrible
# Fix fitted and residuals
# predict
# fixed

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
  biomass='FLQuant',
  index='FLQuants'
),
	validity=validFLsp
)

setGeneric("FLsp", function(model, ...){
		standardGeneric("FLsp")
})

setMethod('FLsp', signature(model='ANY'),
  function(model, ...)
  {
		#browser()
    res <- FLModel(model, ..., class='FLsp')
    return(res)
  }
)

setMethod('FLsp', signature(model='missing'),
  function(...)
  {
    #browser()
    res <- FLModel(..., class='FLsp')
    model(res) <- sp
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
	model  <- index ~ indexhat(catch,index,r,k)

	return(list(logl=logl, model=model, initial=initial))
}

#*******************************************************************************
# Accessor things
indexhat <- function(catch, index, r, k)
{
	res <- .Call("flspCpp",catch,index[[1]],r,1,k)
	ihat <- FLQuant(res[["qhat"]]*res[["B"]])
  indexhat_flqs <- FLQuants()
#  for (i in 1:length(index))
#    indexhat_flqs[[i]] <- FLQuant(indexhat_array[i,],dimnames=dimnames(catch))
	indexhat_flqs[[1]] <- ihat
  return(indexhat_flqs)
}




