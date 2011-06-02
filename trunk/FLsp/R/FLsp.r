# To Do
# Multiple indices
# Fitting function
# methods for returning qhat, ll, sigma2 and so on
# a plot


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
		.Call("flspCpp",catch,index[[1]],r,1,k)[["ll"]]


  initial <- structure(function(catch)
	{
		# The function to provide initial values, probably never used
    return(FLPar(r=1, k=c(max(catch)/10)))
	},

  # lower and upper limits for optim()
	lower=rep(1e-9, 2),
	upper=rep(Inf, 2))

	model  <- biomass ~ .Call("flspCpp",catch,index[[1]],r,1,k)[["B"]]

	return(list(logl=logl, model=model, initial=initial))
}
