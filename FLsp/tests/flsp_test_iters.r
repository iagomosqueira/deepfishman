# Can we fit multiple iters of catches and index?

library(FLsp)
data(nzrl)
head(nzrl)

set.seed(0)

catch <- FLQuant(nzrl$catch, dimnames=list(year=nzrl$year))
index <- FLQuant(nzrl$cpue, dimnames=list(year=nzrl$year))

# one iter
#nzrl <- FLsp(catch=catch,index=index)
#nzrl <- fitsp(nzrl)
#test <- predict(nzrl)


# multi iters
niter <- 2 #4
catch <- propagate(catch,niter)
index <- propagate(index,niter)

# identical iters
#source("C:\\Projects\\Deepfishman\\deepfishman\\trunk\\FLsp\\R\\fitsp.r")
nzrl <- FLsp(catch=catch,index=index)
test <- iter(nzrl,1)

#ind <- lapply(test@index,function(x) iter(x,1))


nzrl <- fitsp(nzrl)

# doesn't crash - but stuffs returns NA
test <- predict(nzrl)
# what is predict actually doing?
# calls ihat function
# iter function not working - returns all iters of index (FLQuants)
obj <- iter(nzrl,1)
obj@index
# must have fixed this for FLaspm

# Is it being called with weird values by solver?
# can't be because one iter is fine

# block out hessian calc using taped call
# still crashes

# seems to perform at least 2 solver iterations
# browse through
# calling predict on the 3rd iter
# commenting out the calls to predict and it works



# Stop
#*******************************************************************************
# Make iters slightly different
catch <- catch * exp(rnorm(prod(dim(catch)),0,0.1))
index <- index * exp(rnorm(prod(dim(index)),0,0.1))

# Create the FLsp object
nzrl <- FLsp(catch=catch,index=index)
test <- fitsp(nzrl)

#nzrl <- FLsp(catch=iter(catch,4),index=iter(index,4))
#test <- fitsp(nzrl)
# Nothing wrong with 4th iter - then why crash?


