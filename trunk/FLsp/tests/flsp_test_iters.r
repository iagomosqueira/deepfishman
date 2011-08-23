# Can we fit multiple iters of catches and index?

library(FLsp)
data(nzrl)
head(nzrl)

set.seed(0)

catch <- FLQuant(nzrl$catch, dimnames=list(year=nzrl$year))
index <- FLQuant(nzrl$cpue, dimnames=list(year=nzrl$year))

niter <- 4
catch <- propagate(catch,niter)
index <- propagate(index,niter)

# Make iters slightly different
catch <- catch * exp(rnorm(prod(dim(catch)),0,0.1))
index <- index * exp(rnorm(prod(dim(index)),0,0.1))

# Create the FLsp object
nzrl <- FLsp(catch=catch,index=index)
test <- fitsp(nzrl)

#nzrl <- FLsp(catch=iter(catch,4),index=iter(index,4))
#test <- fitsp(nzrl)
# Nothing wrong with 4th iter - then why crash?


