# Speed up Bcurrent
library(FLsp)
# Load the New Zealand Rock Lobster data set
data(nzrl)


catch <- FLQuant(nzrl$catch, dimnames=list(year=nzrl$year))
index <- FLQuant(nzrl$cpue, dimnames=list(year=nzrl$year))
# Create the FLsp object
nzrl <- FLsp(catch=catch,index=index)

nzrl <- fitsp(nzrl)


# with the method
system.time(bcm <- bcurrent(nzrl))

#direct
system.time(bcd <- .Call("flspCpp",c(nzrl@catch),c(nzrl@index[[1]]),c(params(nzrl)['r',1]),1,c(params(nzrl)['k',1]))[["B"]][dim(nzrl@catch)[2]])

# The delay is not the C code
# First call is to biomass
system.time(bcm <- biomass(nzrl))
# It's slow

# To get iters - not too bad - and only called once
system.time(dims(nzrl))
# To call the C func
system.time(evalC(iter(nzrl,1))[["B"]])
# slow is it iter or evalC?
system.time(for(i in 1:100) temp <- evalC(nzrl)[["B"]]) # 2.39, now 0.53
system.time(for(i in 1:100) temp <- iter(nzrl,1))
# Iter is slow and evalC uses it again
# Calls iter again

temp <- iter(nzrl,1)


bcurrent(nzrl)
logl(nzrl)

# Get C code to cope with iters? Too hard - iters could be index, catch or params
# Speed up iters method?
# ANd only call it once

# Can we write quick_iter for evalC
# extra argument iter to it
# then .Call has
# easy to pull out 1 iter from that...
object@catch
object@index[[1]]
object@params['r']
1
object@params['k']


