
library(FLCore)
library(FLaspm)
#data(IOalbacore)                 # this doesn't work?
load("../data/IOalbacore.RData")


# input required life history parameters

# age for model
amax <- 9
amin <- 1

# steepness
hh <- 0.75

# natural mortality
M <- 0.34

# mean weight (tonnes)
mw <- 0.014/1e3

# age at commercial selectivity
as <- 5

# age at maturity
am <- 5

#   create the FLaspm object
alb <- FLaspm(catch=IOalbacore$catch,
              index=IOalbacore$index,
              M=M,hh=hh,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(alb) <- aspm.Edwards()

# calculate initial starting values
alb@params <- calc.initial(alb)

# check initial fit
plot(alb)

# now run optimisation        # note that alb@fitted_index does not contain headers for each list item  
alb <- fmle(alb)              # can we change the fitted_flags when we optimise?
                              # do we have to find the initial values again?
                              # currently re-calculates initival values for each fit...
                              
# check fit
plot(alb)
plot(harvest(alb),type='l')

# likelihood profile
profile(alb)

# more than one index

index1 <- IOalbacore$index * rlnorm(dim(alb@index[[1]])[2],0,sqrt(log(1 + 0.1)^2))
index2 <- IOalbacore$index * rlnorm(dim(alb@index[[1]])[2],0,sqrt(log(1 + 0.1)^2))

alb.mi <- FLaspm(catch=IOalbacore$catch,
              index=FLQuants(index1=index1,index2=index2),
              M=M,hh=hh,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(alb.mi) <- aspm.Edwards()

alb.mi <- fmle(alb.mi)

plot(alb.mi)
plot(harvest(alb.mi),type='l')

# likelihood profile

# multiple iterations across hh
hh.it <- runif(10,0.7,0.8)

alb.it <- FLaspm(catch=IOalbacore$catch,
              index=IOalbacore$index,
              M=M,hh=hh.it,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(alb.it) <- aspm.Edwards()

alb.it <- fmle(alb.it)

boxplot(data~year,as.data.frame(exp.biomass(alb.it)),outline=F)
boxplot(data~year,as.data.frame(harvest(alb.it)),outline=F)

# likelihood profile









