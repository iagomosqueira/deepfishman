
library(FLCore)
library(FLaspm)
#data(IOalbacore)        # still doesn't work; something wrong with my computer? get error message when loading package: 'files data/IOalbacore not found' ?????
load('../data/IOalbacore.RData')

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

model(alb) <- aspm.Edwards.C()

# calculate initial starting values
alb@params <- calc.initial(alb)

profile(alb)

# check initial fit
plot(alb)

# now run optimisation        # note that alb@fitted_index does not contain headers for each list item  
alb <- fmle(alb)              
                              


profile(alb,maxsteps=50)
profile(alb,maxsteps=50,fixed=list(sigma2=params(alb)[[2]]))

                              
# check fit
plot(alb)
par(mfrow=c(2,1))
plot(harvest(alb),type='l')
plot(exp.biomass(alb),type="l")

dat <- as.data.frame(FLQuants(h=harvest(alb),b=exp.biomass(alb)))
xyplot(data ~ year | qname,data=dat,type="l",scale=list(relation="free"))

# likelihood profile
profile(alb)

#alb <- fmle(alb,method="SANN")
#alb <- fmle(alb,start=list(B0=params(alb)["B0",],sigma2=params(alb)["sigma2",]),always_eval_initial=FALSE)

# Try AD
albad <- FLaspm(catch=IOalbacore$catch,
              index=IOalbacore$index,
              M=M,hh=hh,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(albad) <- aspm.Edwards.C.AD()
albad <- fmle(albad)
profile(albad)
#albad <- fmle(albad,start=list(B0=120,sigma2=0.12),always_eval_initial=FALSE)
# Do we get a different hessian matrix
alb@hessian
albad@hessian

# Not sure what all this is...
chol(albad@hessian[,,1])
chol(alb@hessian[,,1])
chol2inv(albad@hessian[,,1])
chol2inv(alb@hessian[,,1])
# Look at variance covariance matrix
alb@vcov
albad@vcov
# B0-sigma2 reference has opposite sign


# more than one index

index1 <- IOalbacore$index * rlnorm(dim(IOalbacore$index)[2],0,sqrt(log(1 + 0.1)^2))
index2 <- IOalbacore$index * rlnorm(dim(IOalbacore$index)[2],0,sqrt(log(1 + 0.1)^2))

alb.mi <- FLaspm(catch=IOalbacore$catch,
              index=FLQuants(index1=index1,index2=index2),
              M=M,hh=hh,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(alb.mi) <- aspm.Edwards.C()

alb.mi <- fmle(alb.mi)

profile(alb.mi,maxsteps=50)
profile(alb.mi,maxsteps=50,fixed=list(sigma2=params(alb)[[2]]))

plot(alb.mi)
plot(harvest(alb.mi),type='l')

# likelihood profile

# multiple iterations across hh
hh.it <- runif(10,0.7,0.8)

alb.it <- FLaspm(catch=IOalbacore$catch,
              index=IOalbacore$index,
              M=M,hh=hh.it,sel=as, mat=am,wght=mw,amax=amax, amin=amin)

model(alb.it) <- aspm.Edwards.C()

alb.it <- fmle(alb.it)

boxplot(data~year,as.data.frame(exp.biomass(alb.it)),outline=F)
boxplot(data~year,as.data.frame(harvest(alb.it)),outline=F)
profile(alb.it)

# likelihood profile









