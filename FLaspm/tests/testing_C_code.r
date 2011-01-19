#****************************************************************************
# Testing the common FLaspm methods for both methods on the NZ Orange Roughy data
# Checking that C code gives the same pop.dyn results
#****************************************************************************

#****************************************************************************
# Preliminaries
library(FLCore)
library(FLaspm)

#****************************************************************************

# Data and parameter values
data(NZOR)
catch <- NZOR[["catch"]]
index <- NZOR[["index"]]
# Setting up the test data
amin    <- 1
amax    <- 70
# steepness of BH recruitment
hh <- FLQuant(0.95)
# natural mortality ,not age specific
M <- FLQuant(0.05)
# age-length
Linf  <- 42.5 #cm
k     <- 0.059
t0    <- -0.346
# length-weight
alpha <- 0.0963 # grams
beta  <- 2.68
am <- 23
as <- 23

# Set up FLQuants for weights, maturity and selectivity
# Could make a constructor function for this
w <- FLQuant(age_to_weight(amin:amax,Linf,k,t0,alpha,beta), dimnames=list(age=amin:amax)) / 1e6 # convert to tonnes
# age at maturity
# maturity vector
m <- FLQuant(0, dimnames=list(age=amin:amax))
m[(am - amin + 1):length(amin:amax),] <- 1
# age at commercial selectivity (age of recruitment to fishery)
s <- FLQuant(0, dimnames=list(age=amin:amax))
s[(as - amin + 1):length(amin:amax),] <- 1

#****************************************************************************
# Create the FLaspm object
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(ed) <- aspm.Edwards()

fr <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()

#****************************************************************************
# Test Edwards - R vs C
B0 <- 411000* exp(0.5*c(M))
sigma2 <- 0.1 # for Charlie
ed@params['B0',] <- B0
ed@params['sigma2',] <- sigma2


# Test these methods
# pop.dyn
# exp.biomass
# n
# mat.biomass
# harvest (f or h)
ed.pop.dyn <- pop.dyn(ed)

# Test Ed
edC <- .Call("aspm_ad", ed@catch, ed@index, B0, sigma2,
                ed@hh, ed@M, ed@mat, ed@sel,
                ed@wght, ed@amin, ed@amax, dim(ed@catch)[2], 1)

# Bexp
ed.pop.dyn[["bexp"]]
edC[["bexp"]]
all.equal(c(ed.pop.dyn[["bexp"]]),edC[["bexp"]])

#Bmat
ed.pop.dyn[["bmat"]]
edC[["bmat"]]
all.equal(c(ed.pop.dyn[["bmat"]]),edC[["bmat"]])

# h (or f)
ed.pop.dyn[["harvest"]]
edC[["harvest"]] # ?
all.equal(c(ed.pop.dyn[["harvest"]]),edC[["harvest"]])

# n
ed.pop.dyn[["n"]]
edC[["n"]]
all.equal(c(ed.pop.dyn[["n"]]),c(edC[["n"]]))

# qhat
# calced internally for R code atm
exp(mean(log(ed@index[[1]]/as.vector(ed.pop.dyn[["bexp"]])),na.rm=T))
edC[["qhat"]]

# index_hat
aspm.index.Edwards(ed@catch,ed@index,B0,ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax)
aspm.index.Edwards.C(ed@catch,ed@index,B0,ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax)
all.equal(aspm.index.Edwards(ed@catch,ed@index,B0,ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax),aspm.index.Edwards.C(ed@catch,ed@index,B0,ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax))

# logl
ed@logl(B0,sigma2, ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
edC[["logl"]]

#****************************************************************************
# Test Edwards C as a model
# C Code now has a different intitial value given by searching over the profile
# We can do this because it is fast. R code still uses the old way. Could change
# it but it will be slow
# Create the FLaspm object
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(ed) <- aspm.Edwards()

ed.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.C) <- aspm.Edwards.C()

# Need to get pop.dyn working for C code
ed.res <- fmle(ed)
ed.res@params

ed.C.res <- fmle(ed.C)
ed.C.res@params

#****************************************************************************
# Test Francis - R vs C
B0 <- 411000* exp(0.5*c(M))
fr@params['B0',] <- B0
fr.pop.dyn <- pop.dyn(fr)

# Test FrC
frC <- .Call("aspm_ad", fr@catch, fr@index, B0, 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)

# Bexp
fr.pop.dyn[["bexp"]]
frC[["bexp"]]
all.equal(c(fr.pop.dyn[["bexp"]]),frC[["bexp"]])

#Bmat
fr.pop.dyn[["bmat"]]
frC[["bmat"]]
all.equal(c(fr.pop.dyn[["bmat"]]),frC[["bmat"]])

# h (or f)
fr.pop.dyn[["harvest"]]
frC[["harvest"]] # ?
all.equal(c(fr.pop.dyn[["harvest"]]),frC[["harvest"]])

# n
fr.pop.dyn[["n"]]
frC[["n"]]
all.equal(c(fr.pop.dyn[["n"]]),c(frC[["n"]]))

# qhat
# calced internally for R code atm
bmid <- fr.pop.dyn[["bexp"]] * exp(-0.5*(c(fr@M)+fr.pop.dyn[["harvest"]]))
qhat <- apply(fr@index[[1]]/bmid,c(1,6),sum,na.rm=T) / sum(!is.na(fr@index[[1]]))
qhat
frC[["qhat"]]
all.equal(c(qhat),frC[["qhat"]])

# indexhat
# Not calced for logl but used for residuals I think
ih <- aspm.index.Francis(fr@catch,fr@index,B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax)
ihc <- aspm.index.Francis.C(fr@catch,fr@index,B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax)

# why can't I compare the FLQuants?
# units are not the same - so leave it
all.equal(c(ih[[1]]),c(ihc[[1]]))


#logl
logl <- fr@logl(B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)
logl
frC[["logl"]]
all.equal(c(logl),as.numeric(frC[["logl"]][1]))

# Try logl profile
B0seq <- seq(from=1, to = 800000, length=20)
loglr <- rep(NA,length(B0seq))
loglc <- rep(NA,length(B0seq))
for (i in 1:length(B0seq))
{
  loglr[i] <- fr@logl(B0seq[i],fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)
  loglc[i] <- .Call("aspm_ad", fr@catch, fr@index, B0seq[i], 0, fr@hh, fr@M, fr@mat, fr@sel,
       fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)[["logl"]][1]
}

plot(B0seq,loglr,type="l")
lines(B0seq,loglc,col=3)

# low B0 = 1
B0 <- 2e5#1
test <- .Call("aspm_ad", fr@catch, fr@index, B0, 0, fr@hh, fr@M, fr@mat, fr@sel,
       fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)

fr@params['B0',] <- B0
fr.pop.dyn <- pop.dyn(fr)
fr.pop.dyn[["harvest"]]


#****************************************************************************
# Test Francis C as a model
fr <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(fr) <- aspm.Francis()

fr.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(fr.C) <- aspm.Francis.C()

# Need to get pop.dyn working for C code
fr.res <- fmle(fr)
fr.res@params

fr.C.res <- fmle(fr.C)
fr.C.res@params

profile(fr.res,maxsteps=50, range=0.3)
profile(fr.C.res,maxsteps=50, range = 0.1)

#********************************************************************************
# Profile
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed) <- aspm.Edwards()
ed.res <- fmle(ed)
profile(ed.res,maxsteps=30,range=0.3)

ed.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.C) <- aspm.Edwards.C()
ed.C.res <- fmle(ed.C)
profile(ed.C.res,maxsteps=30,range=0.3)
profile(ed.C.res,range=list(B0=seq(from=300000,to=550000,length=30),
                                        sigma2=exp(seq(from=log(1e-7),to=log(0.01),length=200))))


fr <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(fr) <- aspm.Francis()
fr.res <- fmle(fr)
profile(fr.res,maxsteps=30,range=0.2)

fr.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(fr.C) <- aspm.Francis.C()
fr.C.res <- fmle(fr.C)
profile(fr.C.res,maxsteps=30,range=0.3)

#*******************************************************************************
# Edwards AD
ed <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed) <- aspm.Edwards()

B0seq <- seq(from = 300000, to = 500000, length=10)
#s2seq <- seq(from=0.001, to=2, length=10)
s2seq <- exp(seq(from=log(1e-8), to=log(0.1), length=10))
grads <- expand.grid(B0=B0seq, s2 = s2seq, B0grad_approx=NA, B0grad_ad=NA, s2grad_approx=NA,s2grad_ad=NA, logl=NA)

for (i in 1:nrow(grads))
{
  edC <- .Call("aspm_ad", ed@catch, ed@index, grads[i,"B0"], grads[i,"s2"],
                ed@hh, ed@M, ed@mat, ed@sel,
                ed@wght, ed@amin, ed@amax, dim(ed@catch)[2], 1)

  # AD grads
  grads[i,"B0grad_ad"] <- edC[["logl"]][2]
  grads[i,"s2grad_ad"] <- edC[["logl"]][3]
  grads[i,"logl"] <- edC[["logl"]][1]
  # get approx grad
  loglr <- ed@logl(grads[i,"B0"],grads[i,"s2"], ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
  tiny <- 1+1e-9
  loglr_bumpB0 <- ed@logl(grads[i,"B0"]*tiny,grads[i,"s2"], ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
  grads[i,"B0grad_approx"] <- (loglr_bumpB0 - loglr) / (grads[i,"B0"] * (tiny - 1))
  loglr_bumpsigma2 <- ed@logl(grads[i,"B0"],grads[i,"s2"]*tiny, ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
  grads[i,"s2grad_approx"] <- (loglr_bumpsigma2 - loglr) / (grads[i,"s2"] * (tiny - 1))
}
# Looks OK
ll_mat <- matrix(grads[,"logl"],nrow=length(B0seq),ncol=length(s2seq))
image(x=B0seq,y=s2seq,z=-ll_mat)

# Weird gradient on sigma2 at low values (1e-8) but probably OK
#4262100 1e-08

sigma2 <- 0.52631580#0.2#1e-1
B0 <- 15340.0#5e5#4262100
loglr <- ed@logl(B0,sigma2, ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
tiny <- 1+1e-9
loglr_bumpsigma2 <- ed@logl(B0,sigma2*tiny, ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
(loglr_bumpsigma2 - loglr) / (sigma2 * (tiny - 1))





# Add gradient
grad.Edwards <- function(B0, sigma2, hh, M, mat, sel, wght, amin, amax, catch, index)
{
  out <- .Call("aspm_ad", catch, index, B0, sigma2,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 1)
  grads <-   -1*c(out[["logl"]]["logl_grad_B0"], out[["logl"]]["logl_grad_sigma2"])
  return(grads)
}

B0 <- 5e5
sigma2 <- 0.2
grad.Edwards(B0, sigma2, ed.C@hh, ed.C@M, ed.C@mat, ed.C@sel, ed.C@wght, ed.C@amin, ed.C@amax, ed.C@catch, ed.C@index)

ed.CAD <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.CAD) <- aspm.Edwards.C()
ed.CAD@gr <- grad.Edwards
#ed.CAD <- fmle(ed.CAD,control=list(trace=1,parscale=c(1e8,1)))
ed.CAD <- fmle(ed.CAD,autoParscale=FALSE) # doesn't go anywhere, just that initial guesses were decent
ed.CAD <- fmle(ed.CAD)
params(ed.CAD)
profile(ed.CAD,maxsteps=100)
# log the sigma axis - looks good!
profile(ed.CAD,range=list(B0=seq(from=100000, to= 800000, length=50),
                          sigma2 = exp(seq(from=log(1e-4),to=log(10), length=50))), log="y")

# Get just C model up and running
ed.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.C) <- aspm.Edwards.C()
ed.C <- fmle(ed.C)
#ed.C <- fmle(ed.C, control=list(trace=1,parscale=c(1e8,1)))
params(ed.C)
profile(ed.C,maxsteps=100)

# With gradient, doesn't make any difference in final params
# Path looks a little different OK
#*******************************************************************************
# Try with the new AD model
ed.CAD <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.CAD) <- aspm.Edwards.C.AD()
ed.CAD <- fmle(ed.CAD)
params(ed.CAD)
profile(ed.CAD,maxsteps=100)

# vs just C
ed.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.C) <- aspm.Edwards.C()
ed.C <- fmle(ed.C)
params(ed.C)
profile(ed.C,maxsteps=100)
# slightly different path but end result is the same

# Setting good initial values is crucial
# Then the autoParscale becomes appropriate (else the parscale is set using values
# that are far from the logl and therefore may be wrong)
# Using the gradient makes a difference in the journey but may not make a difference
# to the end result (I guess it depends on the data)

#*******************************************************************************
# Francis AD
B0seq <- seq(from = 300000, to = 500000, length=10)
grads <- data.frame(B0=B0seq,B0grad_approx = NA, B0grad_ad = NA, loglR = NA, loglC=NA)

for (i in 1:nrow(grads))
{
  frC <- .Call("aspm_ad", fr@catch, fr@index, grads[i,"B0"], 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)

  # AD grads
  grads[i,"B0grad_ad"] <- frC[["logl"]][2]
  grads[i,"loglC"] <- frC[["logl"]][1]
  # get approx grad
  loglr <- fr@logl(grads[i,"B0"],fr@hh, fr@M, fr@mat, fr@sel, fr@wght, fr@amin, fr@amax, fr@catch, fr@index)
  grads[i,"loglR"] <- loglr
  tiny <- 1+1e-9
  loglr_bumpB0 <- fr@logl(grads[i,"B0"]*tiny,fr@hh, fr@M, fr@mat, fr@sel, fr@wght, fr@amin, fr@amax, fr@catch, fr@index)
  grads[i,"B0grad_approx"] <- (loglr_bumpB0 - loglr) / (grads[i,"B0"] * (tiny - 1))
}

# terrible

B0 <- 433333
#B0 <- 600000
frC <- .Call("aspm_ad", fr@catch, fr@index, B0, 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)

ll_orig <- frC[["logl"]][1]
tiny <- 1+1e-9
ll_bump <- .Call("aspm_ad", fr@catch, fr@index, B0*tiny, 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)[["logl"]][1]
(ll_bump - ll_orig) / (B0 * (tiny - 1))
frC[["logl"]][2]
# Yes!


# Look at qhat - looks fine @ B0 = 300000 to 600000
qhat_orig <- frC[["qhat"]]
tiny <- 1+1e-9
qhat_bump <- .Call("aspm_ad", fr@catch, fr@index, B0*tiny, 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)[["qhat"]]
(qhat_bump - qhat_orig) / (B0 * (tiny - 1))





# Look at f[0] or f[11]
measure <- "bexp"
i <- 12
f0_orig <- c(frC[[measure]])[i]
tiny <- 1+1e-9
f0_bump <- c(.Call("aspm_ad", fr@catch, fr@index, B0*tiny, 0,
                fr@hh, fr@M, fr@mat, fr@sel,
                fr@wght, fr@amin, fr@amax, dim(fr@catch)[2], 2)[[measure]])[i]
(f0_bump - f0_orig) / (B0 * (tiny - 1))

# fs now look OK
# so does bexp

# look at f[end] and bexp (only AD params in qhat calc)


#*******************************************************************************
# Try using optim with Edwards
logl_optim <- function (x, hh, M, mat, sel, wght, amin, amax, catch, index)
{
    B0 = x[1]
    sigma2 = x[2]
    indexhat.quants <- aspm.index.Edwards(catch, index, B0, hh,
        M, mat, sel, wght, amin, amax)
    log.indexhat.quants <- lapply(indexhat.quants, function(index) window(log(index),
        start = dims(index)$minyear, end = dims(index)$maxyear))
    total.logl <- 0
    for (index.count in 1:length(index)) total.logl <- total.logl +
        sum(dnorm(log(index[[index.count]]), log.indexhat.quants[[index.count]],
            sqrt(sigma2), TRUE), na.rm = T)
    return(-1 * total.logl)
}


initial <- c(B0 = 100*max(ed@catch), sigma2 = 1)
# test logl_optim
logl_optim(c(429690,0.02957), ed@hh, ed@M, ed@mat, ed@sel, ed@wght, ed@amin, ed@amax, ed@catch, ed@index)
# looks same as fmle - good, can we solve?
lower <- c(1,1e-8)
upper <- c(Inf,Inf)
test <- optim(par = initial, fn = logl_optim, lower=lower, upper=upper, method="L-BFGS-B",
                control=list(trace=1, parscale=c(1.93e8,3.68e-1)),
                hh = ed@hh, M=ed@M,mat=ed@mat,
              sel=ed@sel, wght=ed@wght, amin=ed@amin, amax=ed@amax, catch=ed@catch, index=ed@index)
test[["par"]]
# as fmle - good

grad <- function(x, hh, M, mat, sel, wght, amin, amax, catch, index)
{
  B0 <- x[1]
  sigma2 <- x[2]
  out <- .Call("aspm_ad", catch, index, B0, sigma2,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 1)
  grads <- as.numeric(out[["logl"]][c("logl_grad_B0","logl_grad_sigma2")])
  #grads <- -1 * out[["logl"]][c("logl_grad_sigma2")]
  #cat("B0 and sigma2: ", c(B0, sigma2), "\n")
  #cat("grads: ", grads, "\n")
  return(-1*grads)
}

#lower <- c(c(catch)[1],1e-8)
test.gr <- optim(par = initial, fn = logl_optim, lower=lower, upper=upper, method="L-BFGS-B",
                control=list(trace=1,parscale=c(1e8,1)),
                gr=grad,
                hh = ed@hh, M=ed@M,mat=ed@mat,
              sel=ed@sel, wght=ed@wght, amin=ed@amin, amax=ed@amax, catch=ed@catch, index=ed@index)


#********************************************************************************

ed.C <- FLaspm(catch=catch, index=FLQuants(index1 = index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
model(ed.C) <- aspm.Edwards.C()
ed.C <- fmle(ed.C)
params(ed.C)
profile(ed.C,range=0.1,maxsteps=50)

edC <- .Call("aspm_ad", ed@catch, ed@index, 15340, 0.5263,
                ed@hh, ed@M, ed@mat, ed@sel,
                ed@wght, ed@amin, ed@amax, dim(ed@catch)[2], 1)


#********************************************************************************
# Stop
#********************************************************************************
# Fitting test
# If B0 is set too low then are serious issues.
# The likelihood for a B0 lower than B0_min_for_survival is unclear. At the moment
# it is completely flat. This causes problems for the solver.
# So initial B0 has to be above this minimum
# Could try setting some minimum B0 but I don't know of a way to quickly get B0_min_for_survival

# Using default initial and bounds
fr.res <- fmle(fr)
fr.res@params

ed.res <- fmle(ed)
ed.res@params


# Or set the initial and boundary values by hand
# Try setting B0 has 50 x max catch
B0 <- 100*max(fr@catch)
frstart <- list(B0 = B0)
frlower <- 1e-9
frupper <- 1e10

fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=FALSE)
fr.res <- fmle(fr,start=frstart, lower=frlower,upper=frupper, control=list(trace=1),autoParscale=TRUE)
fr.res@params

B0 <- 100*max(ed@catch)
sigma2 <- 1
edstart <- list(B0 = B0, sigma2=sigma2)
edlower <- c(1,1e-9)
edupper <- c(Inf,Inf)

# Same as default
ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1))
# This should be the same
ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1,parscale=c(2e8,0.4)))
# just shite
ed.res <- fmle(ed,start=edstart, lower=edlower,upper=edupper, control=list(trace=1,parscale=c(1e5,1)))
ed.res@params

#********************************************************************************
# Post fitting
# Profile
# Method exists?

# Is the fit any good?
profile(fr.res,maxsteps=30,range=0.3)
profile(ed.res,maxsteps=30,range=0.3)
profile(ed.res,maxsteps=30,range=0.2)
# ed looks a bit dodgy - take a closer look

# fitted logl
ed@logl(params(ed.res)["B0"],params(ed.res)["sigma2"],ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax,ed@catch,ed@index)
# better logl? Yes. Problem with fit
ed@logl(400000,0.024,ed@hh,ed@M,ed@mat,ed@sel,ed@wght,ed@amin,ed@amax,ed@catch,ed@index)

#********************************************************************************
# What happens if B0 too low, i.e. impossible catches
# Francis
fr <- FLaspm(catch=catch, index=FLQuants(index), M=M,hh=hh,sel=s, mat=m, wght=w, fpm=1, amax=amax, amin=amin)
# Set the Francis model
model(fr) <- aspm.Francis()
B0 <- 10#300000#20000
fr@params["B0",] <- B0
#test <- aspm.pdyn.Francis(fr@catch,B0,c(fr@hh),c(fr@M),c(fr@mat),c(fr@sel),c(fr@wght),fr@amin,fr@amax)
test <- pop.dyn(fr)
# harvest maxes out
# and abundance collapses

fr@logl(B0,fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index)

# problem in pdyn and logl
# if bmid = 0, then index / bmid = inf
# set to 1? very hacky...
# how to fail gracefully...

# Plot profile to reveal this
B0seq <- seq(from=10, to=800000,length=100)
frll <- rep(NA,length(B0seq))
for (i in 1:length(B0seq))
    frll[i] <- c(fr@logl(B0seq[i],fr@hh,fr@M,fr@mat,fr@sel,fr@wght,fr@amin,fr@amax,fr@catch,fr@index))
plot(B0seq,frll,type="l",xlab = "B0 start", ylab="logl")



#********************************************************************************

#.Call(aspm_ad, CatchSEXP, SEXP indexSEXP, SEXP B0SEXP, SEXP sigma2SEXP,
#                SEXP steepnessSEXP, SEXP MSEXP, SEXP matSEXP, SEXP selSEXP,
#                SEXP wghtSEXP, SEXP aminSEXP, SEXP amaxSEXP, SEXP nyrsSEXP


#********************************************************************************
# Test creator

amin    <- 1
amax    <- 70
# steepness of BH recruitment
hh <- FLQuant(0.95)
# natural mortality ,not age specific
M <- FLQuant(0.05)
# age-length
Linf  <- 42.5 #cm
k     <- 0.059
t0    <- -0.346
# length-weight
a <- 0.0963 # grams
b  <- 2.68
matage <- 23
selage <- 23

# Pass nothing
test <- FLaspm()

# Pass the model only
test <- FLaspm(model=aspm.Francis())
model(test)

# Everything but model
test <- FLaspm(amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh,M=M)
test <- FLaspm(amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh)
test@M

# Everything
test <- FLaspm(model=aspm.Francis(),amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh,M=M)

# Everything and catch and indices
fr <- FLaspm(model=aspm.Francis(),amin=0, amax=100, matage=matage, selage=selage, Linf=Linf,k=k,t0=t0,a=a,b=b,hh=hh,M=M, catch=catch,index=index)
#test@fpm <- 1
res <- fmle(fr)
res@params
profile(res)
profile(res,maxsteps=30)
profile(res,maxsteps=30,range=0.3)


plot(res)
# Ace
