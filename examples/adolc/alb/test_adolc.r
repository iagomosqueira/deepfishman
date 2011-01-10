
source("prelims.r")
source("../../../code/adolc_tapeless/FLaspm_adolc.r")
dyn.load("../../../code/adolc_tapeless/fit_adolc_tapeless.dll")

# TESTING

ll_func <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  -pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['nLogLk']]
}

grad_func <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  -pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['grad']]
}

grad_func_1 <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  -pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['grad']][1]
}


# test grad func

B0 <- 120000/1e3
sigma2 <- 0.1

grad <- grad_func(B0,sigma2,alb.catch,alb.index,hh,M,m,s,mw,amin,amax)
tiny <- 1e-10

ll  <- ll_func(B0,sigma2,alb.catch,alb.index,hh,M,m,s,mw,amin,amax)
ll1 <- ll_func(B0+tiny,sigma2,alb.catch,alb.index,hh,M,m,s,mw,amin,amax)
ll2 <- ll_func(B0,sigma2+tiny,alb.catch,alb.index,hh,M,m,s,mw,amin,amax)

grad.crude <- c((ll1-ll),(ll2-ll))/tiny
grad
grad.crude
grad / grad.crude
# okay

# fit using mle
start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)
fixed <- list(catch=alb.catch,index=alb.index,hh=hh,M=M,mat=m,sel=s,wght=mw,amin=amin,amax=amax)

res.obj <- mle(ll_func,start=start,lower=lower,upper=upper,fixed=fixed,method="L-BFGS-B")
coef(res.obj)[names(start)]

# fit using optim and grad function for B0 only
res.obj <- optim(120,fn = ll_func,                 sigma2 = sigma2,catch = alb.catch,index = alb.index,hh = hh,M = M,mat = m,wght = mw,sel = s,amin = amin,amax = amax,method = "L-BFGS-B",lower = c(80),upper = c(150),hessian = T)
res.obj$par
res.obj <- optim(120,fn = ll_func,gr = grad_func_1,sigma2 = sigma2,catch = alb.catch,index = alb.index,hh = hh,M = M,mat = m,wght = mw,sel = s,amin = amin,amax = amax,method = "L-BFGS-B",lower = c(80),upper = c(150),hessian = T)
res.obj$par

# RUN FLASPM

#   create the FLaspm object
alb <- FLaspm(catch=alb.catch,
  index=alb.index,
  M=M,hh=hh,sel=s, mat=m, wght=mw, fpm=1, amax=amax, amin=amin)

model(alb) <- aspm()

# initial guess
B0 <- 120000/1e3
sigma2 <- 0.1

start <- list(B0 = B0, sigma2 = sigma2)
lower <- rep(1e-9,2)
upper <- rep(1e10,2)

# call fmle
alb.res <- fmle(alb,start=start,lower=lower,upper=upper,seq.iter=FALSE)
params(alb.res)

# now try with grad
alb@gr <- function(x) {
    
  B0 <- x[1]
  sigma2 <- x[2]
  pdyn(B0,sigma2,alb@catch,alb@index,alb@hh,alb@M,alb@mat,alb@sel,alb@wght,alb@amin,alb@amax)[['grad']] * -1
}

alb.res <- fmle(alb,start=start,lower=lower,upper=upper,seq.iter=FALSE)
params(alb.res)

# try again
B0 <- 200000/1e3
sigma2 <- 0.5

start <- list(B0 = B0, sigma2 = sigma2)
alb.res <- fmle(alb,start=start,lower=lower,upper=upper,seq.iter=FALSE)
params(alb.res)

# try again
B0 <- 80000/1e3
sigma2 <- 0.9

start <- list(B0 = B0, sigma2 = sigma2)
alb.res <- fmle(alb,start=start,lower=lower,upper=upper,seq.iter=FALSE)
params(alb.res)

# accessor functions
biomass(alb.res,'Bexp')
harvest.rate(alb.res)
