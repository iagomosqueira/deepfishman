
# Population dynamics
dyn.load('../cde/sra.dll')
dyn.load('../cde/sra_iter.dll')
#dyn.unload('../cde/sra.dll')

pdyn <- function(B0,catch,hh,M,mat,sel,wght,amin,amax) {
  
  ymin     <- 1
  ymax     <- dim(catch)[1]
  nit      <- dim(catch)[2]
  
  B0     <- as.numeric(B0)
  catch  <- as.vector(catch)
  hh     <- as.numeric(hh)
  M      <- as.vector(M)
  mat    <- as.vector(mat)
  sel    <- as.vector(sel)
  wght   <- as.vector(wght)
  amin   <- as.numeric(amin)
  amax   <- as.numeric(amax)

  out <- .Call('run',B0,catch,hh,M,mat,sel,wght,amin,amax,ymin,ymax,nit)
  
  return(out)
}

fit.func <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {
  
  #browser()
  ymin     <- 1
  ymax     <- length(catch)
  #browser()
  B0     <- as.numeric(B0)
  catch  <- as.vector(catch)
  index  <- as.vector(index)
  hh     <- as.numeric(hh)
  M      <- as.vector(M)
  mat    <- as.vector(mat)
  sel    <- as.vector(sel)
  wght   <- as.vector(wght)
  amin   <- as.numeric(amin)
  amax   <- as.numeric(amax)
  #browser()
  out <- .Call('fit',B0,catch,index,hh,M,mat,sel,wght,amin,amax,ymin,ymax)
  
  return(out)
}

logl <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  fit.func(B0,catch,index,hh,M,mat,sel,wght,amin,amax)[['nLogLk']]
}
  
fit.sra <- function(catch,index,hh,M,mat,sel,wght,amin,amax) {

  fit <- optim(SSB0,fn = logl,catch,index,hh,M,mat,sel,wght,amin,amax,method = "L-BFGS-B",lower = c(800),upper = c(1500),hessian = T)
  return(fit$par)
}  

msy.sra <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {

  tmp <- fit.func(B0,catch,index,hh,M,mat,sel,wght,amin,amax)
  alp <- tmp$srpar['alpha']
  bet <- tmp$srpar['beta']
  
  nag <- amax-amin+1
  
  PF <- numeric(nag)
  PF[1] <- 1
  
  msy.obj <- function(H,alpha,beta) {
  
    for(a in 2:nag)
      PF[a] <- PF[a-1]*exp(-M[a-1])*(1-sel[a-1]*H);
    PF[nag-1] <- PF[nag-1] + PF[nag-1]*exp(-M[nag-1])*(1-sel[nag-1]*H);
    
    SPRF <- sum(PF * mat * wght)
    
    RF <- alpha - beta/SPRF
    NF <- RF * PF
    
    Y <- sum(NF * wght * sel * H * exp(-M[a]/2))
    
    return(Y)
  }
  
  opt <- optimise(msy.obj,c(0.1,0.9),alp,bet,maximum=T)

  return(list(H=opt$maximum,MSY=opt$objective))
}

# overload fwd()
#setMethod("fwd",signature(x="FLQuant",y="FLQuant"),
#  function(x,y) {
#  
#  }
#)

MSE <- function(stk) {
  
  tmp <- stk[['theta']][(proj_strt+1):(proj_end-1),] - stk[['catch']][(proj_strt+1):(proj_end-1),]
  apply(tmp,2,function(x) 1/mean(x^2))

}



