
# Population dynamics
dyn.load('sra.dll')
dyn.load('sra_iter.dll')
#dyn.unload('sra.dll')

pdyn <- function(B0,catch,hh,M,mat,sel,wght,amin,amax,srr) {
  
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
  srr    <- srr[ymin:ymax,]
  srr    <- as.vector(srr)

  out <- .Call('run',B0,catch,hh,M,mat,sel,wght,amin,amax,ymin,ymax,nit,srr)
  
  return(out)
}

logl <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  ymin     <- 1
  ymax     <- length(catch)

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
  
  ll <- .Call('fit',B0,catch,index,hh,M,mat,sel,wght,amin,amax,ymin,ymax)
  
  return(ll)

}
  
fit.sra <- function(catch,index,hh,M,mat,sel,wght,amin,amax) {
  
  fit <- optim(SSB0,fn = logl,catch=catch,index=index,hh=hh,M=M,mat=mat,sel=sel,wght=wght,amin=amin,amax=amax,method="L-BFGS-B",lower = 500,upper = 2000)
  return(fit$par)
}

ipred.sra <- function(catch,index,hh,M,mat,sel,wght,amin,amax,year) {

  B0 <- fit.sra(catch,index,hh,M,mat,sel,wght,amin,amax)
  ipred <- pdyn(B0,matrix(catch,ncol=1),hh,M,mat,sel,wght,amin,amax,matrix(1,length(catch),1))$bexp * catchability
  
  if(isTRUE(all.equal(B0,500)) || isTRUE(all.equal(B0,SSB0))  || isTRUE(all.equal(B0,2000)) ) { return(catch[year-1] * ITAR/CTAR) #return(stop(simpleError('opt failure'))) #ipred[] <- as.numeric(NA)}
  } else { return(ipred[year])
  }
}  

msy.sra <- function(B0,catch,hh,M,mat,sel,wght,amin,amax) {

  tmp <- pdyn(B0,catch,hh,M,mat,sel,wght,amin,amax,matrix(1,nrow=nrow(catch),ncol=ncol(catch)))
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

# efficiency: 1/MSE
efficiency <- function(stk) {
  
  xx <- stk[['theta']][(proj_strt+1):(proj_end-1),] - stk[['catch']][(proj_strt+1):(proj_end-1),]
  yy <- stk[['bexp']][(proj_strt+1):(proj_end-1),]
  xx[yy<1e-3] <- NA
  apply(xx,2,function(x) 1/mean(x^2,na.rm=T))

}

# data entropy                                                       
H <- function(stk,year,ryr,cv) {

  # residual error
  xx <- exp(stk[['index']][(year-ryr-1):(year-1),] - stk[['bexp']][(year-ryr-1):(year-1),] * catchability)
  yy <- stk[['bexp']][(year-ryr-1):(year-1),]
  xx[yy<1e-3] <- NA
  apply(xx,2,function(x) -sum(dlnorm(x,0,sqrt(log(cv+1)))*log2(dlnorm(x,0,sqrt(log(cv+1)))),na.rm=T))
}

entropy <- function(stk) {

  xx <- stk[['entropy']][(proj_strt+1):(proj_end-1),]
  apply(xx,2,mean,na.rm=T)
}




