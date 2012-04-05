
# Population dynamics
dyn.load('../cde/sra.dll')
#dyn.unload('../cde/sra.dll')

pdyn <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {
  
  #browser()
  ymin     <- dims(catch)$minyear
  ymax     <- dims(catch)$maxyear
  ymin_obj <- dims(index)$minyear
  ymax_obj <- dims(index)$maxyear
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
  out <- .Call('fit',B0,catch,index,hh,M,mat,sel,wght,amin,amax,ymin,ymax,ymin_obj,ymax_obj)
  
  return(out)
}

pdyn.index <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {

  FLQuant(pdyn(B0,catch,index,hh,M,mat,sel,wght,amin,amax)[['Ipred']],dimnames=dimnames(index))
}

logl <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  pdyn(B0,catch,index,hh,M,mat,sel,wght,amin,amax)[['nLogLk']]
}
  
fit.sra <- function(catch,index,hh,M,mat,sel,wght,amin,amax) {

  fit <- optim(SSB0,fn = logl,catch,index,hh,M,mat,sel,wght,amin,amax,method = "L-BFGS-B",lower = c(800),upper = c(1500),hessian = T)
  return(fit$par)
}  

msy.sra <- function(B0,catch,index,hh,M,mat,sel,wght,amin,amax) {

  tmp <- pdyn(B0,catch,index,hh,M,mat,sel,wght,amin,amax)
  alp <- tmp$srpar['alpha']
  bet <- tmp$srpar['beta']
  
  nag <- amax-amin+1
  
  PF <- numeric(nag)
  PF[1] <- 1
  
  msy.obj <- function(F,alpha,beta) {
  
    for(a in 2:nag)
      PF[a] = PF[a-1]*exp(-M[a-1])*(1-sel[a-1]*F);
    PF[nag-1] = PF[nag-1] + PF[nag-1]*exp(-M[nag-1])*(1-sel[nag-1]*F);
    
    SPRF <- sum(PF * mat * wght)
    
    RF <- alpha - beta/SPRF
    NF <- RF * PF
    
    Y <- sum(NF * wght * sel * F)
    
    return(Y)
  }
  
  opt <- optimise(msy.obj,c(0.1,0.9),alp,bet,maximum=T)

  return(list(F=opt$maximum,MSY=opt$objective))
}
