
# Population dynamics
dyn.load('sra.dll')

pdyn <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {
  
  #browser()
  ymin     <- dims(catch)$minyear
  ymax     <- dims(catch)$maxyear
  ymin_obj <- dims(index)$minyear
  ymax_obj <- dims(index)$maxyear
  #browser()
  B0     <- as.numeric(B0)
  sigma2 <- as.numeric(sigma2)
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
  out <- .Call('fit',B0,sigma2,catch,index,hh,M,m,sel,wght,amin,amax,ymin,ymax,ymin_obj,ymax_obj)
}

pdyn.index <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {

  FLQuant(pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['Ipred']],dimnames=dimnames(index))
}

logl <- function(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax) {
    
  pdyn(B0,sigma2,catch,index,hh,M,mat,sel,wght,amin,amax)[['nLogLk']]
}
  
fit.sra <- function() {

  # use DEoptim?

}