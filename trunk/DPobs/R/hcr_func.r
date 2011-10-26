# The HCR function

hcr <- function(B,r,k,Blim)
{
  #browser()
  catch <- rep(NA,length(B))
  Bmsy <- k/2
  Msy <- r*k/4

  # equation for lower bit of HCR
  m <- Msy / (Bmsy - Blim)
  x <- Msy - m * Bmsy
  catch[B > Bmsy] <- (B * (Msy / Bmsy))[B > Bmsy]
  catch[B <= Bmsy] <- (B * m + x)[B <= Bmsy]
  catch[catch<0] <- 0
  return(catch)
}

#Blim <- 300
#Bvec <- seq(from=0,to=1000,length=50)
#r <- rep(0.5,length(Bvec))
#k <- rep(1000,length(Bvec))

#C <- hcr(Bvec,r,k,Blim)
#plot(Bvec,C,type="l")
