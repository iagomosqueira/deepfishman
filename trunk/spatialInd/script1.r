
library(FLCore)

dyn2 <- function(N,r_star,F) {

  # calcuate migration coefficients
  mig <- vector('numeric',length=dim(N)[5])
  #mig[1] <- as.numeric((r_star[,,,,2] - r_star[,,,,1])*abs(N[,,,,1]-N[,,,,2]))
  #mig[2] <- as.numeric((r_star[,,,,1] - r_star[,,,,2])*abs(N[,,,,1]-N[,,,,2]))

  #if(r_star[,,,,1] > r_star[,,,,2]) {
  #  # movement into 1
  #  mig[1] <- as.numeric((r_star[,,,,1] - r_star[,,,,2])*(N[,,,,1]-N[,,,,2]))
  #  mig[2] <- as.numeric((r_star[,,,,1] - r_star[,,,,2])*(N[,,,,2]-N[,,,,1]))
  #} else if(r_star[,,,,2] > r_star[,,,,1]) {
  #  # movement into 2
  #  mig[1] <- as.numeric((r_star[,,,,1] - r_star[,,,,2])*(N[,,,,2]-N[,,,,1]))
  #  mig[2] <- as.numeric((r_star[,,,,1] - r_star[,,,,2])*(N[,,,,1]-N[,,,,2]))
  #} else {
  #  mig[1] <- 0
  #  mig[2] <- 0
  #}
  
  r_tmp <- r_star[,,,,1] - r_star[,,,,2]
  if(r_tmp > 0) {
    # movement into 1
    mig[1] <- as.numeric(N[,,,,1]*r_tmp - N[,,,,2]*r_tmp)
    mig[2] <- as.numeric(N[,,,,2]*r_tmp - N[,,,,1]*r_tmp)
  } else if(r_tmp < 0) {
    # movement into 2
    mig[1] <- as.numeric(N[,,,,2]*r_tmp - N[,,,,1]*r_tmp)
    mig[2] <- as.numeric(N[,,,,1]*r_tmp - N[,,,,2]*r_tmp)
  } else {
    mig[1] <- 0
    mig[2] <- 0
  }

  N[,,,,1] <- max(N[,,,,1] + r_star[,,,,1]*N[,,,,1] - F[,,,,1]*N[,,,,1] - mig[1],0)
  N[,,,,2] <- max(N[,,,,2] + r_star[,,,,2]*N[,,,,2] - F[,,,,2]*N[,,,,2] - mig[2],0)    # try dynamics first then movement
  
  return(N)

}

dyn <- function(N,r,b,F) {

  r_star <- r_calc(N,r,b)

  N[,,,,1] <- max(N[,,,,1] + r_star[,,,,1]*N[,,,,1] - F[,,,,1]*N[,,,,1],0)
  N[,,,,2] <- max(N[,,,,2] + r_star[,,,,2]*N[,,,,2] - F[,,,,2]*N[,,,,2],0) 

  r_star <- r_calc(N,r,b)    
  r_tmp <- r_star[,,,,1] - r_star[,,,,2]
  
  # calcuate migration coefficients
  mig <- vector('numeric',length=dim(N)[5])

  if(r_tmp > 0) {
    # movement into 1
    mig[1] <- as.numeric(N[,,,,1]*r_tmp - N[,,,,2]*r_tmp)
    mig[2] <- as.numeric(N[,,,,2]*r_tmp - N[,,,,1]*r_tmp)
  } else if(r_tmp < 0) {
    # movement into 2
    mig[1] <- as.numeric(N[,,,,2]*r_tmp - N[,,,,1]*r_tmp)
    mig[2] <- as.numeric(N[,,,,1]*r_tmp - N[,,,,2]*r_tmp)
  } else {
    mig[1] <- 0
    mig[2] <- 0
  }

  N[,,,,1] <- max(N[,,,,1] - mig[1],0)
  N[,,,,2] <- max(N[,,,,2] - mig[2],0) 
  
  return(N)

}


r_calc <- function(N,r,b) {

  r_star <- r - b*N
  return(r_star)
  
}

nyr <- 60

N      <- FLQuant(dim=c(1,nyr,1,1,2))
F      <- FLQuant(dim=c(1,nyr,1,1,2))
#r_star <- FLQuant(dim=c(1,nyr,1,1,2))

r <- c(0.1,0.2)
b <- 0.01

N[,1] <- r/b
#r_star[,1] <- r_calc(N[,1],r,b)

F[] <- 0
F[,,,,1] <- 0.0
F[,10:12,,,2] <- 0.6

for(t in 2:nyr) {

  #r_star[,t-1] <- r_calc(N[,t-1],r,b)
  #N[,t] <- dyn(N[,t-1],r_star[,t-1],F[,t-1])
  N[,t] <- dyn(N[,t-1],r,b,F[,t-1])
  
}

windows()
plot(N,type='l')
#windows()
#plot(r_star,type='l')



