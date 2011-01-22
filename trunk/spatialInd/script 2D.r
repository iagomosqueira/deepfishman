
# 2 area model with no diffusion and advection coefficient based on difference in realised per capita growth
# rate, assuming viscosity of 1 and spatial unit of 1. see p19 of MacCall 1990

library(FLCore)

# population dynamics
dyn <- function(N,r_star,F,V) {

 # calcuate migration coefficients as difference between
 # realised suitabilities, scaled by the viscosity
 mig<- matrix(ncol=dim(N)[5],nrow=dim(N)[5])
 
 for(i in 1:dim(N)[5]) {
  for(j in 1:dim(N)[5]) {
  
    mig[i,j] <- V[i,j]*(r_star[,,,,j] - r_star[,,,,i])
  
  }
 }

 for(i in 1:dim(N)[5]) {
  
  # growth and death
  N[,,,,i] <- max(N[,,,,i] + r_star[,,,,i]*N[,,,,i] - F[,,,,i]*N[,,,,i],0)
  
  for(j in 1:dim(N)[5]) {
  
    # migration
    N[,,,,i] <- max(N[,,,,i] - mig[i,j]*N[,,,,i],0)
  }
 }
  
  return(N)

}

# per captia growth rate
r_calc <- function(N,r,b) {

  r_star <- r - b*N
  return(r_star)
  
}

# start simulation

nyr <- 100

N      <- FLQuant(dim=c(1,nyr,1,1,4))
F      <- FLQuant(dim=c(1,nyr,1,1,4))
r_star <- FLQuant(dim=c(1,nyr,1,1,4))

# viscosity (calculate geographic distance between squares)
V <- matrix(ncol=dim(N)[5],nrow=dim(N)[5])
V[] <- 1

r <- c(0.1,0.2,0.1,0.1)
b <- 0.01

N[,1] <- r/b
r_star[,1] <- r_calc(N[,1],r,b)

F[] <- 0
F[,10:12,,,2] <- 0.6

for(t in 2:nyr) {

  N[,t] <- dyn(N[,t-1],r_star[,t-1],F[,t-1],V)
  
  r_star[,t] <- r_calc(N[,t],r,b)
  
}

windows()
plot(N,type='l',ylab="Abundance")

windows()
plot(r_star,type='l',ylab="r*")



