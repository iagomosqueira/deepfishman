
# 2 area model with no diffusion and advection coefficient based on difference in realised per capita growth
# rate, assuming viscosity of 1 and spatial unit of 1. see p19 of MacCall 1990

library(FLCore)

# population dynamics
dyn <- function(N,r_star,F) {

 # calcuate migration coefficients
  mig<- vector('numeric',length=dim(N)[5])
 
 #diff between realised suitabilities, neg. will move out , positive will move in
 
 mig[1]<- (r_star[,,,,2] - r_star[,,,,1]) # this is the spatial gradient in realised suitability between the two areas
 mig[2]<- (r_star[,,,,1] - r_star[,,,,2]) 

  N[,,,,1] <- max(N[,,,,1] + r_star[,,,,1]*N[,,,,1] - F[,,,,1]*N[,,,,1] - mig[1]*N[,,,,1],0)
  N[,,,,2] <- max(N[,,,,2] + r_star[,,,,2]*N[,,,,2] - F[,,,,2]*N[,,,,2] - mig[2]*N[,,,,2],0)
  
  return(N)

}

# per captia growth rate
r_calc <- function(N,r,b) {

  r_star <- r - b*N
  return(r_star)
  
}

# start simulation

nyr <- 100

N      <- FLQuant(dim=c(1,nyr,1,1,2))
F      <- FLQuant(dim=c(1,nyr,1,1,2))
r_star <- FLQuant(dim=c(1,nyr,1,1,2))

r <- c(0.1,0.2)
b <- 0.01

N[,1] <- r/b
r_star[,1] <- r_calc(N[,1],r,b)

F[] <- 0
F[,,,,1] <- 0.
F[,10:12,,,2] <- 0.6

for(t in 2:nyr) {

 
  N[,t] <- dyn(N[,t-1],r_star[,t-1],F[,t-1])
  
  r_star[,t] <- r_calc(N[,t],r,b)
  
}


windows()
plot(N,type='l',ylab="Abundance")

windows()
plot(r_star,type='l',ylab="r*") # should be = 0 , when an equilibrium is reached?



