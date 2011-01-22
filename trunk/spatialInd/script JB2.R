rm(list=ls())

# MacCall 1990 example 
# Reproduce Table 2.1 along 1-D spatial transect discretised into 5 grid cells
# local population growth, death and advection (DDHS).

#set-up parameters 

V<-1.25 #viscosity coefficient
b <- 0.4 # density dependence coeff
nyr <- 100 #max num years

#set-up grid of physical suitability
r <- rep(1.0,5)

#set-up matrices for storing variables

r_star <- F <- N  <- matrix(0,5,nyr)


# initial values

N[,1] <- c(0.5,1.5,2.5,0.5,1.0) # initial values of N

r_star[,1] <- r - b*N[,1]  # realised suitability at t=1


#loop through time

for (t in 2:nyr) {
	
#gradient in suitability at t=1, multiplied by viscosity coeff

grad<-(r_star[2:5,t-1]-r_star[1:4,t-1])*V

#fluxes on LHS, leftmost boundary condition = 0 flux

#constraint to limit fluxes below 0.5*N[,t-1] 

lhs<-pmax(c(0.0,grad),-0.5*N[,t-1])


#fluxes on RHS   , rightmost boundary condition = 0 flux

rhs<-pmin(c(grad,0.0),0.5*N[,t-1])

#new abundance


N[,t] <- N[,t-1] + (lhs  - rhs) 
# fluxes from movement only according to Table 2.1, works out to be the same


#simultaneous population dynamics full model, Euler method

#N[,t] <-N[,t-1] + r_star[,t-1]*N[,t-1] - F[,t-1]*N[,t-1] + (lhs - rhs)*N[,t-1]
 
 # equ'n p. 19



r_star[,t] <- r - b*N[,t] # realised suitability at t=1

#population growth, separately?

#N[,t] <-N[,t] + r_star[,t-1]*N[,t-1] - F[,t-1]*N[,t-1] 



}



par(mfrow=c(2,3),mai=rep(0.5,4))
for ( i in 1:5) plot(1:nyr,N[i,],type='l',ylab="Abundance",xlab="Time",ylim=c(0,10),col=4)

# all cells occupied, but different equilibrium abudance than in Table 2.1
# possibly because per capita growth not taken into account in Table? just movement?



par(mfrow=c(2,3),mai=rep(0.5,4))
for ( i in 1:5) plot(1:nyr,r_star[i,],type='l',ylab="r*",xlab="Time",ylim=c(-1,1),col=4)

# equilibrium reached when r*=0, all =0
