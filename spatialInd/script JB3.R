rm(list=ls())

# MacCall 1990 example 
# local population growth, death and advection (DDHS).

#set-up parameters 

V<-1.25 #viscosity coefficient
b <- 0.1 # density dependence coeff
nyr <- 100 #max num years
deltaT<-1/365 #small time step
maxiter<-365*nyr 
tvec<-deltaT*1:maxiter


#set-up grid of physical suitability
r<-rep(1,5) #equal suitability, shoudl result in equal abundance
#r <- c(0.0,0.0,0.9,0.8,0.0) # should only occupy 2 cells


#set-up matrices for storing variables

r_star <- F <- N  <- matrix(0,5,nyr*365)


# initial values

N[,1] <- c(0.5,1.5,2.5,0.5,1.0) # initial values of N

r_star[,1] <- r - b*N[,1]  # realised suitability at t=1


#loop through time

for (t in 2:maxiter) {
	
#gradient in suitability at t=1, multiplied by viscosity coeff

grad<-(r_star[2:5,t-1]-r_star[1:4,t-1])*V

#fluxes on LHS, leftmost boundary condition = 0 flux

#constraint to limit fluxes below 0.5*N[,t-1]  - needed? shoudl be OK if deltaT is small 

lhs<-c(0.0,grad)


#fluxes on RHS   , rightmost boundary condition = 0 flux

rhs<-c(grad,0.0)

#new abundance


#population dynamics full model, Euler method, see MacCall equ'n p. 19

N[,t] <-N[,t-1] + r_star[,t-1]*N[,t-1]*deltaT - F[,t-1]*N[,t-1]*deltaT + (lhs - rhs)*N[,t-1]*deltaT
 
 

r_star[,t] <- r - b*N[,t] # realised suitability at t=1




}



par(mfrow=c(2,3),mai=rep(0.5,4))
for ( i in 1:5) plot(tvec,N[i,],type='l',ylab="Abundance",xlab="Time",ylim=c(0,max(N[])),col=4)


X11()

par(mfrow=c(2,3),mai=rep(0.5,4))
for ( i in 1:5) plot(tvec,r_star[i,],type='l',ylab="r*",xlab="Time",ylim=c(-1,1),col=4)

# equilibrium reached when r*=0, all =0
