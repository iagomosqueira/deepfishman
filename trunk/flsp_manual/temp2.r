nzrl <- read.csv("../FLsp/data/NZRL.csv")
saa <- read.csv("../FLsp/data/SAA.csv")
nnh <- read.csv("../FLsp/data/NNH.csv")

save(nzrl, file="../FLsp/data/nzrl.RData")
save(saa, file="../FLsp/data/saa.RData")
save(nnh, file="../FLsp/data/nnh.RData")

#*******************************************************************************
# Uncertainty
# Load the library
library(FLsp)
# Load the New Zealand Rock Lobster data set
#data(nzrl)
nzrl <- read.csv("../FLsp/data/NZRL.csv")
# This is a dataframe with year, catch and cpue
# Make FLQuant objects of the catch and cpue series
catch <- FLQuant(nzrl$catch, dimnames=list(year=nzrl$year))
index <- FLQuant(nzrl$cpue, dimnames=list(year=nzrl$year))
# Create the FLsp object
nzrl <- FLsp(catch=catch,index=index)
nzrl <- fitsp(nzrl)

# hessian
nzrl@hessian
#vcov
#vcov.matrix <- solve(res.betIO[['hessian']])*(-1)
hess <- nzrl@hessian[,,1]
#solve(hess)*(-1)
vcov.matrix <- solve(-1 * hess)

# standard errors simply sqrt of diagonal
# elements of variance-covariance matrix
SE <- sqrt(diag(vcov.matrix))
#This matrix assumes the (multivariate) normality of the parameter distribution
# so with this matrix and the MLE (which forms the mean of this distribution)
# we can use multivariate normal simulation methods to simulate a Monte Carlo sample of the key parameters.

## the mean of the distribution (MLE)
mu <- c(nzrl@params['r',],nzrl@params['k',])

## derive the correlation matrix from covariance matrix
cor.matrix <- vcov.matrix
diag(cor.matrix)[] <- 1
cor.matrix[1,2] <- cor.matrix[1,2] / sqrt(prod(diag(vcov.matrix)))
cor.matrix[2,1] <- cor.matrix[1,2]

## the cholesky decomposition of the correlation matrix
Psi <- t(chol(cor.matrix))

# cholesky is like the sqrt of a matrix
# t(chol(cor.matrix)) %*% chol(cor.matrix)
# so Psi is like sqrt of correlation matrix

## Individual SD of parameters (same as SE - is that right?)
param.sd <- sqrt(diag(vcov.matrix))

## number of Monte Carlo samples
nsam <- 1000

## set up theta matrix of samples
theta <- matrix(nrow=nsam,ncol=length(param.sd))
theta[,1] <- rnorm(nsam)
theta[,2] <- rnorm(nsam)

## use apply to "entrain" the correlation into the samples
theta <- t(apply(theta,1,function(x,Psi){x <- Psi%*%x},Psi))

## rescale to the appropriate mean and variance

theta[,1] <- nzrl@params['r',]+theta[,1]*param.sd[1]
theta[,2] <- nzrl@params['k',]+theta[,2]*param.sd[2]
colnames(theta) <- c("r","k")

# check statistical properties
cor(theta)
cor.matrix

## check covariance
cov(theta)
vcov.matrix

# negative rs?

# Why use correlation matrix?
# cholesky of vcov
ch <- t(chol(vcov.matrix))
# uncorrelated samples
nsam <- 1000
## set up theta matrix of samples
theta <- matrix(nrow=nsam,ncol=2)
theta[,1] <- rnorm(nsam)
theta[,2] <- rnorm(nsam)
theta <- t(apply(theta,1,function(x,ch){x <- ch%*%x},ch))
theta[,1] <- nzrl@params['r',]+theta[,1]
theta[,2] <- nzrl@params['k',]+theta[,2]
apply(theta,2,mean)

