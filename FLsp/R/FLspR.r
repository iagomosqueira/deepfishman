# R implementation of the Polcheck SP model

# By+1 = By + g(By) - Cy
# g(B) = (r / p) B (1 - (B/K)^p) - Pella-Tomlinson
# Iy = q By

# Observation error - assuming error is multiplicative and log-normal
# Assume that B0 = K
# L = prod (exp (- vhat^2 / (2 sigmahat^2_v)) / (sqrt(2 pi) sigmahat)
# vhat = log (C/E)_y - log (C/E)hat_y   (where q By = Iy = Cy / Ey)
# sigma^2hat = sum (vhat^2 / n)
# qhat = exp ((1/n) sum (log (Iy / Bhaty)) )

# index and catch same length
#projectB <- function(r,k,catch,p=1)
#{
#	B <- rep(NA,length(catch)+1)
#	B[1] <- k
#	for (i in 2:(length(catch)+1))
#		B[i] <- B[i-1] + (r/p) * B[i-1] * (1 - (B[i-1] / k)^p) - catch[i-1]
#	B[B<0] <- 0
#	return(B)
#}
#
projectB <- function(r,k,catch,p=1)
{
  niters <- length(r)
	B <- array(NA,dim=c(niters,length(catch)+1))
	B[,1] <- k
	for (i in 2:(length(catch)+1))
		B[,i] <- B[,i-1] + (r/p) * B[,i-1] * (1 - (B[,i-1] / k)^p) - catch[i-1]
	B[B<0] <- 0
	return(B)
}



ll_obs <- function(r,k,catch,index,p=1)
{
  #browser()
	n <- length(catch)
	# Get biomass
	Bhat <- projectB(r,k,catch,p)
	# qhat
	qhat <- exp(sum(log(index / Bhat[1:n])) / n)
	indexhat <- qhat * Bhat
	vhat <- log(index / indexhat[1:n])
	sigma2 <- sum(vhat^2 / n)
	l <- prod(exp(-(vhat^2)/(2*sigma2))	/ sqrt(2*pi*sigma2))
	ll <- log(l)
	return(list(B=Bhat,qhat=qhat, sigma2=sigma2, ll=ll))
}

ll_obs_q <- function(r,k,qhat,catch,index,p=1)
{
  #browser()
	n <- length(catch)
	# Get biomass
	Bhat <- projectB(r,k,catch,p)
	# qhat

#	indexhat <- qhat * Bhat
#	vhat <- log(index / indexhat[1:n])
#	sigma2 <- sum(vhat^2 / n)
#	l <- prod(exp(-(vhat^2)/(2*sigma2))	/ sqrt(2*pi*sigma2))
#	ll <- log(l)

# As Punt F Research 1995
CoverE <- qhat * ((Bhat[1:(length(Bhat)-1)] + Bhat[2:(length(Bhat))])/2)
sigma2 <- sum((log(index) - log(CoverE))^2) / n
ll <- -(n*log(sqrt(sigma2))) * (1 / (2*sigma2)) * sum((log(index) - log(CoverE))^2)
	return(list(B=Bhat,qhat=qhat, sigma2=sigma2, ll=ll))
}


ll_obj <- function(logpars,catch,index,mylower,myupper)
{
	r <- exp(logpars[1])
	k <- exp(logpars[2])
	ll <- -1*ll_obs(r,k,catch,index)$ll

	#cat("r: ", r, "\n")
	#cat("k: ", k, "\n\n")
	#cat("ll: ", ll, "\n\n")

	# bounds
	#browser()
	if (!is.nan(r))
		if (r > myupper[1] | r < mylower[1]) ll <- Inf
	if (!is.nan(k))
		if (k > myupper[2] | k < mylower[2]) ll <- Inf
	if (is.nan(ll)) ll <- Inf

	return(ll)
}

ll_obj_q <- function(logpars,catch,index,mylower,myupper)
{
#  browser()
	r <- exp(logpars[1])
	k <- exp(logpars[2])
	qhat <- exp(logpars[3])
	ll <- -1*ll_obs_q(r,k,qhat,catch,index)$ll

	#cat("r: ", r, "\n")
	#cat("k: ", k, "\n\n")
	#cat("ll: ", ll, "\n\n")
	
	# bounds
	#browser()
	if (!is.nan(r))
		if (r > myupper[1] | r < mylower[1]) ll <- Inf
	if (!is.nan(k))
		if (k > myupper[2] | k < mylower[2]) ll <- Inf
	if (is.nan(ll)) ll <- Inf

	return(ll)
}

