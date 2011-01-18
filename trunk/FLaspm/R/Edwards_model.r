
# Population dynamics
# Uses Charlie's formulation that uses h.
aspm.pdyn.Edwards <- function(catch,B0,hh,M,mat,sel,wght,amin,amax) {
    C <- as.vector(catch)
    nyr <- length(C)
    nag <- amax-amin+1
    yr <- as.numeric(dimnames(catch)[['year']])
    iyr <- as.numeric(dimnames(index)[['year']])
    dm <- dimnames(catch)
    #ind <- as.vector(index)
    ys <- iyr[1]
    yf <- iyr[length(iyr)]
    y1 <- (ys-yr[1])+1
    y2 <- (yf-yr[1])+1
    n <- array(dim=c(nag,nyr))
    bmat <- vector("numeric",length=nyr)
    bexp <- vector("numeric",length=nyr)
    h <- vector("numeric",length=nyr)
    p <- vector("numeric",length=nag)

    # set up eqm population
    p[1] <- 1
    for(a in 2:nag)
	p[a] <- p[a-1]*exp(-M)
    p[nag] <- p[nag]/(1-exp(-M))
    rho <- sum(p * mat * wght)
    R0 <- B0 / rho
    n[,1] <- R0 * p
    bmat[1] <- sum(n[,1] * mat * wght)
    bexp[1] <- sum(n[,1] * sel * wght)
    h[1] <- C[1] / bexp[1]
    h[1] <- max(h[1],0)
    h[1] <- min(h[1],0.999)

    # set up S-R parameters
    alp <- (4*hh*R0)/(5*hh-1)
    #bet <- B0*(1-hh)/(5*hh-1)
    # Should be based on virgin mature population, not virgin exploitable
    bet <- bmat[1]*(1-hh)/(5*hh-1)

    for(y in 2:nyr) {

	n[1,y] <- alp * bmat[y-1]/(bet + bmat[y-1])

	# adult dynamics
	for(a in 2:nag)
	    n[a,y] <- n[a-1,y-1]*exp(-M)*(1-sel[a-1]*h[y-1])
	n[nag,y] <- n[nag,y] + n[nag,y-1]*exp(-M)*(1-sel[nag]*h[y-1])
	bexp[y] <- sum(n[,y] * sel * wght)
	h[y] <- C[y] / bexp[y]
	h[y] <- max(h[y],0)
	h[y] <- min(h[y],0.999)
	bexp[y] <- C[y] / h[y]
	bmat[y] <- sum(n[,y] * mat * wght)
    }

    return(list(bexp=FLQuant(bexp,dimnames=dm),
		    bmat=FLQuant(bmat,dimnames=dm),
		    n = FLQuant(n,dimnames=list(age=amin:amax,year=dm$year)),
		    harvest=FLQuant(h,dimnames=dm, units="h") ))
}

aspm.index.Edwards <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
{
    #browser()
    yr <- as.numeric(dimnames(catch)[['year']])
    dm <- dimnames(catch)


    # Strip out the FLQuant to speed it up
    mat <- c(mat)
    sel <- c(sel)
    hh   <- c(hh)
    M   <- c(M)

    index.hat <- FLQuants()

    for (index.count in 1:length(index))
    { 
	pdyn <- aspm.pdyn.Edwards(catch,B0,hh,M,mat,sel,wght,amin,amax)
	bexp <- pdyn[["bexp"]]
	# nuisance q
	ind <- as.vector(index[[index.count]])
	iyr <- as.numeric(dimnames(index[[index.count]])[['year']])
	ys <- iyr[1]
	yf <- iyr[length(iyr)]
	y1 <- (ys-yr[1])+1
	y2 <- (yf-yr[1])+1
	q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
	index.hat[[index.count]] <- FLQuant(q*bexp,dimnames=dm)
    }
    return(index.hat)
}

aspm.Edwards <- function()
{
    # set the likelihood function
    logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
    {
	#browser()
	# Get the FLQuants object with the estimated indices
	indexhat.quants <- aspm.index.Edwards(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
	log.indexhat.quants <- lapply(indexhat.quants,function(index) window(log(index),start=dims(index)$minyear,end=dims(index)$maxyear))
	total.logl <- 0
	for (index.count in 1:length(index))
	    total.logl <- total.logl + sum(dnorm(log(index[[index.count]]),log.indexhat.quants[[index.count]], sqrt(sigma2), TRUE),na.rm=T)
	return(total.logl)
    }
  
    # initial parameter values
    initial <- structure(function(catch){
	    return(FLPar(B0=100*max(catch), sigma2=1))
	},
	# lower and upper limits for optim()
	lower=c(1, 1e-8),
	upper=c(Inf, Inf)
    )

    model <- index ~ aspm.index.Edwards(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
    return(list(logl=logl,model=model,initial=initial))
} # }}}

#*******************************************************************************
# C Code
aspm.Edwards.C <- function()
{
  # set the likelihood function
  logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
  {
    total.logl <- .Call("aspm_ad", catch, index, B0, sigma2,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 1)[["logl"]]["logl"]
    return(total.logl)
  }

    # initial parameter values
  initial <- structure(function(hh,M,mat,sel,wght,amin,amax,catch,index){
  #browser()
    # Let's do something more sophisticated to get the start values
    B0seq <- seq(from = c(catch)[1], to = 100*max(catch),length=50)
    s2seq <- exp(seq(from = log(1e-8), to = log(10), length = 50))
    llgrid <- expand.grid(B0=B0seq,s2=s2seq,ll=NA)
    
    for (i in 1:nrow(llgrid))
    {
    
      llgrid[i,"ll"] <- .Call("aspm_ad", catch, index, llgrid[i,"B0"], llgrid[i,"s2"],
                            hh, M, mat, sel,
                            wght, amin, amax, dim(catch)[2], 1)[["logl"]]["logl"]
    }
    B0_max <- llgrid[which.max(llgrid[,"ll"]),"B0"]
    s2_max <- llgrid[which.max(llgrid[,"ll"]),"s2"]
    
    #browser()
    cat("Got initial guess\n")
    cat("Initial B0: ", B0_max, "\n")
    cat("Initial sigma2: ", s2_max, "\n")
    return(FLPar(B0=B0_max, sigma2 = s2_max))
    },

  #initial <- structure(function(catch){
    #return(FLPar(B0=100*max(catch), sigma2=1))
    #},
    # lower and upper limits for optim()
    lower=c(1, 1e-8),
    upper=c(Inf, Inf)
  )

  model <- index ~ aspm.index.Edwards.C(catch,index,B0,hh,M,mat,sel,wght,amin,amax)

  return(list(logl=logl,model=model,initial=initial))
} # }}}

aspm.index.Edwards.C <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
{
  # sigma2 not needed so set to 1 in .Call
  indexhat_array <- .Call("aspm_ad", catch, index, B0, 1,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 1)[["indexhat"]]
  indexhat_flqs <- FLQuants()
  for (i in 1:length(index))
    indexhat_flqs[[i]] <- FLQuant(indexhat_array[i,],dimnames=dimnames(catch))
  return(indexhat_flqs)
}

#*******************************************************************************
# model function with gradients from AD
aspm.Edwards.C.AD <- function()
{
  # set the likelihood function
  logl <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
  {
    total.logl <- .Call("aspm_ad", catch, index, B0, sigma2,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 1)[["logl"]]["logl"]
    return(total.logl)
  }

    # initial parameter values
  initial <- structure(function(hh,M,mat,sel,wght,amin,amax,catch,index){
    # Let's do something more sophisticated to get the start values
    # Get a profile surface and pick the parameter values that give the max LL
    B0seq <- seq(from = c(catch)[1], to = 100*max(catch),length=50)
    s2seq <- exp(seq(from = log(1e-8), to = log(10), length = 50))
    llgrid <- expand.grid(B0=B0seq,s2=s2seq,ll=NA)

    for (i in 1:nrow(llgrid))
    {

      llgrid[i,"ll"] <- .Call("aspm_ad", catch, index, llgrid[i,"B0"], llgrid[i,"s2"],
                            hh, M, mat, sel,
                            wght, amin, amax, dim(catch)[2], 1)[["logl"]]["logl"]
    }
    B0_max <- llgrid[which.max(llgrid[,"ll"]),"B0"]
    s2_max <- llgrid[which.max(llgrid[,"ll"]),"s2"]

    #browser()
    cat("Got initial guess\n")
    cat("Initial B0: ", B0_max, "\n")
    cat("Initial sigma2: ", s2_max, "\n")
    return(FLPar(B0=B0_max, sigma2 = s2_max))
    },

  #initial <- structure(function(catch){
    #return(FLPar(B0=100*max(catch), sigma2=1))
    #},
    # lower and upper limits for optim()
    lower=c(1, 1e-8),
    upper=c(Inf, Inf)
  )

gr <- function(B0,sigma2,hh,M,mat,sel,wght,amin,amax,catch,index)
  {
    grads <- -1 * .Call("aspm_ad", catch, index, B0, sigma2,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 1)[["logl"]][c("logl_grad_B0","logl_grad_sigma2")]
    return(grads)
  }


  model <- index ~ aspm.index.Edwards.C(catch,index,B0,hh,M,mat,sel,wght,amin,amax)

  return(list(logl=logl,model=model,initial=initial,gr=gr))
} # }}}


