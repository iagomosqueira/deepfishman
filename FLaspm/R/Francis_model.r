#********************************************************************************
# Francis model REF
# Using F
# Assumes:
    # mortality (fishing and natural) is constant throughout year, therefore catches in that year are calculated using exploitable biomass at start of year and f
    # survey indices are midyear
    # Recruitment at start of year is based on mature biomass at end of year (check)
#********************************************************************************
# Population dynamics
# Has to estimate f using ratner_search from bexp and catch
aspm.pdyn.Francis <- function(catch,B0,hh,M,mat,sel,wght,amin,amax) {
    #cat("Current B0: ", B0, "\n")
    #browser()
    # Set up stuff
    C <- as.vector(catch)
    nyr <- length(C)
    nag <- amax-amin+1
    yr <- as.numeric(dimnames(catch)[['year']])
    dm <- dimnames(catch)
    n <- array(dim=c(nag,nyr))
    bmat <- vector("numeric",length=nyr)
    bexp <- vector("numeric",length=nyr)
    f <- vector("numeric",length=nyr)
    p <- vector("numeric",length=nag)

    # Equilibrium population
    # number of fish per recruit using n2 = n1 exp(-m); n3 = n2 exp(-m) = n1 exp(-2m)
    p <- exp(-M * (0:(nag-1)))
    # plus group
    p[nag] <- p[nag]/(1-exp(-M))
    # Calc virgin recruits at beginning of year
    R0 <- B0 / sum(p * sel * wght)
    # and the other virgin abundances at the start of the year
    n[,1] <- R0 * p
    # virgin exploitable and mature biomass at start of year
    bexp[1] <- sum(n[,1] * sel * wght)
    bmat[1] <- sum(n[,1] * mat * wght)
    #browser()
    # will need to estimate f from
    # C = B exp(0.5 z) (1 - exp(-z)) * (f / z)
    #f[] <-0
    # from WG report
    #f <- c(0.67, 0.669, 0.305, 0.142, 0.091, 0.088, 0.109, 0.075)
    maxf <- 1e9
    f[1] <- ratner_search(func=fobj,x=c(0,1,maxf),m=M,catch=C[1],biomass=bexp[1])
    # What happens if f maxes out?
    #    if(all.equal(maxf,f[1])==TRUE)
    #    {
    #	warning("Initial B0 is too low")
    #    }

    # set up S-R parameters
    # Here we use BH formulation: R = (a * Bmat) / (b + Bmat)
    # Francis uses this?
    # Recruitment at start of year is calculated using mature biomass at end of previous year
    # hence using (bmat * exp(-M))
    alp <- (4*hh*R0) / (5*hh-1)
    bet <- ((bmat[1] * exp(-M)) * (1 - hh)) / (5*hh - 1)
    #browser()
    # Loop through the years
    #browser()
    for(y in 2:nyr)
    {
	# recruitment using biomass at end of last year
	bmat_end_last <- bmat[y-1]*exp(-M-f[y-1])
	n[1,y] <-  (bmat_end_last * alp) / (bet + bmat_end_last)
	# adults
	n[2:nag,y] <- n[1:(nag-1),y-1] * exp(-M - f[y-1]*sel[1:(nag-1)])
	# plus group
	n[nag,y] <- n[nag,y] + n[nag,y-1]*exp(-M - f[y-1]*sel[nag])
	# But if f in last year maxed out then no suvivors this year
	# Could just set maxf as something huge... 1e9
	#if(all.equal(maxf,f[y-1])==TRUE)
	#{
	#    n[,y] <- n[,y] * (1-sel)
	#}
	# Should quit at this point 
	# biomasses 
	bexp[y] <- sum(n[,y] * sel * wght)
	bmat[y] <- sum(n[,y] * mat * wght)
	# f[y] is estimated
	# not behaving properly if bexp = 0
	# should just be 0, but then stock can grow back...
	f[y] <- ratner_search(func=fobj,x=c(0,1,maxf),m=M,catch=C[y],biomass=bexp[y])
    }
    #browser()
    B0_too_low_warning=F
    #if(maxf %in% f) B0_too_low_warning=T
    # use all.equal
    #if(all.equal(maxf,max(f))==TRUE) B0_too_low_warning=T

    return(list(bexp=FLQuant(bexp,dimnames=dm),
		    # bexp_mid = FLQuant(bexp*exp(-0.5*(M+f)),dimnames=dm),
		    bmat=FLQuant(bmat,dimnames=dm),
		    n = FLQuant(n,dimnames=list(age=amin:amax,year=dm$year)),
		    harvest=FLQuant(f,dimnames=dm,units="f")))
		    #B0_too_low_warning=B0_too_low_warning))

}

# Clean this up
# May not actually need this
aspm.index.Francis <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
{
    #browser()
    # Francis calculates qhat and sigma using the midyear biomass
    #midBexp <- apply(sweep(midN,1,W*sel,"*"),2,sum)
    #n <- sum(!is.na(index)) # number of biomass indices
    #qhat <- sum(index / midBexp,na.rm=T)/n # 0.096, why factor of 1000 out? diff units?
    ## was looking for qhat = 1.03e-4
    #chat2 <- sum((index/(qhat*Bmid)-1)^2,na.rm=T) / (n-2)
    #chat <- sqrt(chat2) # 0.23
# looking for 0.24
    #browser()
    yr <- as.numeric(dimnames(catch)[['year']])
    dm <- dimnames(catch)

    # Strip out the FLQuant to speed it up
    mat <- c(mat)
    sel <- c(sel)
    hh   <- c(hh)
    M   <- c(M)

    index.hat <- FLQuants()
    pdyn <- aspm.pdyn.Francis(catch,B0,hh,M,mat,sel,wght,amin,amax)
    bexp <- pdyn[["bexp"]]
    for (index.count in 1:length(index))
    { 
	# nuisance q
	ind <- as.vector(index[[index.count]])
	iyr <- as.numeric(dimnames(index[[index.count]])[['year']])
	ys <- iyr[1]
	yf <- iyr[length(iyr)]
	y1 <- (ys-yr[1])+1
	y2 <- (yf-yr[1])+1
	# clean all these years up by making the object have same years for each index
	#browser()
	bmid <- bexp*exp(-0.5*(M+pdyn[["harvest"]]))
	    nonnaindexyears <- !is.na(index[[index.count]])
	    n <- dim(index[[index.count]][nonnaindexyears])[2]
	    qhat <- apply(index[[index.count]]/bmid,c(1,6),sum,na.rm=T) / n
	    #    chat2 <- apply((index[[index.count]] / sweep(bmid,1,qhat,"*") - 1)^2,c(1,6),sum,na.rm=T) / (n-2)
	#	q <- exp(mean(log(ind[y1:y2]/as.vector(bexp[,y1:y2])),na.rm=T))
	#return(q)
	#cat("q: ", q, "\n")
	#index.hat[[index.count]] <- FLQuant(q*bexp,dimnames=dm)
	index.hat[[index.count]] <- sweep(bmid,1,qhat,"*")
	#cat("index hat: ", c(index.hat[[index.count]]), "\n")
	#cat("bexp: ", c(bexp), "\n")
	# return predicted index
    }
    #return(FLQuant(q*bexp,dimnames=dm))
    return(index.hat)
}

aspm.Francis <- function()
{
    # set the likelihood function
    logl <- function(B0,hh,M,mat,sel,wght,amin,amax,catch,index)
    {
    
    
    mat <- c(mat)
    sel <- c(sel)
    hh   <- c(hh)
    M   <- c(M)
    #browser()
	# Get the FLQuants object with the estimated indices
	# Actually, do we need indexhat?
	# Francis just uses qhat and chat, no need to calc the indexhat
	#indexhat.quants <- aspm.index.Francis(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
	# What is likelihood function?
	# need q and chat calculation
	#log.indexhat.quants <- lapply(indexhat.quants,function(index) window(log(index),start=dims(index)$minyear,end=dims(index)$maxyear))
	pdyn <- aspm.pdyn.Francis(catch,B0,hh,M,mat,sel,wght,amin,amax)
	#browser()
	bexp <- pdyn[["bexp"]]
	bmid <- bexp*exp(-0.5*(M+pdyn[["harvest"]]))
	# if overfished bmid goes to 0 which kills qhat calculation later on.
	# Set to something small
  #bmid[bmid==0] <- 1e-9
  # But this doesn't work for the chat2 calc, because if only value of bmid is 0
  # qhat is still massive, but bmid/qhat is tiny
  # has weird effect that being just under min B0 gives worse logl than when B0 is very much less than minB0
  # This is risky: if ANY bmid == 0, set all to 0
  if(any(bmid==0)) bmid[] <- 1e-9
  # Gives a flat likelihood - impossible to solve over
  total.logl <- 0
	#browser()
	for (index.count in 1:length(index))
	{
	    nonnaindexyears <- !is.na(index[[index.count]])
	    # number of non NA years in index
	    n <- dim(index[[index.count]][nonnaindexyears])[2]
	    qhat <- apply(index[[index.count]]/bmid,c(1,6),sum,na.rm=T) / n
	    chat2 <- apply((index[[index.count]] / sweep(bmid,1,qhat,"*") - 1)^2,c(1,6),sum,na.rm=T) / (n-2)
	    total.logl <- total.logl + (-n*log(sqrt(chat2)) -n*log(qhat) -apply(log(bmid[nonnaindexyears]),c(1,6),sum))
	}
	# check for warning flag of too low B0 and return massive logl
	#if(pdyn[["B0_too_low_warning"]]) total.logl <- 1000 
	# But fucks up
	return(total.logl)
    }
  
    # initial parameter values
    # Want to start away from B0 crash, but real problems here
    initial <- structure(function(catch){
	    return(FLPar(B0=100*max(catch)))
	},
	# lower and upper limits for optim()
	lower=1,
	upper=1e12
    )

    model <- index ~ aspm.index.Francis(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
    return(list(logl=logl,model=model,initial=initial))
} # }}}





#*******************************************************************************
# Functions for estimating f
#*******************************************************************************

# estimate f objective function
logfobj <- function(logf,m,catch_obs,biomass)
{
    f <- exp(logf)
    z <- f + m
    #catch_hat <- biomass*exp(0.5*z)*(1-exp(-z))*(f/z)
    catch_hat <- biomass*(1-exp(-z))*(f/z)
    return((catch_hat-catch_obs)^2)
}

# test
fobj <- function(f,m,catch,biomass)
{
    z <- f + m
    catch_hat <- biomass*(1-exp(-z))*(f/z)
    return((catch_hat-catch)^2)
}


ratner_search <- function(func,x,abstol=1e-9,...)
{
    #browser()
    phi <- (1+sqrt(5))/2
    resphi <- 2 - phi
    # termination criteria
    if (abs(x[1] - x[3]) < abstol)
	return ((x[1] + x[3])/2)

    # calc new f using golden search
    xnew <- x[2] + resphi * (x[3]- x[2])
    #print(xnew)
    #browser()
    # evaluate obj func at fmid and fnew
    #cat("do call xnew ",do.call(func,list(xnew,...)) ,"\n"); 
    #cat("do call x2 ",do.call(func,list(x[2],...)) ,"\n"); 
    if (do.call(func,list(xnew,...)) < do.call(func,list(x[2],...)))
	ratner_search(func=func,x=c(x[2],xnew,x[3]),abstol=abstol,...)
    else
	ratner_search(func=func,x=c(xnew,x[2],x[1]),abstol=abstol,...)
}

#*******************************************************************************
# C bits

aspm.Francis.C <- function()
{
  # set the likelihood function
  # no sigma2 so set to 1 in .Call
  logl <- function(B0,hh,M,mat,sel,wght,amin,amax,catch,index)
  {
    total.logl <- .Call("aspm_ad", catch, index, B0, 1,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 2)[["logl"]]["logl"]
    return(total.logl)
  }

  initial <- structure(function(hh,M,mat,sel,wght,amin,amax,catch,index){
        B0seq <- seq(from = c(catch)[1], to = 100*max(catch),length=50)
        llout <- rep(NA,length(B0seq))
        for (i in 1:length(B0seq))
          llout[i] <- .Call("aspm_ad", catch, index, B0seq[i], 1,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 2)[["logl"]]["logl"]
        B0_max <- B0seq[which.max(llout)]
        cat("Got initial value of B0: ", B0_max, "\n")
        return(FLPar(B0 = B0_max))
  },
    # initial parameter values
#  initial <- structure(function(catch){
#    return(FLPar(B0=100*max(catch)))
#    },
    # lower and upper limits for optim()
    # Could run profile to get the min B
    lower=1,
    upper=Inf
  )

  model <- index ~ aspm.index.Francis.C(catch,index,B0,hh,M,mat,sel,wght,amin,amax)

  return(list(logl=logl,model=model,initial=initial))
} # }}}


aspm.index.Francis.C <- function(catch,index,B0,hh,M,mat,sel,wght,amin,amax)
{
  # sigma2 not needed so set to 1 in .Call
  indexhat_array <- .Call("aspm_ad", catch, index, B0, 1,
                hh, M, mat, sel,
                wght, amin, amax, dim(catch)[2], 2)[["indexhat"]]
  indexhat_flqs <- FLQuants()
  for (i in 1:length(index))
    indexhat_flqs[[i]] <- FLQuant(indexhat_array[i,],dimnames=dimnames(catch))
  return(indexhat_flqs)
}

