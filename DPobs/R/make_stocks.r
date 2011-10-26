# Conditioning the stocks

# Set up the stock using FLH
# Project forward under some fishing scenario
# Save it


library(FLsp)
library(FLH)
library(ggplot2)

#*******************************************************************************
# Make a codoid
Linf <- 120
maxage <- 25
k <- exp(0.5235792+log(Linf)*-0.4540248) # where did these numbers come from?
# selectivity - sort of dome shaped
a1 <- 0.5
sL <- 0.5
sR <- 5
# maturity
mat95 <- 6
# SRR
s <- 0.75
v <- 1e3
# Make the stock
cod1 <- genBRP(age=1:maxage, Linf=Linf, k=k, a1=a1, sL=sL, sR=sR,
                mat95=mat95, s=s, v=v)

# Check it out
ggplot(as.data.frame(catch.sel(cod1))) + geom_line(aes(age,data)) + ylab("Catch selectivity") + xlab("Age")

#*******************************************************************************
# Set up fishing history
maxt <- 10
fmsy <- c(refpts(cod1)[,'harvest'])[4]
maxFfactor <- 2
maxF <- maxFfactor * fmsy

# One way trip
#F1 <- maxF*sin(seq(from=0,to=pi/2,length=maxt))[-1]

# Some contrast
#F1 <- maxF*sin(seq(from=0,to=pi,length=maxt))[-1]

# Constant at a high level of F - if maxF = 3*fmsy, almost dead
F1 <- rep(maxF,maxt-1)

# Add some noise to F to make it interesting
Fnoise_sd <- 0
F1 <- F1 * rlnorm(length(F1),0,Fnoise_sd)

# Plot fishing history
ggplot(data.frame(time=2:maxt, F=F1)) + geom_line(aes(x=time,y=F)) +
	geom_line(aes(x=time,y=FMSY),data=data.frame(time=2:maxt, FMSY=fmsy)) +
	scale_y_continuous(limits = c(0,maxF*1.3))
	
#*******************************************************************************
# Project stock forward

# Noise level of SRR
srr_sd <- 0.0

# Set up fwd()
stk1 <- as(cod1, 'FLStock')
stk1 <- window(stk1,end=maxt)
ctrl_F1 <- fwdControl(data.frame(year=2:maxt, quantity="f",val=F1))
srr_residuals <- FLQuant(rlnorm(dims(stk1)$year * dims(stk1)$iter,0,srr_sd),
									dimnames=list(age=1,year=dimnames(stock.n(stk1))$year,iter=dimnames(stock.n(stk1))$iter))
stk1 <- fwd(stk1, ctrl=ctrl_F1, sr=list(model=model(cod1), params=params(cod1)), sr.residuals=srr_residuals)

print(plot(stk1))

#*******************************************************************************

