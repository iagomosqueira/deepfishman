
library(FLaspm)
data(NZOR)
catch <- NZOR[["catch"]]
index <- NZOR[["index"]]

# Parameter values
amin    <- 1
amax    <- 70
hh <- 0.95
M <- 0.05
mat_age <- 23
sel_age <- 23
Linf  <- 42.5 #cm
k     <- 0.059
t0    <- -0.346
alpha <- 0.0963 # grams
beta  <- 2.68
mass_at_age <- age_to_weight(amin:amax,Linf,k,t0,alpha,beta)
w <- FLQuant(mass_at_age, dimnames=list(age=amin:amax)) / 1e6

index.cv <- 0.15 # upper estimate 0.19
niters <- 500

or_true <- FLaspm(catch=catch, index=index, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, amax=amax, amin=amin)
model(or_true) <- aspm.Francis.CAD()
or_true <- fmle(or_true,method="BFGS")
B0_true <- params(or_true)['B0',]
Bmid_true <- B0_true * exp(-0.5*M)

Bmidtrial <- seq(from=350000,to=600000,by=25000)
B0trial <- Bmidtrial*exp(0.5*M)
B0store <- array(NA,dim=c(length(B0trial),niters))

#set.seed(0)
#na.years <- dimnames(index)$year[which(is.na(index))]
# one more year
#na.years <- as.character(1978:1982)
# two more
#na.years <- as.character(1978:1981)
# 3 more
#na.years <- as.character(1978:1980)
# 4 more
#na.years <- as.character(1978:1979)
# 5 more
#na.years <- as.character(1978)
  # Fit a Johnson's Su distribution to get the cumulative distribution function
  # Try fitting Johnson distribution with real data


jcpars <- NULL
Bmidplot <- seq(from=300,to=600,length=100)
jp <- array(NA,dim=c(7,length(Bmidplot)))
initialparms <- c(g=-1,delta=1,xi=0,lambda=1)
for (i in 1:7)
{
  load(paste("C:/Projects/Deepfishman/deepfishman/trunk/B0storeAD_plus",i-1,".Rdata",sep=""))
  Bmidstore <- B0store * exp(-0.5*M)
  propBmid <- apply(Bmidstore>=c(Bmid_true),1,sum) / niters
  #plot(Bmidtrial,propBmid,type="l")
  #lines(x=c(Bmid_true,Bmid_true),y=c(0,1),lty=2)
  jcpars <- rbind(jcpars,optim(par=initialparms,fn=Johnsonll,b=Bmidtrial/1000,m=niters,p=propBmid)$par)
  jp[i,] <- JohnsonPDF(jcpars[i,],Bmidplot)
}

#jc <- JohnsonCum(jcpars[i,],Bmidplot)
#plot(Bmidplot,jc,type="l")
#lines(x=Bmidtrial/1000,y=propBmid,lty=2)
# Looks OK.
# So what is PDF?


plot(Bmidplot,jp[7,],type="n",xlab="Initial biomass (000 t)", ylab= "Probability")
for (i in 1:7)
  lines(x=Bmidplot,y=jp[i,],col=i)



# Risk analysis
# Pick two PDFs - original and 1 more year
# Project forward with different TAC schedules
# See probability of collapse
# see how TAC schedule would not sharper cuts to get same prob of collapse for original PDF

iters <- 200
sigmaR <- 1.2
fcollapse <- 1
# Set up future catch schedules as an array
TAC_schedules <- array(NA,dim=c(5,6),dimnames=list(schedule=1:5,year=1990:1995))
TAC_schedules[1,] <- c(28637,25637,22637,19637,16637,13637)
TAC_schedules[2,] <- c(28637,23637,18637,13637,8637,7500)
TAC_schedules[3,] <- c(28637,21637,14637,7637,7500,7500)
TAC_schedules[4,] <- c(28637,19637,10637,7500,7500,7500)
TAC_schedules[5,] <- c(28637,16637,7500,7500,7500,7500)
#prop_collapse <- rep(NA,dim(TAC_schedules)[1])
prop_collapse <- array(NA,dim=c(nrow(jcpars),dim(TAC_schedules)[1]))
expb <- array(NA,dim=c(dim(TAC_schedules)[1],18))
f <-array(NA,dim=c(dim(TAC_schedules)[1],18))

# Deterministic Projection
for (schedule in 1:dim(TAC_schedules)[1])
{
  # Set up new catch
  catch_proj <- window(catch,start=1978,end=1995)
  catch_proj[,ac(1990:1995)] <- TAC_schedules[schedule,] * 1.3
  # stretch index to fit catch - should add check for this in constructor
  index_proj <- window(index,start=1978,end=1995)
  # Set up FLaspm object with the new catch
  or <- FLaspm(catch=catch_proj, index=index_proj, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin)
  model(or) <- aspm.Francis.CAD()
  params(or)["B0",] <- c(B0_true)
  # Project forward
  or.proj <- calc.pop.dyn(or)
  expb[schedule,] <- c(or.proj[["bexp"]])
  f[schedule,] <- c(or.proj[["harvest"]])
}
# Does anyone collapse? Check f > 1
apply(f,1,function(x)any(x>1))
# First two schedules collapse
# Which means that wider the pdf of B0, the more likely they will survive...
# Makes risk analysis a bit weird.


for (i in 1:nrow(jcpars))
{
cat("i: ", i, "\n")
for (schedule in 1:dim(TAC_schedules)[1])
{
  # Set up new catch
  catch_proj <- window(catch,start=1978,end=1995)
  catch_proj[,ac(1990:1995)] <- TAC_schedules[schedule,] * 1.3
  # stretch index to fit catch - should add check for this in constructor
  index_proj <- window(index,start=1978,end=1995)
  # Set up random SRR residuals
  sr_res <- FLQuant(1,dimnames=dimnames(catch_proj))
  sr_res <- propagate(sr_res,iters)
  sr_res[,ac(1990:1995)] <- rlnorm(iters * (1995-1990+1),meanlog=0,sdlog=sigmaR)
  # Set up FLaspm object with the new catch and sr_res
  or <- FLaspm(catch=catch_proj, index=index_proj, M=M,hh=hh,sel=sel_age,
              mat=mat_age, wght=w, fpm=1, amax=amax, amin=amin,sr_res=sr_res)
  model(or) <- aspm.Francis.CAD()
  # Set up B0 values from previously fitted Johnson distribution
  # Remember we scaled B0 by 1000 for the fit
  params(or) <- propagate(params(or),iters)
  params(or)["B0",] <- rJohnson(iters,jcpars[i,])*1000
  # Project forward
  or.proj <- calc.pop.dyn(or)
  # Who has collapsed? (F >= 1)at any time during the projection
  prop_collapse[i,schedule] <- length(which(apply(or.proj[["harvest"]]>=fcollapse,6,sum)>0))/iters
}
}


#save(TAC_schedules,prop_collapse,jcpars,file="C:/Projects/Deepfishman/deepfishman/trunk/prop_collapse.Rdata")

b1 <- rJohnson(10000,jcpars[1,])*1000
b2 <- rJohnson(10000,jcpars[2,])*1000
par(mfrow=c(2,1))
test <- hist(b1)
hist(b2,breaks=test$breaks)
# b2 has tighter pdf - right

