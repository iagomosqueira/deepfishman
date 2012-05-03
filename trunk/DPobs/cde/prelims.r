
library(FLAssess)
library(FLH)
#library(ggplot2)
library(xtable)
library(mvtnorm)
library(DEoptim)

source("sra.r")

# Parameters for FLH
Linf <- 120
maxage <- 25
k <- exp(0.5235792+log(Linf)*-0.4540248) # where did these numbers come from?
# selectivity
a1 <- 0.5
sL <- Inf
sR <- Inf
# maturity
mat95 <- 6
# SRR
slope <- 0.75
SSB0 <- 1e3
# Make the stock
gen <- genBRP(age=1:maxage, Linf=Linf, k=k, a1=a1, sL=sL, sR=sR, mat95=mat95, s=slope, v=SSB0,minfbar=5,maxfbar=5)

# run historic dynamics
hist_nyr <- 40
nit <- 20
cv.sr <- 0.2
fmsy <- c(refpts(gen)[,'harvest'])[4]
maxFfactor <- 2
maxF <- maxFfactor * fmsy
F1 <- maxF*sin(seq(from=0,to=2 * pi/3,length=hist_nyr))[-1]
stk_hist <- as(gen, 'FLStock')
stk_hist <- window(stk_hist,end=hist_nyr)
stk_hist <- propagate(stk_hist,nit)
ctrl_F1 <- fwdControl(data.frame(year=2:hist_nyr, quantity="f",val=F1))
sr_resid <- FLQuant(rlnorm(nit*hist_nyr,0,sqrt(log(cv.sr+1))),dimnames=list(age=1,year=dimnames(stock.n(stk_hist))$year,iter=dimnames(stock.n(stk_hist))$iter))
stk_hist <- fwd(stk_hist, ctrl=ctrl_F1, sr=list(model=model(gen), params=params(gen)), sr.residuals=sr_resid)

# get information
wt <- as.vector(stock.wt(gen))
M <- as.vector(m(gen))
mat <- as.vector(mat(gen))
sel <- as.vector(catch.sel(gen))
amin <- range(gen)['min']
amax <- range(gen)['max']

# Now we set up for management simulation
proj_nyr  <- 60
proj_strt <- dims(stk_hist)$maxyear
#stk_init  <- stf(stk_hist,proj_nyr)
proj_end  <- proj_strt + proj_nyr
sr_resid <- window(sr_resid,start=1,end=proj_end)
sr_resid[,(proj_strt+1):proj_end] <- rlnorm(nit*proj_nyr,0,sqrt(log(cv.sr+1)))

# set catchability
catchability <- 1e-4

# get historic catch and index
catch_hist <- matrix(catch(stk_hist)[,,drop=TRUE],hist_nyr,nit)
index_hist <- matrix(quantSums(sweep(stock.n(stk_hist) * stock.wt(stk_hist),1,sel,"*"))[,,drop=TRUE],hist_nyr,nit) * catchability

# get historic dynamics according to sra operating model
# with no iterations
#out<-pdyn(SSB0,catch_hist,slope,M,mat,sel,wt,amin,amax)
#stk_opt <- vector('list',6)
#names(stk_opt) <- c('catch','index','H','ssb','bexp','theta')
#stk_opt[[1]] <- matrix(NA,hist_nyr+proj_nyr,1)
#stk_opt[[2]] <- matrix(NA,hist_nyr+proj_nyr,1)
#stk_opt[[3]] <- matrix(NA,hist_nyr+proj_nyr,1)
#stk_opt[[4]] <- matrix(NA,hist_nyr+proj_nyr,1)
#stk_opt[[5]] <- matrix(NA,hist_nyr+proj_nyr,1)
#stk_opt[[6]] <- matrix(NA,hist_nyr+proj_nyr,1)
#
#stk_opt[[1]][1:(hist_nyr+1),1] <- out$catch
#stk_opt[[2]][1:(hist_nyr+1),1] <- out$bexp * catchability
#stk_opt[[3]][1:(hist_nyr+1),1] <- out$H
#stk_opt[[4]][1:(hist_nyr+1),1] <- out$ssb
#stk_opt[[5]][1:(hist_nyr+1),1] <- out$bexp

# management targets
msy <- msy.sra(SSB0,catch_hist,slope,M,mat,sel,wt,amin,amax)
CTAR <- msy$MSY
BTAR <- msy$MSY/msy$H
ITAR <- BTAR * catchability

# describe 'optimum' situation
#hcr_opt <- function(bexp,year) {
#
#  GB <- as.vector(bexp[year,])
#
#  TAC <- (CTAR * GB)/BTAR
#
#  return(TAC)
#}
#
#for(y in (proj_strt+1):(proj_end-1)) {
#  # control rule
#  stk_opt$catch[y,] <- hcr_opt(stk_opt$bexp,y)
#  # project
#  out <- pdyn(SSB0,as.matrix(stk_opt$catch[1:y,]),slope,M,mat,sel,wt,amin,amax)
#  stk_opt[[2]][y,]   <- out$bexp[y,] * catchability
#  stk_opt[[3]][y,]   <- out$H[y,]
#  stk_opt[[4]][y+1,] <- out$ssb[y+1,]
#  stk_opt[[5]][y+1,] <- out$bexp[y+1,]
#}
#
#dfr_opt <- data.frame(Time=1:proj_end,B=stk_opt$ssb[,1]/SSB0,H=stk_opt$H[,1])
#dfr_opt <- dfr_opt[-proj_end,]
#
#print(ggplot(dfr_opt)
#  + geom_line(aes(x=Time,y=B),linetype=1)
#  + geom_line(aes(x=Time,y=H),linetype=2)
#  + ylab("SSB/K and Harvest") + xlab("Time")
#  )
#
# PERFORMANCE SIMULATIONS

eff_seq <- seq(1000,10000,1000)#10^(2:6)
ryr_seq <- c(5,10,15,20,25,30)

# simulate observation error
load('../dat/cv_pred_func.Rdata')

index_epsilon <- array(dim=c(length(eff_seq),hist_nyr+proj_nyr,nit))

for(e in 1:length(eff_seq)) {

  cv.index <- cv.pred(eff_seq[e])

  for(y in 1:(hist_nyr+proj_nyr)) {
    for(i in 1:nit) {
      index_epsilon[e,y,i] <- rlnorm(1,0,sqrt(log(cv.index+1)))
    }
  }
}

sr_epsilon <- array(dim=c(hist_nyr+proj_nyr,nit))

for(y in 1:(hist_nyr+proj_nyr)) {
  for(i in 1:nit) {
    sr_epsilon[y,i] <- sr_resid[,y,,,,i]
  }
}

save(sr_epsilon,index_epsilon,eff_seq,ryr_seq,CTAR,BTAR,ITAR,wt,M,mat,sel,amin,amax,catchability,catch_hist,index_hist,slope,SSB0,hist_nyr,proj_strt,proj_end,proj_nyr,nit,file='prelims_short.Rdata')
