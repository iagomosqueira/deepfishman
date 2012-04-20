# efficiency of estimators


library(FLAssess)
library(FLH)
library(ggplot2)
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

# run historic dynamics and 
hist_nyr <- 40
fmsy <- c(refpts(gen)[,'harvest'])[4]
maxFfactor <- 2
maxF <- maxFfactor * fmsy
F1 <- maxF*sin(seq(from=0,to=2 * pi/3,length=hist_nyr))[-1]
stk_hist <- as(gen, 'FLStock')
stk_hist <- window(stk_hist,end=hist_nyr)
ctrl_F1 <- fwdControl(data.frame(year=2:hist_nyr, quantity="f",val=F1))
sr_resid <- FLQuant(1,dimnames=list(age=1,year=dimnames(stock.n(stk_hist))$year,iter=dimnames(stock.n(stk_hist))$iter))
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
stk_init  <- stf(stk_hist,proj_nyr)
proj_end  <- proj_strt + proj_nyr
sr_resid <- window(sr_resid,start=1,end=proj_end)
sr_resid[,(proj_strt+1):proj_end] <- 1

# get catch and index
catchability <- 1e-4
catch_hist <- matrix(catch(stk_hist)[,,drop=TRUE],hist_nyr,1)

# get historic biomass according to sra operating model
nit <- 200
out<-pdyn(SSB0,catch_hist,slope,M,mat,sel,wt,amin,amax)
stk <- vector('list',6)
names(stk) <- c('catch','index','H','ssb','bexp','theta')
stk[[1]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk[[2]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk[[3]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk[[4]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk[[5]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk[[6]] <- matrix(NA,hist_nyr+proj_nyr,nit)

for(i in 1:nit) {
  stk[[1]][1:(hist_nyr+1),i] <- out$catch
  stk[[2]][1:(hist_nyr+1),i] <- out$bexp * catchability
  stk[[3]][1:(hist_nyr+1),i] <- out$H
  stk[[4]][1:(hist_nyr+1),i] <- out$ssb
  stk[[5]][1:(hist_nyr+1),i] <- out$bexp
}


# management targets
msy <- msy.sra(SSB0,stk[['catch']][1:hist_nyr,1],stk[['index']][1:hist_nyr,1],slope,M,mat,sel,wt,amin,amax)
CTAR <- msy$MSY
BTAR <- msy$MSY/msy$H
ITAR <- BTAR * catchability

# define 'optimum' catch theta
stk_opt <- stk

hcr_opt <- function(bexp,year) {

  GB <- as.vector(bexp[year-1,]) 
  
  TAC <- (CTAR * GB)/BTAR 
  
  return(TAC)
}

for(y in (proj_strt+1):(proj_end-1)) {
  # control rule
  stk_opt$catch[y,] <- hcr_opt(stk_opt$bexp,y)
  # project
  out <- pdyn(SSB0,stk_opt$catch[1:y,],slope,M,mat,sel,wt,amin,amax)
  stk_opt[[2]][y,]   <- out$bexp[y,] * catchability
  stk_opt[[3]][y,]   <- out$H[y,]
  stk_opt[[4]][y+1,] <- out$ssb[y+1,]
  stk_opt[[5]][y+1,] <- out$bexp[y+1,]
}

dfr_opt <- data.frame(Time=1:proj_end,B=stk_opt$ssb[,1]/SSB0,H=stk_opt$H[,1])
dfr_opt <- dfr_opt[-proj_end,]

print(ggplot(dfr_opt)
  + geom_line(aes(x=Time,y=B),linetype=1) 
  + geom_line(aes(x=Time,y=H),linetype=2)
  + ylab("SSB/K and Harvest") + xlab("Time")
  )

# now work out performance statistic E[(TAC-theta)^2] for HCR
eff_seq <- seq(500,10000,500)#10^(2:6)
ryr_seq <- 1:10

# simulate observation error
load('../dat/cv_pred_func.Rdata')

index_epsilon <- array(dim=c(length(eff_seq),hist_nyr+proj_nyr,nit))

for(e in 1:length(eff_seq)) {

  cv.index <- cv.pred(eff_seq[e])
  
  for(y in 1:(hist_nyr+proj_nyr)) {
    for(i in 1:nit) {
      index_epsilon[e,y,i] <- rlnorm(1,0,sqrt(log(cv.index+1)))
    }}
      
}

# results array; dim=c(eff,hcr,iter)
rarr <- array(dim=c(length(eff_seq),length(ryr_seq),nit))
dimnames(rarr) <- list(effort=eff_seq,ryr=ryr_seq,iter=1:nit)

# HCR I: mean of historic survey index

hcr_emp <- function(index,year,ryr=1) {

  GB <- apply(index[(year-ryr-1):(year-1),],2,mean) 
  
  TAC <- (CTAR * GB)/ITAR 
  
  return(TAC)
}

for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    stk_emp <- stk

    for(y in (proj_strt+1):(proj_end-1)) {
    
      cat(eff_seq[e],ryr_seq[r],y,'\n')
  
      # control rule
      stk_emp$catch[y,] <- hcr_emp(stk_emp$index,y,ryr_seq[r])
      stk_emp$theta[y,] <- hcr_opt(stk_emp$bexp,y)
    
      # project
      out <- pdyn(SSB0,stk_emp$catch[1:y,],slope,M,mat,sel,wt,amin,amax)
      stk_emp[[2]][y,]   <- out$bexp[y,] * catchability * index_epsilon[e,y,]
      stk_emp[[3]][y,]   <- out$H[y,]
      stk_emp[[4]][y+1,] <- out$ssb[y+1,]
      stk_emp[[5]][y+1,] <- out$bexp[y+1,]
    }
  
    rarr[e,r,] <- MSE(stk_emp)
  }  
}


rmat <- apply(rarr,1:2,mean)
windows(width=18)
par(mfrow=c(1,3))
image(y=ryr_seq,x=eff_seq,z=rmat,ylab='Retrospective years of data',xlab='Survey effort',main='Efficiency',cex.lab=1.5,cex.axis=1.5)
contour(y=ryr_seq,x=eff_seq,z=rmat,add=T,labels='')

eff_max <- which.max(apply(rmat,1,mean)) 
ryr_max <- which.max(apply(rmat,2,mean))

abline(v=eff_seq[eff_max],lwd=2,lty=3)
abline(h=ryr_seq[ryr_max],lwd=2,lty=3)

# include marginals
plot(ryr_seq,rmat[eff_max,],ylab='Efficiency',xlab='Retrospecitve years',main='',type='l')
plot(eff_seq,rmat[,ryr_max],ylab='Efficiency',xlab='Survey effort',main='',type='l')
savePlot(file='../res/hcr_emp.pdf',type='pdf')
save(rarr,rmat,hcr_emp,file='../res/hcr_emp.Rdata')


# HCR II: regression of historic survey index

hcr_reg <- function(index,year,ryr=1) {

  tmp <- c(index[(year-ryr-1):(year-1)])
  cf  <- coef(lm(tmp~c(1:(ryr+1))))
  GB  <- sum(cf[1],cf[2] * (ryr+1),na.rm=TRUE) 
  
  TAC <- (CTAR * GB)/ITAR 
  
  return(TAC)
}

for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    for(i in 1:nit) {

      stk_reg <- stk

      for(y in (proj_strt+1):(proj_end-1)) {
    
        # control rule
        stk_reg$catch[y] <- hcr_reg(stk_reg$index,y,ryr_seq[r])
        stk_reg$theta[y] <- hcr_opt(stk_reg$bexp,y)
      
        # project
        out <- pdyn(SSB0,stk_reg$catch[1:y],stk_reg$index[1:y],slope,M,mat,sel,wt,amin,amax)
        stk_reg[[2]][y]   <- out$Bexp[y] * catchability * index_epsilon[e,y,i]
        stk_reg[[3]][y]   <- out$H[y]
        stk_reg[[4]][y+1] <- out$B[y+1]
        stk_reg[[5]][y+1] <- out$Bexp[y+1]
      }
    
      rarr[e,r,i] <- 1/mean(((stk_reg[['theta']][(proj_strt+1):(proj_end-1)] - stk_reg[['catch']][(proj_strt+1):(proj_end-1)])[stk_reg[['H']][(proj_strt+1):(proj_end-1)]<1])^2)
    }
  }  
}

rmat <- apply(rarr,1:2,mean)
windows(width=18)
par(mfrow=c(1,3))
image(y=ryr_seq,x=eff_seq,z=rmat,ylab='Retrospective years of data',xlab='Survey effort',main='Efficiency',cex.lab=1.5,cex.axis=1.5)
contour(y=ryr_seq,x=eff_seq,z=rmat,add=T,labels='')

eff_max <- which.max(apply(rmat,1,mean)) 
ryr_max <- which.max(apply(rmat,2,mean))

abline(v=eff_seq[eff_max],lwd=2,lty=3)
abline(h=ryr_seq[ryr_max],lwd=2,lty=3)

# include marginals
plot(ryr_seq,rmat[eff_max,],ylab='Efficiency',xlab='Retrospecitve years',main='',type='l')
plot(eff_seq,rmat[,ryr_max],ylab='Efficiency',xlab='Survey effort',main='',type='l')
savePlot(file='../res/hcr_reg.pdf',type='pdf')
save(rarr,rmat,hcr_emp,file='../res/hcr_reg.Rdata')


<<label=hcr_sc2,echo=FALSE>>=
hcr <- function(catch,index,year) {

  GB <- numeric(dim(catch)[6])

  for(i in 1:length(GB)) {
    # get estimated exploitable biomass from SRA
    fit.sra <- optim(SSB0,fn = logl,catch=iter(catch,i),index=iter(index,i),hh=slope,M=m(gen),mat=mat(gen),sel=landings.sel(gen),wght=stock.wt(gen),amin=range(gen)['min'],amax=range(gen)['max'],method = "L-BFGS-B",lower = c(800),upper = c(1500),hessian = T)
    out<-pdyn(B0=fit.sra$par,catch=iter(catch,i),index=iter(index,i),hh=slope,M=m(gen),mat=mat(gen),sel=landings.sel(gen),wght=stock.wt(gen),amin=range(gen)['min'],amax=range(gen)['max'])
    GB[i] <- out$Bexp[year-1]
  }
    
  #if(GB<=BLIM) { TAC <- 0
  #} else { if(GB>BLIM & GB<BTAR) { TAC <- (CTAR * (GB - BLIM))/(BTAR - BLIM) 
  #} else { if(GB>=BTAR) { TAC <- (CTAR * GB)/BTAR
  #}}}

  TAC <- (CTAR * GB)/BTAR

  ctrl <- fwdControl(data.frame(year=year,val=TAC[1],quantity="catch"))
  ctrl@trgtArray <- array(dim=c(1,3,length(TAC)),dimnames=list(1,c("min","val","max"),iter=1:length(TAC)))
  ctrl@trgtArray[,"val",] <- TAC 
  
  return(ctrl)
  
}
@

<<echo=FALSE,term=FALSE>>=
index    <- index_init
stk_proj <- stk_init

for(y in (proj_strt+1):proj_end) {
  # control rule
  ctrl <- hcr(catch(stk_proj)[,1:(y-1)],index[,1:(y-1)],y)
  # project
  stk_proj <- fwd(stk_proj,ctrl=ctrl,sr=list(model=sr_model,params=sr_params),sr.residuals=sr_resid)
  # update survey index
  index[,y] <- quantSums(stock.n(stk_proj)[,y] * stock.wt(stk_proj)[,y] * sel) * q.surv
}

index_sc2 <- index
stk_sc2   <- stk_proj
@

<<label=hcr_sc3,echo=FALSE>>=
hcr <- function(catch,index,year) {
  
  GB <- as.numeric(index[,year-1]) 
  
  #if(GB<=ILIM) { TAC <- 0
  #} else { if(GB>ILIM & GB<ITAR) { TAC <- (CTAR * (GB - ILIM))/(ITAR - ILIM) 
  #} else { if(GB>=ITAR) { TAC <- (CTAR * GB)/ITAR
  #}}}

  TAC <- (CTAR * GB)/ITAR

  ctrl <- fwdControl(data.frame(year=year,val=TAC[1],quantity="catch"))
  ctrl@trgtArray <- array(dim=c(1,3,length(TAC)),dimnames=list(1,c("min","val","max"),iter=1:length(TAC)))
  ctrl@trgtArray[,"val",] <- TAC 
  
  return(ctrl)
  
}
@

<<echo=FALSE,term=FALSE>>=
index    <- index_init
stk_proj <- stk_init

for(y in (proj_strt+1):proj_end) {
  # control rule
  ctrl <- hcr(catch(stk_proj)[,1:(y-1)],index[,1:(y-1)],y)
  # project
  stk_proj <- fwd(stk_proj,ctrl=ctrl,sr=list(model=sr_model,params=sr_params),sr.residuals=sr_resid)
  # survey index
  index[,y] <- quantSums(stock.n(stk_proj)[,y] * stock.wt(stk_proj)[,y] * sel) * q.surv
}

index_sc3 <- index
stk_sc3   <- stk_proj
@

\begin{figure}
\centering
<<fig=TRUE,echo=FALSE>>=

dfr_det <- data.frame(Time=1:proj_end, B1=c(quantSums(sweep(stk_sc1@stock.n * stk_sc1@stock.wt,1,catch.sel(gen),"*"))),B2=c(quantSums(sweep(stk_sc2@stock.n * stk_sc2@stock.wt,1,catch.sel(gen),"*"))),B3=c(quantSums(sweep(stk_sc3@stock.n * stk_sc3@stock.wt,1,catch.sel(gen),"*"))))

print(ggplot(dfr_det)
  + geom_line(aes(x=Time,y=B1),linetype=1) 
  + geom_line(aes(x=Time,y=B2),linetype=2)
  + geom_line(aes(x=Time,y=B3),linetype=3)
  + geom_line(aes(x=Time,y=BTAR),data=data.frame(Time=1:proj_end, BTAR=BTAR))
  + geom_line(aes(x=proj_strt,y=c(0,max(dfr_det$B1,dfr_det$B2,dfr_det$B3))))
  + ylab("Exploitable Biomass") + xlab("Time")
  )
@
\caption{Performance of control rules against exploitable biomass target. Vertical line represents start of projection period. 
Horizontal line represents the target biomass $B^{TAR}$.}
\label{fig:hcr_det_biomass}
\end{figure}

\begin{figure}
\centering
<<fig=TRUE,echo=FALSE>>=

dfr_det <- data.frame(dfr_det, C1=c(catch(stk_sc1)),C2=c(catch(stk_sc2)),C3=c(catch(stk_sc3)))

print(ggplot(dfr_det)
  + geom_line(aes(x=Time,y=C1),linetype=1) 
  + geom_line(aes(x=Time,y=C2),linetype=2)
  + geom_line(aes(x=Time,y=C3),linetype=3)
  + geom_line(aes(x=Time,y=CTAR),data=data.frame(Time=1:proj_end, CTAR=CTAR))
  + geom_line(aes(x=proj_strt,y=c(0,max(dfr_det$C1,dfr_det$C2,dfr_det$C3))))
  + ylab("Catch") + xlab("Time")
  )
@
\caption{Performance of control rules against catch target. Vertical line represents start of projection period. 
Horizontal line represents the target catch $C^{TAR}$.}
\label{fig:hcr_det_catch}
\end{figure}

\section{Stochastic results}

Sources of stochasticity were as follows:
\begin{enumerate}
\item recruitment: multiplicative log-normal noise around the predicted recruitment
\item observation: the survey catch rate was subjected to a degree of noise equivalent to that given in Figure~\ref{fig:effCV}
for a specified level of effort
\end{enumerate}


<<echo=FALSE>>=
srr_sd <- 0.3
nit <- 100
eff_seq <- c(100,1000) #seq(100,1000,200)

# historic sr residuals
sr_resid <- FLQuant(rlnorm(maxt*nit,0,srr_sd),dimnames=list(age=1,year=dimnames(stock.n(stk_hist))$year,iter=1:nit))
stk_hist <- propagate(stk_hist,nit)
stk_hist <- fwd(stk_hist,ctrl=ctrl_F1,sr=list(model=sr_model,params=sr_params),sr.residuals=sr_resid)
  
# set up FLStock and get components
stk_init <- stf(stk_hist,proj_nyr)
  
# future sr residuals
sr_resid <- window(sr_resid,start=1,end=proj_end)
sr_resid[,(proj_strt+1):proj_end] <- rlnorm(proj_nyr*nit,0,srr_sd)

# initialise index
index_init <- quantSums(sweep(stock.n(stk_init) * stock.wt(stk_init),1,sel,"*"))

# simulate observation error
load('../dat/cv_pred_func.Rdata')

index_epsilon <- array(dim=c(length(eff_seq),proj_nyr,nit))

for(e in 1:length(eff_seq)) {

  cv.index <- cv.pred(e)
  
  for(y in 1:proj_nyr)
    for(i in 1:nit)
      index_epsilon[e,y,i] <- rlnorm(1,0,sqrt(log(cv.index+1)))
      
}
rm(fit,cv.pred)

# set up FLStocks object
#stk_sclist <- vector('list',3)
#names(stk_sclist) <- c('Sc1','Sc2','Sc3')

# set up results matrix
res_bexp <- data.frame(effort=eff_seq,Sc1=NA,Sc2=NA,Sc3=NA)

@

<<echo=FALSE>>=
<<hcr_sc1>>

index <- index_init
stk_proj <- stk_init

for(e in 1:length(eff_seq)) {

  for(y in (proj_strt+1):proj_end) {
    # control rule
    ctrl <- hcr(catch(stk_proj),index,y)
    # project
    stk_proj <- fwd(stk_proj,ctrl=ctrl,sr = list(model=sr_model,params=sr_params),sr.residuals=sr_resid)
    # exploitable biomass
    index[,y] <- quantSums(stock.n(stk_proj)[,y] * stock.wt(stk_proj)[,y] * sel) 
  }
  
  # difference between results and ideal information case 
  res_bexp[e,'Sc1'] <- dfr_det[proj_end,'B1'] - c(iterMeans(quantSums(sweep(stk_proj@stock.n[,proj_end] * stk_proj@stock.wt[,proj_end],1,catch.sel(gen),"*"))))
}
@

<<echo=FALSE>>=
<<hcr_sc2>>

index <- index_init
stk_proj <- stk_init

for(e in 1:length(eff_seq)) {
   
  for(y in (proj_strt+1):proj_end) {

    # control rule
    ctrl <- hcr(catch(stk_proj)[,1:(y-1)],index[,1:(y-1)],y)
    # project
    stk_proj <- fwd(stk_proj,ctrl=ctrl,sr=list(model=sr_model, params=sr_params), sr.residuals=sr_resid)
    # update survey index
    index[,y] <- quantSums(stock.n(stk_proj)[,y] * stock.wt(stk_proj)[,y] * sel) * q.surv * index_epsilon[e,y-proj_strt,i]
  }
  
  # difference between results and ideal information case 
  res_bexp[e,'Sc2'] <- dfr_det[proj_end,'B2'] - c(iterMeans(quantSums(sweep(stk_proj@stock.n[,proj_end] * stk_proj@stock.wt[,proj_end],1,sel,"*"))))
}

@

<<echo=FALSE>>=
<<hcr_sc3>>

index <- index_init
stk_proj <- stk_init

for(e in 1:length(eff_seq)) {
   
  for(y in (proj_strt+1):proj_end) {

    # control rule
    ctrl <- hcr(catch(stk_proj)[,1:(y-1)],index[,1:(y-1)],y)
    # project
    stk_proj <- fwd(stk_proj,ctrl=ctrl,sr=list(model=sr_model,params=sr_params),sr.residuals=sr_resid)
    # update survey index
    index[,y] <- quantSums(stock.n(stk_proj)[,y] * stock.wt(stk_proj)[,y] * sel) * q.surv * index_epsilon[e,y-proj_strt,]
  }
  
  # difference between results and ideal information case 
  res_bexp[e,'Sc3'] <- dfr_det[proj_end,'B3'] - c(iterMeans(quantSums(sweep(stk_proj@stock.n[,proj_end] * stk_proj@stock.wt[,proj_end],1,sel,"*"))))
}

@

\begin{figure}
\centering
<<fig=TRUE,echo=FALSE>>=
dfr_stoch <- data.frame(Effort=eff_seq/1000, B1=res_bexp[,'Sc1'],B2=res_bexp[,'Sc2'],B3=res_bexp[,'Sc3'])
print(ggplot(dfr_stoch)
  + geom_line(aes(x=Effort,y=B1),linetype=1) 
  + geom_line(aes(x=Effort,y=B2),linetype=2)
  + geom_line(aes(x=Effort,y=B3),linetype=3)
  + ylab("Relative Exploitable Biomass") + xlab("Effort")
  )
@
\caption{Performance of control rules against exploitable biomass target, for different survey effort values.}
\label{fig:hcr_stoch_biomass}
\end{figure}

\begin{figure}
\centering
<<fig=TRUE,echo=FALSE>>=
dfr_stoch <- data.frame(dfr_stoch, B1_norm=res_bexp[,'Sc1']/mean(res_bexp[,'Sc1']),B2_norm=res_bexp[,'Sc2']/mean(res_bexp[,'Sc2']),B3_norm=res_bexp[,'Sc3']/mean(res_bexp[,'Sc3']))
print(ggplot(dfr_stoch)
  + geom_line(aes(x=Effort,y=B1_norm),linetype=1) 
  + geom_line(aes(x=Effort,y=B2_norm),linetype=2)
  + geom_line(aes(x=Effort,y=B3_norm),linetype=3)
  + ylab("Relative Exploitable Biomass") + xlab("Effort")
  )
@
\caption{Renormalised performance of control rules against exploitable biomass target, for different survey effort values.}
\label{fig:hcr_stoch_biomass_norm}
\end{figure}

\begin{table}
\centering
\begin{tabular}{lc}
\hline
Scenario & Regression coefficient\\
\hline
1: Perfect knowledge & \Sexpr{signif(coef(lm(dfr_stoch$B1_norm~dfr_stoch$Effort))[2])}\\
2: Model based HCR   & \Sexpr{signif(coef(lm(dfr_stoch$B2_norm~dfr_stoch$Effort))[2])}\\
3: Empirical HCR     & \Sexpr{signif(coef(lm(dfr_stoch$B3_norm~dfr_stoch$Effort))[2])}\\
\hline
\end{tabular}
\caption{Regression of renormalised performance statistic against survey effort}
\label{tab:res}
\end{table}


<<echo=FALSE,term=FALSE>>=
dev.off()
@

\newpage
\bibliography{lib}

R
FLR
Sweave
Polacheck


\end{document}

