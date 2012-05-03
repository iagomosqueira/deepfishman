

# efficiency of estimators
library(DEoptim)
load('prelims_short.Rdata')
source('sra.r')
load('../dat/cv_pred_func.Rdata')

# results array; dim=c(eff,ryr,iter)
rarr <- array(dim=c(length(eff_seq),length(ryr_seq),nit))
dimnames(rarr) <- list(effort=eff_seq,ryr=ryr_seq,iter=1:nit)

earr <- array(dim=c(length(eff_seq),length(ryr_seq),nit))
dimnames(earr) <- list(effort=eff_seq,ryr=ryr_seq,iter=1:nit)

# get historic dynamics over iterations
stk_init <- vector('list',5)
names(stk_init) <- c('catch','index','bexp','theta','entropy')
stk_init[[1]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk_init[[2]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk_init[[3]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk_init[[4]] <- matrix(NA,hist_nyr+proj_nyr,nit)
stk_init[[5]] <- matrix(NA,hist_nyr+proj_nyr,nit)

out <- pdyn(SSB0,catch_hist,slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:hist_nyr,])
stk_init$catch[1:(hist_nyr+1),] <- out$catch
stk_init$index[1:(hist_nyr+1),] <- out$bexp * catchability
stk_init$bexp[1:(hist_nyr+1),]  <- out$bexp

# HCR 0: known stock dynamics
hcr_opt <- function(catch,srr,year) {

  GB <- pdyn(SSB0,catch[1:(year-1),],slope,M,mat,sel,wt,amin,amax,srr[1:(year-1),])$bexp[year,] * catchability
  
  TAC <- (CTAR * GB)/ITAR 
  
  return(TAC)
}

stk <- stk_init

for(y in (proj_strt+1):(proj_end-1)) {

  # control rule
  stk$catch[y,] <- hcr_opt(stk$catch,sr_epsilon,y)

  # project
  out <- pdyn(SSB0,stk$catch[1:y,],slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:y,])
  stk$index[y+1,]  <- out$bexp[y+1,] * catchability
  stk$bexp[y+1,]   <- out$bexp[y+1,]
  
  plot(1:(y+1),out$bexp[1:(y+1),1] * catchability,xlim=c(0,100),type='l')
  points(1:(y+1),stk$index[1:(y+1),1])
  points((proj_strt+1):y,stk$catch[(proj_strt+1):y,1] * ITAR/CTAR,pch=19,col=2,type='b')

}

save(stk,hcr_opt,file='../res/hcr_opt.Rdata')


# HCR I: mean of historic survey index
hcr_mav <- function(index,year,ryr=1) {

  GB <- apply(index[(year-ryr-1):(year-1),],2,mean) 
  
  TAC <- (CTAR * GB)/ITAR 
  
  return(TAC)
}

for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    stk <- stk_init
    stk$index[1:(proj_strt+1),] <- stk$index[1:(proj_strt+1),] * index_epsilon[e,1:(proj_strt+1),]
    
    #cat(eff_seq[e],ryr_seq[r],'\n')
   
    for(y in (proj_strt+1):(proj_end-1)) {  
  
      # control rule
      stk$catch[y,]   <- hcr_mav(stk$index,y,ryr_seq[r])
      stk$theta[y,]   <- hcr_opt(stk$catch,sr_epsilon,y)
      stk$entropy[y,] <- H(stk,y,ryr_seq[r],cv.pred(eff_seq[e]))
    
      # project
      out <- pdyn(SSB0,stk$catch[1:y,],slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:y,])
      stk$index[y+1,] <- out$bexp[y+1,] * catchability * index_epsilon[e,y+1,]
      stk$bexp[y+1,]  <- out$bexp[y+1,]
      
      #plot(1:(y+1),out$bexp[1:(y+1),1] * catchability,xlim=c(0,100),type='l')
      #points((y+1-ryr_seq[r]):(y+1),stk_emp$index[(y+1-ryr_seq[r]):(y+1),1])
      
    }
  
    rarr[e,r,] <- efficiency(stk)
    earr[e,r,] <- entropy(stk)
  }  
}

res<-data.frame(effort=rep(eff_seq,length(ryr_seq*nit)),ryr=rep(rep(ryr_seq,each=length(eff_seq)),nit),efficiency=as.vector(rarr),entropy=as.vector(earr))
save(res,earr,rarr,stk,hcr_mav,file='../res/hcr_mav.Rdata')

# HCR II: regression of historic survey index
hcr_reg <- function(index,year,ryr=1) {

  tmp <- index[(year-ryr-1):(year-1),]
  GB  <- apply(tmp,2,function(x) {
                      cf  <- coef(lm(x~c(1:(ryr+1))))
                      sum(cf[1],cf[2] * (ryr+2),na.rm=T) 
                      })
  
  TAC <- (CTAR * GB)/ITAR 
  
  return(TAC)
}

for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    stk <- stk_init
    stk$index[1:(proj_strt+1),] <- stk$index[1:(proj_strt+1),] * index_epsilon[e,1:(proj_strt+1),]
    
    cat(eff_seq[e],ryr_seq[r],'\n')
    
    for(y in (proj_strt+1):(proj_end-1)) {

      # control rule
      stk$catch[y,] <- hcr_reg(stk$index,y,ryr_seq[r])
      stk$theta[y,] <- hcr_opt(stk$catch,sr_epsilon,y)
      stk$entropy[y,] <- H(stk,y,ryr_seq[r],cv.pred(eff_seq[e]))
    
      # project
      out <- pdyn(SSB0,stk$catch[1:y,],slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:y,])
      stk$index[y+1,] <- out$bexp[y+1,] * catchability * index_epsilon[e,y,]
      stk$bexp[y+1,]  <- out$bexp[y+1,]
    }
  
    rarr[e,r,] <- efficiency(stk)
    earr[e,r,] <- entropy(stk)

  }  
}

res<-data.frame(effort=rep(eff_seq,length(ryr_seq*nit)),ryr=rep(rep(ryr_seq,each=length(eff_seq)),nit),efficiency=as.vector(rarr),entropy=as.vector(earr))
save(res,earr,rarr,stk,hcr_reg,file='../res/hcr_reg.Rdata')

# HCR III: smoothed historic index (exponential)
hcr_smt <- function(index,year,ryr=1) {

   xx <- index[(year-ryr-1):(year-1),]
   alp <- 0.5
   gcf <- rev(c(dgeom(0:(dim(xx)[1]-2),alp),(1-alp)^(dim(xx)[1]-1)))
   
   GB  <- apply(xx,2,function(x) sum(gcf * x))
   
   TAC <- (CTAR * GB)/ITAR

   return(TAC)
}

for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    stk <- stk_init
    stk$index[1:(proj_strt+1),] <- stk$index[1:(proj_strt+1),] * index_epsilon[e,1:(proj_strt+1),]
    
    cat(eff_seq[e],ryr_seq[r],'\n')
    
    for(y in (proj_strt+1):(proj_end-1)) {
  
      # control rule
      stk$catch[y,] <- hcr_smt(stk$index,y,ryr_seq[r])
      stk$theta[y,] <- hcr_opt(stk$catch,sr_epsilon,y)
      stk$entropy[y,] <- H(stk,y,ryr_seq[r],cv.pred(eff_seq[e]))
    
      # project
      out <- pdyn(SSB0,stk$catch[1:y,],slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:y,])
      stk$index[y+1,] <- out$bexp[y+1,] * catchability * index_epsilon[e,y,]
      stk$bexp[y+1,]  <- out$bexp[y+1,]
      
      #plot(1:(y+1),out$bexp[1:(y+1),10] * catchability,xlim=c(0,100),type='l')
      #points((y+1-ryr_seq[r]):(y+1),stk$index[(y+1-ryr_seq[r]):(y+1),10])
      #points(proj_strt:y,stk$catch[proj_strt:y,10] * ITAR/CTAR,pch=19,col=2,type='b')
    }
  
    rarr[e,r,] <- efficiency(stk)
    earr[e,r,] <- entropy(stk)

  }  
}

res<-data.frame(effort=rep(eff_seq,length(ryr_seq*nit)),ryr=rep(rep(ryr_seq,each=length(eff_seq)),nit),efficiency=as.vector(rarr),entropy=as.vector(earr))
save(res,earr,rarr,stk,hcr_smt,file='../res/hcr_smt.Rdata')

# HCR IV: smoothed historic index (double exponential)
hcr_dex <- function(index,year,ryr=1) {

   xx <- index[(year-ryr-1):(year-1),]
   alp <- 0.5
   bet <- 0.5
   
   s <- numeric(ryr+1)
   b <- numeric(ryr+1)
   
   GB <- apply(xx,2,function(x) {
                  s[1] <- x[1]
                  b[1] <- x[2] - x[1]
                  for(t in 2:(ryr+1)) {
                    s[t] <- alp*x[t] + (1-alp)*(s[t-1] + b[t-1])
                    b[t] <- bet*(s[t] - s[t-1]) + (1-bet)*b[t-1]
                  }
                  return(s[ryr+1]+b[ryr+1])
                  })
   
   TAC <- (CTAR * GB)/ITAR

   return(TAC)

}

for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    stk <- stk_init
    stk$index[1:(proj_strt+1),] <- stk$index[1:(proj_strt+1),] * index_epsilon[e,1:(proj_strt+1),]
    
    cat(eff_seq[e],ryr_seq[r],'\n')
    
    for(y in (proj_strt+1):(proj_end-1)) {
  
      # control rule
      stk$catch[y,] <- hcr_dex(stk$index,y,ryr_seq[r])
      stk$theta[y,] <- hcr_opt(stk$catch,sr_epsilon,y)
      stk$entropy[y,] <- H(stk,y,ryr_seq[r],cv.pred(eff_seq[e]))
    
      # project
      out <- pdyn(SSB0,stk$catch[1:y,],slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:y,])
      stk$index[y+1,] <- out$bexp[y+1,] * catchability * index_epsilon[e,y,]
      stk$bexp[y+1,]  <- out$bexp[y+1,]
      
      #plot(1:(y+1),out$bexp[1:(y+1),10] * catchability,xlim=c(0,100),type='l')
      #points((y+1-ryr_seq[r]):(y+1),stk$index[(y+1-ryr_seq[r]):(y+1),10])
      #points(proj_strt:y,stk$catch[proj_strt:y,10] * ITAR/CTAR,pch=19,col=2,type='b')
    }
  
    rarr[e,r,] <- efficiency(stk)
    earr[e,r,] <- entropy(stk)

  }  
}

res<-data.frame(effort=rep(eff_seq,length(ryr_seq*nit)),ryr=rep(rep(ryr_seq,each=length(eff_seq)),nit),efficiency=as.vector(rarr),entropy=as.vector(earr))
save(res,earr,rarr,stk,hcr_dex,file='../res/hcr_dex.Rdata')


# HCR V: model based control rule
hcr_mdl <- function(catch,index,year,ryr=1) {

  GB <- numeric(dim(catch)[2])
  
  #plot(1:(100-1),index[1:(100-1),10])
  index[1:(year-ryr-2),] <- NA
  #points((year-ryr-1):(year-1),index[(year-ryr-1):(year-1),10],pch=19,col=2)
  #lines(1:(year),tryCatch(ipred.sra(catch[1:(year-1),10],index[1:(year-1),10],slope,M,mat,sel,wt,amin,amax),error = function(e) return(index[1:year,10])))
  #cat('off we go\n')
  # get predicted CPUE from SRA
  for(i in 1:length(GB)) {
      GB[i] <- tryCatch(ipred.sra(catch[1:(year-1),i],index[1:(year-1),i],slope,M,mat,sel,wt,amin,amax)[year],error = function(e) return(catch[year-1,i] * ITAR/CTAR))
      #cat(i,'ok\n')
      #browser()
  }
  
  TAC <- (CTAR * GB)/ITAR

  return(TAC)
  
}


for(e in 1:length(eff_seq)) {

  for(r in 1:length(ryr_seq)) {

    stk <- stk_init
    stk$index[1:(proj_strt+1),] <- stk$index[1:(proj_strt+1),] * index_epsilon[e,1:(proj_strt+1),]
    
    cat(eff_seq[e],ryr_seq[r],'\n')
    
    for(y in (proj_strt+1):(proj_end-1)) {
    
      # control rule
      stk$catch[y,]   <- hcr_mdl(stk$catch,stk$index,y,ryr_seq[r])
      stk$theta[y,]   <- hcr_opt(stk$catch,sr_epsilon,y)
      stk$entropy[y,] <- H(stk,y,ryr_seq[r],cv.pred(eff_seq[e]))
    
      # project
      out <- pdyn(SSB0,stk$catch[1:y,],slope,M,mat,sel,wt,amin,amax,sr_epsilon[1:y,])
      stk$index[y+1,] <- out$bexp[y+1,] * catchability * index_epsilon[e,y+1,]
      stk$bexp[y+1,]  <- out$bexp[y+1,]
    
      #plot(1:(y+1),out$bexp[1:(y+1),10] * catchability,xlim=c(0,100),type='l')
      #points((y+1-ryr_seq[r]):(y+1),stk$index[(y+1-ryr_seq[r]):(y+1),10])
      #points(proj_strt:y,stk$catch[proj_strt:y,10] * ITAR/CTAR,pch=19,col=2,type='b')

    }
    
    rarr[e,r,] <- efficiency(stk_mdl)
    earr[e,r,] <- entropy(stk_mdl)
  }  
}

res<-data.frame(effort=rep(eff_seq,length(ryr_seq*nit)),ryr=rep(rep(ryr_seq,each=length(eff_seq)),nit),efficiency=as.vector(rarr),entropy=as.vector(earr))
save(res,earr,rarr,stk,hcr_mdl,file='../res/hcr_mdl.Rdata')

windows();boxplot(t(stk_emp$bexp),main='empirical')
windows();boxplot(t(stk_reg$bexp),main='regression')
windows();boxplot(t(stk_mdl$bexp),main='model')


#rmat <- apply(rarr,1:2,mean)
#windows(width=18)
#par(mfrow=c(1,3))
#image(y=ryr_seq,x=eff_seq,z=rmat,ylab='Retrospective years of data',xlab='Survey effort',main='Efficiency',cex.lab=1.5,cex.axis=1.5)
#contour(y=ryr_seq,x=eff_seq,z=rmat,add=T,labels='')
#
#eff_max <- which.max(apply(rmat,1,mean)) 
#ryr_max <- which.max(apply(rmat,2,mean))
#
#abline(v=eff_seq[eff_max],lwd=2,lty=3)
#abline(h=ryr_seq[ryr_max],lwd=2,lty=3)
#
## include marginals
#plot(ryr_seq,rmat[eff_max,],ylab='Efficiency',xlab='Retrospecitve years',main='',type='l')
#plot(eff_seq,rmat[,ryr_max],ylab='Efficiency',xlab='Survey effort',main='',type='l')
#savePlot(file='../res/hcr_mdl.pdf',type='pdf')
#save(rarr,rmat,hcr_mdl,file='../res/hcr_mdl.Rdata')


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

