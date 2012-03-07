
dat.raw <- read.csv('cod_IV.csv')
years <- unique(dat.raw$year)
n.tow <- numeric(length(years))

dat <- new('list')

for(y in 1:length(years)) {

  dat[[y]] <- subset(dat.raw,year==years[y])$n_per_tow
  n.tow[y] <- length(dat[[y]])
}

cv <- function(x, ...) sqrt(var(x, ...))/mean(x, ...)

n.bits <- 1000
boot <- numeric(n.bits)
eff.seq <- seq(100,max(n.tow),10)
dat.cv <- array(NA,dim=c(length(years),length(eff.seq)),dimnames=list(year=years,effort=eff.seq))

for(y in 1:length(years)) {
  
  for(n in 1:length(eff.seq)) {
  
    if(n.tow[y]>=n) {
      for(i in 1:1000) boot[i] <- mean(sample(dat[[y]],size=eff.seq[n],replace=T))
      dat.cv[y,n] <- cv(boot,na.rm=T)
    }
  }
}

save(dat,dat.cv,file='NSIBTS_cod_surv_dat.Rdata')
boxplot(dat.cv,outline=F,ylab=expression(paste('CV [',italic(hat(I)),']')),xlab='Survey effort (annual no. of tows)',main='')
savePlot('effCV.pdf',type='pdf')
dev.off()


