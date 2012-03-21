
library(FLAssess)
library(FLash)

## Assessment upto and including 2001
data(ple4)
ss <- stf(ple4,20)

ii <- FLIndex(FLQuant(dimnames=dimnames(catch(ss)),dim=dim(catch(ss))))
index.q(ii)[] <- 1e-4

index(ii)[] <- index.q(ii) * quantSums(stock.n(ss) * stock.wt(ss) * (1 - exp(-harvest(ss)-m(ss)/2)))
index(ii)


hcr <- function(catch,index,year) {

  TAC <- as.numeric(catch[,as.character(year-1)])

  ctrl <- fwdControl(data.frame(year=year,val=TAC,quantity="catch"))
  
}

for(y in 2009:2028) {

  ctrl <- hcr(catch(ss),index(ii),y)
  ss <- fwd(ss,ctrl=ctrl,sr=list(model="mean", params=FLPar(25000)))
}
plot(FLStocks(ss))    










# set courtship and egg laying in Autumn
black.bird@m.spwn[]      <-0.66
black.bird@harvest.spwn[]<-0.66

# assessment is in year 2002, set catch constraint in 2002 and a first guess for F in 2003
ctrl          <-fwdControl(data.frame(year=2002:2003,val=c(85000,.5),quantity=c("catch","f")))
black.bird    <-fwd(black.bird, ctrl=ctrl, sr=list(model="mean", params=FLPar(25000)))

# HCR specifies F=0.1 if ssb<100000, F=0.5 if ssb>300000
# otherwise linear increase as SSB increases
min.ssb<-100000
max.ssb<-300000
min.f  <-0.1
max.f  <-0.5

# slope of HCR
a.    <-(max.f-min.f)/(max.ssb-min.ssb)
b.    <-min.f-a.*min.ssb

# plot of HCR
plot(c(0.0,min.ssb,max.ssb,max.ssb*2),c(min.f,min.f,max.f,max.f),type="l",ylim=c(0,max.f*1.25),xlim=c(0,max.ssb*2))

## find F through iteration
t.    <-999
i     <-0
while (abs(ctrl@target[2,"val"]-t.)>10e-6 & i<50)
   {
   t.<-ctrl@target[2,"val"]  ## save last val of F

   # calculate new F based on SSB last iter
   ctrl@target[2,"val"]    <-a.*c(ssb(black.bird)[,"2010"])+b.
   ctrl@trgtArray[2,"val",]<-a.*c(ssb(black.bird)[,"2010"])+b.
   black.bird<-fwd(black.bird, ctrl=ctrl, sr=list(model="mean", params=FLPar(25000)))

   # 'av a gander
   points(c(ssb(black.bird)[,"2010"]),c(ctrl@target[2,"val"]),cex=1.25,pch=19,col=i)
   print(c(ssb(black.bird)[,"2010"]))
   print(c(ctrl@target[2,"val"]))
   i<-i+1
   }

# F bounds
black.bird      <-fwd(black.bird, ctrl=ctrl, sr=list(model="mean",params=FLPar(25000)))
plot(FLStocks(black.bird))