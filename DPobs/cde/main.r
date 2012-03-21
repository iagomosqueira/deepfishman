# The main loop

# Load the preconditioned stock
# Set parameters (obs error, srr error, nyrs etc)
# Loop over years and run fowards
# Calc the performance stats

#*******************************************************************************
rm(list=ls())
library(FLCore)
library(FLash)
library(FLAssess)
library(FLsp)

#*******************************************************************************
setwd("c:/Projects/Deepfishman/deepfishman/trunk/DPobs/R")
source("hcr_func.r")

# Load the conditioned stock
source("make_stocks.r")

#*******************************************************************************
# Set some parameters
obs_cv <- 0.6
srr_sd <- 0
niters <- 1
projyrs <- 20

#Blim <- 0

hcryrs <- (dims(stk1)$year+1) : (dims(stk1)$year+projyrs)

#*******************************************************************************
# Set up the stk for projection
# Including future weights and all that gubbins
stk_true <- stf(stk1,nyears=projyrs)
# Clean out harvest to avoid confusion
harvest(stk_true)[,projyrs] <- NA
#stk_true <- window(stk1,end=dims(stk1)$year+projyrs)
stk_true <- propagate(stk_true,niters)

# To store output of the HCR
TAC <- array(NA,dim=c(niters,dims(stk_true)$year))

index_sp <- stock(stk_true)
index_sp <- index_sp * rlnorm(prod(dim(index_sp)),0,obs_cv)

for (yr in hcryrs)
{
  #if (yr==20) browser() # something weird going on...

  cat("yr: ", yr, "\n")
  # Get new TAC
  # Assess using last years catch and index data
  # No error on catch
  catch_sp <- catch(stk_true)[,1:(yr-1)]
  index_sp[,yr-1] <- stock(stk_true)[,yr-1] * rlnorm(niters,0,obs_cv)
  # Plus some obs error!
  #index_sp * rlnorm(prod(dim(index_sp)),0,obs_cv)
  cat("index: ", index_sp, "\n")
  
  flsp <- FLsp(catch=catch_sp,index=window(index_sp, start=1,end=yr-1))
  flsp <- fitsp(flsp,lower=c(1e-6,10),upper=c(10,1e9),control=DEoptim.control(itermax=1500,trace=500))
  Bcurrent <- bcurrent(flsp)
  TAC[,yr] <- hcr(B=c(Bcurrent), r= c(params(flsp)['r',]), k= c(params(flsp)['k',]), Blim=0.9*Bmsy(flsp))
  cat("New TAC: ", TAC[,yr], "\n")

  # Project stock using TAC
  ctrl <- fwdControl(data.frame(year=yr, quantity="catch",val=TAC[,yr]))
  stk_true <- fwd(stk_true, ctrl=ctrl, sr=list(model=model(cod1), params=params(cod1)))

}

plot(stk_true)
# final estimate of msy
Msy(flsp)
