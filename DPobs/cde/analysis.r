
load('../res/hcr_mav.Rdata')

windows()
xyplot(uncertainty~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Moving average HCR')
savePlot('../res/hcr_mav_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Moving average HCR')
savePlot('../res/hcr_mav_plot2.pdf',type='pdf')
xyplot(log(efficiency)~uncertainty,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Moving average HCR')
savePlot('../res/hcr_mav_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='Moving average')
res_all <- res


load('../res/hcr_reg.Rdata')

windows()
xyplot(uncertainty~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Regression-based HCR')
savePlot('../res/hcr_reg_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Regression-based HCR')
savePlot('../res/hcr_reg_plot2.pdf',type='pdf')
xyplot(log(efficiency)~uncertainty,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Regression-based HCR')
savePlot('../res/hcr_reg_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='Linear regression')
res_all <- rbind(res_all,res)


load('../res/hcr_smt.Rdata')

windows()
xyplot(uncertainty~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Smoothed CPUE HCR')
savePlot('../res/hcr_smt_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Smoothed CPUE HCR')
savePlot('../res/hcr_smt_plot2.pdf',type='pdf')
xyplot(log(efficiency)~uncertainty,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Smoothed CPUE HCR')
savePlot('../res/hcr_smt_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='Smoothed (Exponential)')
res_all <- rbind(res_all,res)

load('../res/hcr_dex.Rdata')

windows()
par(cex.lab=2)
xyplot(uncertainty~effort|as.factor(ryr),res,pch=19,panel='panel.loess',par.settings=list(cex.lab=2))#,main='Double exponential HCR')
savePlot('../res/hcr_dex_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Double exponential HCR')
savePlot('../res/hcr_dex_plot2.pdf',type='pdf')
xyplot(log(efficiency)~uncertainty,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Double exponential HCR')
savePlot('../res/hcr_dex_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='Smoothed')
res_all <- rbind(res_all,res)

load('../res/hcr_mdl.Rdata')

windows()
xyplot(uncertainty~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Model-based HCR')
savePlot('../res/hcr_mdl_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Model-based HCR')
savePlot('../res/hcr_mdl_plot2.pdf',type='pdf')
xyplot(log(efficiency)~uncertainty,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Model-based HCR')
savePlot('../res/hcr_mdl_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='Model-based')
res_all <- rbind(res_all,res)

res <- res_all

windows()
xyplot(log(efficiency)~uncertainty|hcr,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='',par.strip.text = list(cex = 1.75))
savePlot('../res/hcr_all_plot.pdf',type='pdf')


save(res,file='../res/hcr_all.Rdata')



