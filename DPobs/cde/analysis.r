
load('../res/hcr_mav.Rdata')

windows()
xyplot(entropy~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Model-based HCR')
savePlot('../res/hcr_mav_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Model-based HCR')
savePlot('../res/hcr_mav_plot2.pdf',type='pdf')
xyplot(log(efficiency)~entropy,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Model-based HCR')
savePlot('../res/hcr_mav_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='mav')
res_all <- res


load('../res/hcr_reg.Rdata')

windows()
xyplot(entropy~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Regression-based HCR')
savePlot('../res/hcr_reg_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Regression-based HCR')
savePlot('../res/hcr_reg_plot2.pdf',type='pdf')
xyplot(log(efficiency)~entropy,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Regression-based HCR')
savePlot('../res/hcr_reg_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='reg')
res_all <- rbind(res_all,res)


load('../res/hcr_smt.Rdata')

windows()
xyplot(entropy~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Smoothed CPUE HCR')
savePlot('../res/hcr_smt_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Smoothed CPUE HCR')
savePlot('../res/hcr_smt_plot2.pdf',type='pdf')
xyplot(log(efficiency)~entropy,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Smoothed CPUE HCR')
savePlot('../res/hcr_smt_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='smt')
res_all <- rbind(res_all,res)

load('../res/hcr_dex.Rdata')

windows()
xyplot(entropy~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Double exponential HCR')
savePlot('../res/hcr_dex_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Double exponential HCR')
savePlot('../res/hcr_dex_plot2.pdf',type='pdf')
xyplot(log(efficiency)~entropy,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Double exponential HCR')
savePlot('../res/hcr_dex_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='dex')
res_all <- rbind(res_all,res)

load('../res/hcr_mdl.Rdata')

windows()
xyplot(entropy~effort|as.factor(ryr),res,pch=19,panel='panel.loess',main='Model-based HCR')
savePlot('../res/hcr_mdl_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Model-based HCR')
savePlot('../res/hcr_mdl_plot2.pdf',type='pdf')
xyplot(log(efficiency)~entropy,res,panel=function(x,y) {panel.xyplot(x,y);panel.loess(x,y,col=2,lwd=3)},main='Model-based HCR')
savePlot('../res/hcr_mdl_plot3.pdf',type='pdf')

res     <- data.frame(res,hcr='mdl')
res_all <- rbind(res_all,res)

windows()
xyplot(entropy~effort|as.factor(ryr),res_all,groups=hcr,panel='panel.xyplot',main='HCR comparison')
savePlot('../res/hcr_all_plot1.pdf',type='pdf')
xyplot(log(efficiency)~effort|as.factor(ryr),res_all,groups=hcr,panel='panel.xyplot',main='HCR comparison')
savePlot('../res/hcr_all_plot2.pdf',type='pdf')
xyplot(log(efficiency)~entropy,res_all,groups=hcr,panel='panel.xyplot',main='HCR comparison')
savePlot('../res/hcr_all_plot3.pdf',type='pdf')



