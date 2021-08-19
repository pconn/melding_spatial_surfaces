# summarize results of simulation study

load('sim_RMSE.RData')

arr_names = vector("list",4)
arr_names[[1]]=c("2","4","8")
arr_names[[2]]=c("All biased","One unbiased","All unbiased")
arr_names[[3]]=c("All biased","All unbiased","One unbiased","Mean","First surface")
arr_names[[4]]=as.character(c(1:500))
RMSE_all = array(NA,dim=c(3,3,5,500),dimnames=arr_names)
RMSE_all[,,,1:100]=RMSE

load('sim_RMSE_bat1.RData')
RMSE_all[,,,101:200]=RMSE
load('sim_RMSE_bat2.RData')
RMSE_all[,,,201:300]=RMSE
load('sim_RMSE_bat3.RData')
RMSE_all[,,,301:400]=RMSE
load('sim_RMSE_bat4.RData')
RMSE_all[,,,401:500]=RMSE

RMSE = apply(RMSE_all,c(1,2,3),'mean',na.rm=T)

SE = sqrt(apply(RMSE_all,c(1,2,3),'var',na.rm=T))/sqrt(500)
SE_df = reshape::melt(SE,varnames=c("Num.surf","Bias.type","Estimator"))

library(ggplot2)
Plot_df = reshape::melt(RMSE,varnames=c("Num.surf","Bias.type","Estimator"))
colnames(Plot_df)[4]="RMSE"
Plot_df$SE = SE_df$value

bar_plot = ggplot(Plot_df)+geom_bar(stat="identity",aes(x=Estimator,y=RMSE,fill=Estimator))+
  facet_grid(Num.surf~Bias.type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  geom_errorbar(aes(x=Estimator,ymin=RMSE-2*SE,ymax=RMSE+2*SE))

pdf("RMSE_bar_plots.pdf")
  bar_plot
dev.off()


