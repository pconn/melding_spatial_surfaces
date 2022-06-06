# plot SSL maps for ms.
library(ggplot2)
library(RColorBrewer)
library(sf)


load("output_rr_bias.RData")
load('./data_integration/data/fitLangevin.RData') #only needed for Grid_sf
Grid_locs = as.matrix(st_coordinates(st_centroid(Grid_sf)))

myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))

get_pi <- function(x) { exp(x) / sum(exp(x))}


#we're going to need to use grid_arrange since facet_grid, etc. doesn't allow different 
#legends for each plot

log_lims = c(-35.2,15.6) #includes all log-scale rel abundance values
real_lims = c(0,0.012)
n_s = nrow(Data_SSL$Y_i)
my_theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),
               plot.title = element_text(hjust = 0, size=10)) 


scale_log <-  scale_fill_gradientn(colours=myPalette(4),limits=log_lims,
                                   breaks = c(-30,-15,0,15),values=scales::rescale(c(-.9,-0.3,.1,.4)))
#breaks = c(-30,-5,-2,0,2,5,15),values=scales::rescale(c(-.8, -.45, -0.3, -.2,-.1,0,0.05)))
scale_real <- scale_fill_gradientn(colours=myPalette(100),limits=real_lims, 
                                   breaks=c(0,0.003,0.006,0.009,0.012),values=scales::rescale(c(-0.9,-0.8,-0.7,1.0)))


Plot_DF = data.frame( Easting = round(Grid_locs[,1],8),
                      Northing = round(Grid_locs[,2], 8),
                      Rel_abund = Data_SSL$Y_i[,1])
POP_log_plot = ggplot(Plot_DF) + aes(x=Easting,y=Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_log + labs(fill = "Log") + ggtitle('POP') 


Plot_DF$Rel_abund = get_pi(Data_SSL$Y_i[,1])
POP_abs_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_real+ labs(fill = expression(pi)) + ggtitle(' ')
  
Plot_DF$Rel_abund = Data_SSL$Y_i[,2]
UD_log_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_log+ labs(fill = "Log") + ggtitle('UD') 

Plot_DF$Rel_abund = get_pi(Data_SSL$Y_i[,2])
UD_abs_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_real+ labs(fill = expression(pi)) + ggtitle(' ')

Plot_DF$Rel_abund = Report$Mu_pred
Biased_log_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_log+ labs(fill = "Log") + ggtitle('Biased') 

Plot_DF$Rel_abund = get_pi(Report$Mu_pred)
Biased_abs_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_real+ labs(fill = expression(pi)) + ggtitle(' ')

load("output_rr_nobias.RData")
Plot_DF$Rel_abund = Report$Mu_pred
Unbiased_log_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_log+ labs(fill = "Log") + ggtitle('Unbiased') 

Plot_DF$Rel_abund = get_pi(Report$Mu_pred)
Unbiased_abs_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_real+ labs(fill = expression(pi))+ ggtitle(' ')

Plot_DF$Rel_abund = 0.5 * get_pi(Data_SSL$Y_i[,1]) + get_pi(Data_SSL$Y_i[,2])
Mean_abs_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_real+ labs(fill = expression(pi))+ ggtitle(' ')

Plot_DF$Rel_abund = log(Plot_DF$Rel_abund) - mean(log(Plot_DF$Rel_abund))
Mean_log_plot = ggplot(Plot_DF) + aes(Easting,Northing,fill=Rel_abund)+geom_raster()+my_theme+
  scale_log+ labs(fill = "Log")+ ggtitle('Mean') 

### now for the awesome big paneled plot

pdf("SSL_abund_plots.pdf",height=8,width=6)
gridExtra::grid.arrange(POP_log_plot,POP_abs_plot,
                        UD_log_plot,UD_abs_plot,
                        Mean_log_plot,Mean_abs_plot,
                        Biased_log_plot,Biased_abs_plot,
                        Unbiased_log_plot,Unbiased_abs_plot,
                        ncol=2)
dev.off()


