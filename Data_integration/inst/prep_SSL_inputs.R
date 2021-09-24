# integrate_SSL.R  integrate POP and UDs for SSL

# load UD info 
library(terra)
library(sf)
library(ggplot2)
library(INLA)

load('./data_integration/data/fitLangevin.RData')

hab_cov = values(hab_cov)

raster_sf = st_as_sf(raster::rasterToPolygons(logUD))
raster_sf <- st_transform(raster_sf, crs=crs(Grid_sf))

grid_outline = st_union(Grid_sf)
raster_centroids = st_centroid(raster_sf)
I_keep = st_intersects(raster_centroids,grid_outline,sparse=FALSE)

#plot raster, study area together (for initial debugging, prob no longer of interest)
grid_plot = ggplot()+geom_sf(data=raster_sf,color='blue') 
grid_plot + geom_sf(data=Grid_sf)

hab_cov = hab_cov[which(I_keep==1),]

X = cbind(hab_cov[,"depth"],hab_cov[,"slope"],hab_cov[,"dist_site"])

logUD = values(logUD)[which(I_keep==1)]
logUD = logUD-mean(logUD)  #standardize so has mean 0...better for optimization
VC_logUD = X %*% (UDvar %*% t(X))
diag(VC_logUD) = diag(VC_logUD)+0.000001  #nugget regularization - makes sure proper var-cov matrix
Cor_logUD = cov2cor(VC_logUD)

#load POP data
load("./data_integration/data/POPfits.rda")
log_POPfits = log(POPfits[,3:1002])
logPOP = rowMeans(log_POPfits)
VC_logPOP = cov(t(log_POPfits))
Cor_logPOP = cov2cor(VC_logPOP)
diag(VC_logPOP)=diag(VC_logPOP)+0.000001

VCinv_logUD = solve(VC_logUD)
VCinv_logPOP = solve(VC_logPOP)

n_s = length(logPOP)
Grid_locs = data.frame(POPfits[,1:2])
mesh = inla.mesh.create( Grid_locs )
n_knots = mesh$n
Eta_index = mesh$idx$loc-1  #which spde REs to apply as random effects for each cell centroid
#note it might be a lot faster to just to model the individual data points instead of the whole grid
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]


Data_SSL <- list("M0"=spde$M0,"M1"=spde$M1,"M2"=spde$M2,
                 "flag"=1,"n_s"=n_s,n_mesh=n_knots,
                 "options"=c(0,0),"Mesh_index"=Eta_index,"matern_pri"=c(0,0,0,0),
                 "n_surf"=2)
Data_SSL$Y_i = cbind(logPOP,logUD)
Data_SSL$Omega_Y = array(dim=c(n_s,n_s,2))   #inverse of Sigma needed for fast dmvnorm-type calcs
Data_SSL$Omega_Y[,,1]=VCinv_logPOP
Data_SSL$Omega_Y[,,2]=VCinv_logUD

save(Data_SSL,file="SSL_TMB_data.RData")


