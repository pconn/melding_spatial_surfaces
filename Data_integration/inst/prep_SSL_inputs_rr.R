# prep_SSL_inputs_rr.R  Prepare inputs for combining 
# Steller sea lion POP and UD surfaces, assuming diagonal covariance

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
raster_sf <- raster_sf[I_keep==1,]
raster_centroids = st_centroid(raster_sf)
Inter = st_intersects(raster_centroids,Grid_sf,sparse=FALSE)
Order = rep(NA,nrow(Inter))
for(i in 1:nrow(Inter)){
  Order[i]=which(Inter[,i]==1)
}
#plot_N_map_xy(N=raster_sf$layer[Order],XY=Grid_locs,leg.title="UD") 
log_UD = raster_sf$layer[Order]
hab_cov = hab_cov[Order,]

n_model = nrow(hab_cov)
log_UD = log_UD[1:n_model]
hab_cov = hab_cov[1:n_model,]

X = cbind(hab_cov[,"depth"],hab_cov[,"slope"],hab_cov[,"dist_site"])

logUD = log_UD-mean(log_UD)  #standardize so has mean 0...better for optimization
VC_logUD = X %*% (UDvar %*% t(X))
diag(VC_logUD)=diag(VC_logUD)+0.000001  #diagonal regularization so invertible
VC_inv_logUD_full = solve(VC_logUD)
VCinv_logUD = diag(length(logUD))
diag(VCinv_logUD)=diag(VC_inv_logUD_full)


#load POP data
load("./data_integration/data/POPests.RData")
n_s <- length(logPOP)

Grid_locs = st_coordinates(st_centroid(Grid_sf))
mesh = inla.mesh.create( Grid_locs )
n_knots = mesh$n
Eta_index = mesh$idx$loc-1  #which spde REs to apply as random effects for each cell centroid
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]


Data_SSL <- list("M0"=spde$M0,"M1"=spde$M1,"M2"=spde$M2,
                 "flag"=1,"n_s"=n_s,n_mesh=n_knots,
                 "options"=c(0,0),"Mesh_index"=Eta_index,"matern_pri"=c(0,0,0,0),
                 "n_surf"=2)
Data_SSL$Y_i = cbind(logPOP,logUD)
Data_SSL$Omega_Y1 = as(VCinv_logPOP,"dgTMatrix")
Data_SSL$Omega_Y2 = as(VCinv_logUD,"dgTMatrix")


#this block lifted from sdmTMB package
set.seed(12345)
Grid_locs = Grid_locs/1000
knots <- stats::kmeans(x = Grid_locs, centers = 1000)  #1000 "center" locations (# knots will be more after triangulation)
loc_centers <- knots$centers
mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
spde <- INLA::inla.spde2.matern(mesh)  

A <- INLA::inla.spde.make.A(mesh, loc = as.matrix(Grid_locs)) #prediction matrix
Data_SSL$A = as(A,"dgTMatrix")
Data_SSL$n_mesh = ncol(A)
Data_SSL$M0 = spde$param.inla$M0
Data_SSL$M1 = spde$param.inla$M1
Data_SSL$M2 = spde$param.inla$M2
Data_SSL$Mesh_index = NULL


save(Data_SSL,file="SSL_TMB_data.RData")


