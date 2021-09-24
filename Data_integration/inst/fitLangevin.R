#remotes::install_github("bmcclintock/momentuHMM@develop")
library(momentuHMM)
library(raster)
library(sp)
library(sf)
library(Brobdingnag)
#remotes::install_github("papayoun/Rhabit")
library(Rhabit)
load("./data_integration/data/ssl_habitat.RData")
load("./data_integration/data/SSL_crawl.RData")

tracks <- prepData(clipcrwData,coordNames=c("mu.x","mu.y"),altCoordNames = "mu")

# convert to km ######
hab_cov$depth <- hab_cov$depth/1000
hab_cov$dist_site = hab_cov$dist_site/1000
hab_cov$dist_land = hab_cov$dist_land/1000

covlist0 <- list(depth = hab_cov$depth,
                 slope = hab_cov$slope,
                 d2site = hab_cov$dist_site,
                 d2land = hab_cov$dist_land)

tracks$mu.x <- tracks$mu.x/1000
tracks$mu.y <- tracks$mu.y/1000

for(i in 1:length(covlist0)) {
  extent(covlist0[[i]]) <- extent(c(xmin(covlist0[[i]]), xmax(covlist0[[i]]), 
                                    ymin(covlist0[[i]]), ymax(covlist0[[i]]))/1000)
  projection(covlist0[[i]]) <- gsub("units=m", "units=km", projection(covlist0[[i]]))
}
extent(hab_cov) <- extent(c(xmin(hab_cov), xmax(hab_cov), 
                            ymin(hab_cov), ymax(hab_cov))/1000)
projection(hab_cov) <- gsub("units=m", "units=km", projection(hab_cov))
####################

# Crop covariates to area of interest
border <- 30
lim <- c(min(tracks$mu.x)-border,max(tracks$mu.x)+border,min(tracks$mu.y)-border,max(tracks$mu.y)+border)
covlist0 <- lapply(covlist0, crop, y=extent(lim))

covlist <- lapply(covlist0, Rhabit::rasterToRhabit)
gradarray <- Rhabit::bilinearGradArray(locs = as.matrix(tracks[c("mu.x","mu.y")]), cov_list = covlist)

covNames <- c("depth","slope","d2site","d2land")
gradNames <- c("x","y")
for(i in 1:length(covNames)){
  for(j in 1:length(gradNames)){
    tracks[[paste0(covNames[i],".",gradNames[j])]] <- gradarray[,j,i]
  }
}

# specify pseudo-design matrix
DM <- list(mu=matrix(c("mu.x_tm1","langevin(depth.x)","langevin(slope.x)","langevin(d2site.x)",0,0,
                       "mu.y_tm1","langevin(depth.y)","langevin(slope.y)","langevin(d2site.y)",0,0,
                                0,                  0,                  0,                   0,1,0,
                                0,                  0,                  0,                   0,0,1,
                                0,                  0,                  0,                   0,1,0),
                     5,6,byrow=TRUE,dimnames=list(c("mean.x","mean.y","sigma.x","sigma.xy","sigma.y"),
                                                  c("mean:mu_tm1","depth","slope","d2site","sigma:(Intercept)","sigma.xy:(Intercept)"))))

# fit model
fitLangevin <- fitCTHMM(tracks,Time.name="date_time",Time.unit="hours",dist=list(mu="rw_mvnorm2"),nbStates=1,DM=DM,
                             Par0=list(mu=c(1,0.76804652,  0.05470275, -0.07548756, -0.37966980,0)),
                             fixPar=list(mu=c(1,rep(NA,3),NA,0)))


# calculate and plot logUD for entire study area
logUD <- fitLangevin$CIbeta$mu$est[2]*hab_cov$depth+fitLangevin$CIbeta$mu$est[3]*hab_cov$slope+fitLangevin$CIbeta$mu$est[4]*hab_cov$dist_site
values(logUD) <- getValues(logUD) - log(sum(Brobdingnag::brob(getValues(logUD)[!is.na(getValues(logUD))])))
par(mfrow=c(1,2))
raster::plot(logUD,main="logUD") 
plot(exp(logUD),main="UD") 

# plot UD and habitat covariate
par(mfrow=c(2,2))
plot(logUD,main="logUD")
plot(hab_cov$depth,main="depth")
plot(hab_cov$slope,main="slope")
plot(hab_cov$dist_site,main="d2site")

# plot tracks and logUD
par(mfrow=c(1,1))
plotSpatialCov(tracks,logUD)

# extract coefficients and VC matrix
UDcoeffs <- fitLangevin$CIbeta$mu$est[,c(2,3,4)] # RSF coefficients (negative = select for smaller values of covariate)
UDvar <- fitLangevin$mod$Sigma[c(2,3,4),c(2,3,4)] # RSF coefficient variance-covariance matrix
rownames(UDvar) <- colnames(UDvar) <- names(UDcoeffs)

save.image("fitLangevin.RData")
