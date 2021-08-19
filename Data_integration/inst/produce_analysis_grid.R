#produce Steller grid
library('sf')
library('sp')

load('POPhexagons.rda')
Hexes_sf = st_as_sf(POPhexagons)
bounding_box = st_bbox(Hexes_sf)


bbox_coords = st_bbox(Hexes_sf)
coord_mat = data.frame(matrix(0,4,2))  # bottom left, top left, top right, bottom right on rows
coord_mat[,1]=c(rep(bbox_coords$xmin,2),rep(bbox_coords$xmax,2))
coord_mat[,2]=c(bbox_coords$ymin,bbox_coords$ymax,bbox_coords$ymax,bbox_coords$ymin)
colnames(coord_mat)=c("x","y")
bbox_pts = st_as_sf(coord_mat,coords=c("x","y"))
st_crs(bbox_pts) <- st_crs(Hexes_sf)
Grid_big = st_make_grid(bbox_pts,cellsize=15000)  #15,000 meters per grid cell length
Grid_big=st_as_sf(Grid_big)

##limit to western DPS
LongLat = st_coordinates(st_centroid(st_transform(Grid_big,crs=4326)))
Which_remove = which(LongLat[,1]<0 & LongLat[,1]>(-144))
Grid_big = Grid_big[-Which_remove,]  #limit to western DPS hexes

#remove cells not intersecting hex grid
Intersects = st_intersects(Grid_big,st_union(Hexes_sf),sparse=FALSE)
Grid_sf = Grid_big[which(Intersects),]

#remove any remaining cells entirely overlapping land & attach proportion that is salt water habitat
AK = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ak_dcw.shp')
AK_union = st_union(AK)
n_cells=nrow(Grid_sf)
Land_area=rep(0,n_cells)
I_intersect = st_intersects(Grid_sf,AK_union)
for(icell in 1:n_cells){
  if(length(I_intersect[[icell]])>0)Land_area[icell]=st_area(st_intersection(Grid_sf[icell,],AK_union))
}
Grid_sf$prop_land=as.numeric(Land_area/st_area(Grid_sf[1,]))
Grid_sf = Grid_sf[-which(Grid_sf$prop_land==1),]

png('Analysis_grid.png')
  plot(Grid_sf)
dev.off()

save(Grid_sf,file="Surface_analysis_grid.RData")







