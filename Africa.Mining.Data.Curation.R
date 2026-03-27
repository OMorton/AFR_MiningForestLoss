library(terra)
library(sf)
sf_use_s2(F)
library(dplyr)
library(ggplot2)
library(stringr)
options(scipen = 999)

#### Prepare masolele data for analysis ####

## take annual raster of forest loss and crop to ivory coast and filter for mining #
years<- c(2001:2020)

for (i in 1:length(years)){
  year<- years[i]
  print(year)
  forest.loss.rast<- rast(paste0("Data/masolele.base.layers/LandUseFollowingDeforestation_dl_30m_20010101_", year, "1231_af_epsg.4326_v1.tif"))
 # filter for only mines (value = 3), and convert value to year of loss
  cote.divoire.annual.mining.loss<- subst(forest.loss.rast, 3, (year-2000), others = NA, filename = paste0("Data/masolele.processed/africa.forest.loss.mining.", year, ".tif"), overwrite = T, datatype='INT1U')
}
tmpFiles(remove = T)
rm(list = ls())
gc()

## combine all annual mine rasters, and retain the lowest value in each cell (i.e. the first year mine is detected)
africa.forest.loss.raster<- rast("Data/masolele.base.layers/LandUseFollowingDeforestation_dl_30m_20010101_20201231_af_epsg.4326_v1.tif")
example.rast<- rast("Data/masolele.processed/2020.tif")
mine.rasters<- list.files(path = "Data/masolele.processed/", pattern = "mining")
mine.rasters.list<- list()
for (i in 1:length(mine.rasters)){
  filename<- mine.rasters[i]
  year<- 2000+i
  print(year)
  #mine.raster<- rast(paste0("Data/masolele.processed/", filename))
  #mine.raster <- extend(mine.raster, ext(-26.7182321694533, 54.3172175592647, -30.7000146663131, 18.015713023017), filename = paste0("Data/masolele.processed/extended/africa.forest.loss.mining.", year, ".extended.tif"), overwrite = T, datatype='INT1U')
 # mine.raster.new<- resample(mine.raster, africa.forest.loss.raster, filename = paste0("Data/masolele.processed/resampled/africa.forest.loss.mining.", year, ".resampled.tif"), overwrite = T, datatype='INT1U')
    mine.raster<- rast(paste0("Data/masolele.processed/resampled/africa.forest.loss.mining.", year, ".resampled.tif"))
  print(ext(mine.raster))
  mine.rasters.list[[i]]<- mine.raster # add to mine raster list
}
mine.rasters<- rast(mine.rasters.list) # convert list of rasters into sprc object for merge

# Merge the rasters using reduce and min
mine.rasters.loss.year <- app(mine.rasters, "min", na.rm = T, filename = "Data/masolele.processed/africa.mining.loss.year.tif", overwrite = T, datatype='INT1U')

mine.rasters.loss.year<- rast("Data/masolele.processed/africa.mining.loss.year.tif")
# mask mine loss by the tmf biome
tmf<- st_read("Data/tropical.moist.forest.biome.shp")
africa.tmf<- crop(vect(tmf), mine.rasters.loss.year)
mine.rasters.loss.year.tmf<- mask(mine.rasters.loss.year, africa.tmf, filename = "Data/masolele.processed/africa.mining.loss.year.tmf.tif", overwrite = T, datatype='INT1U')




#### Convert to points and get mining clusters based on distance ####
country.mining.raster.files<- list.files(path = "Data/country.mining.data/rasters/", pattern = "mining.loss.year.tif")
library(stringr)
library(dbscan)
for (i in 1:length(country.mining.raster.files){
  print(i)
  country.filename<- country.mining.raster.files[i]
  country.mining.raster<- rast(paste0("Data/country.mining.data/rasters/", country.filename))
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  # convert to points
  # convert to points for density clustering - converts from a raster to points so I can do the density clustering of mines in ArcGIS
 # mine.loss.year.points<- st_as_sf(as.points(country.mining.raster, values = T, na.rm = T))
  #st_write(mine.loss.year.points, paste0("Data/country.mining.data/points/", country.n, ".mining.loss.year.points.shp"), append = F)
  # distance based density clustering
  mine.loss.year.points<- st_read(paste0("Data/country.mining.data/points/", country.n, ".mining.loss.year.points.shp"))
  mine.loss.year.points.proj<- st_transform(st_as_sf(mine.loss.year.points), "ESRI:54009")
  # Convert the sf object to a matrix of coordinates (for dbscan)
  coords <- st_coordinates(mine.loss.year.points.proj)

  # Run DBSCAN
  dbscan_result <- dbscan(coords, eps = 1000, minPts = 5)
  
  # Add cluster assignments to the sf object
  mine.loss.year.points$cluster <- dbscan_result$cluster
  #mine.loss.year.points<- rename(mine.loss.year.points, "year" = "min")
  # remove the noise (cluster = 0)
  mine.loss.year.points<- filter(mine.loss.year.points, cluster > 0)
  colnames(mine.loss.year.points)[1]<- "year"
  mines.cluster.first.loss.values<- mine.loss.year.points %>% 
    st_drop_geometry() %>% 
    group_by(cluster) %>% 
    slice_min(year, n = 1) %>% 
    sample_n(1) %>% 
    rename("first.mine.year" = "year", "CLUSTER_ID" = "cluster")
  
  mine.loss.year.points.vect<- vect(mine.loss.year.points)
  # convert back to 30m raster - because we need the mines as polygons not points
  mines.loss.clusters.raster.mine.id<- rasterize(mine.loss.year.points.vect,country.mining.raster, field = "cluster", filename = paste0("Data/tmp.rasters/", country.n, "tmp.raster.tif"), overwrite = T)
  # this converts the point data back into raster, with 'mine.rasters.loss.year' as the template for conversion
  mines.loss.clusters.raster.mine.id.polygons.vect<-as.polygons(mines.loss.clusters.raster.mine.id, aggregate = T, values = T) 
  mines.loss.clusters.raster.mine.id.polygons<- st_as_sf(mines.loss.clusters.raster.mine.id.polygons.vect) # convert into SF polygons
  # add loss year value into mine polygons
  colnames(mines.loss.clusters.raster.mine.id.polygons)[1]<- "CLUSTER_ID"
  mines.loss.clusters.raster.mine.id.polygons<- left_join(mines.loss.clusters.raster.mine.id.polygons, mines.cluster.first.loss.values)
  # export
  st_write(mines.loss.clusters.raster.mine.id.polygons, paste0("Data/country.mining.data/country.mining.cluster.shapefiles/", country.n, ".mining.clusters.shp"), append = F)
  rm(list = ls()[ls() != "country.mining.raster.files"])
  gc()
  tmpFiles(remove = T)
}
mine.loss.year.points.sf<- st_as_sf(mine.loss.year.points)


st_write(mine.loss.year.points, paste0("Data/country.mining.data/points/", country.n, ".mining.loss.year.points.shp"), append = F)



#### Retain only mines in areas with >50% tree cover in 2000 ####


### Loop to retain only mines where >1/3 of the 5km buffer is dense forest ###

country.mine.files<- list.files(path = "Data/country.mining.data/country.mining.cluster.shapefiles/", pattern = ".mining.clusters.shp")
hansen.forest.cover.2000<- rast("Data/hansen.forest.cover/sub_saharan_Africa_2000_tree_cover_30m.tif")
for (j in 1:length(country.mine.files)){
  print(j)
  country.filename<- country.mine.files[j]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  # get country raster
  country.mining.raster<- rast(paste0("Data/country.mining.data/rasters/", country.n, ".mining.loss.year.tif"))
  # crop forest cover raster to this country
  country.forest.cover<- crop(hansen.forest.cover.2000, country.mining.raster, filename = paste0("Data/hansen.forest.cover/country.forest.cover/", country.n, "_2000_tree_cover_30m.tif"), overwrite = T, datatype='INT1U')
  # get country mine shapefiles
  country.mining.clusters<- st_read(paste0("Data/country.mining.data/country.mining.cluster.shapefiles/", country.filename))
  colnames(country.mining.clusters)[1]<- "CLUSTER_ID"
  colnames(country.mining.clusters)[2]<- "first.mine.year"
  country.mining.clusters.vect<- vect(country.mining.clusters)
  mine.5km.buffers<- st_as_sf(buffer(country.mining.clusters.vect, 5000))
  buffer.forest.prop.list<- list()
  for (i in 1:nrow(mine.5km.buffers)){
  buffer<- mine.5km.buffers[i,]
mine.n<- buffer$CLUSTER_ID
print(mine.n)
# crop and mask forest cover raster to match buffer
forest.cover.cropped<- crop(country.forest.cover, buffer)
forest.cover.masked<- mask(forest.cover.cropped, buffer)
# get area of forest in buffer
buffer.forest.area<- as.numeric(expanse(forest.cover.masked, unit = "m"))[2]
# divide by buffer area
buffer.forest.prop<- buffer.forest.area / (as.numeric(expanse(vect(buffer), unit = "m")))
# turn into dataframe
buffer.forest.prop.df<- data.frame(CLUSTER_ID = mine.n, forest.prop = buffer.forest.prop)
buffer.forest.prop.list[[i]]<- buffer.forest.prop.df
}
# combine into one df
buffer.forest.prop.df<- do.call(rbind, buffer.forest.prop.list)

# combine forest prop with mine data
country.mining.clusters$buffer.forest.prop<- buffer.forest.prop.df$forest.prop
# remove mines where surrounding 5km are not >33% forest
country.mining.clusters.forest <- filter(country.mining.clusters, buffer.forest.prop >= .33)
st_write(country.mining.clusters.forest, paste0("Data/country.mining.data/country.forest.mine.shapefiles/",  country.n, ".forest.mining.clusters.shp"), append = F)
tmpFiles(remove = T)
}

#### Create buffers and rings around forest mines for each country ####


rm(list = ls())
gc()
tmpFiles(remove = T)
country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")

for (j in 1:length(country.forest.mine.files){
  print(j)
  country.filename<- country.forest.mine.files[j]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  
mine.clusters<- st_read(paste0("Data/country.mining.data/country.forest.mine.shapefiles/", country.filename))
print(paste0(country.n, " - ", nrow(mine.clusters), " Mines"))
colnames(mine.clusters)[1]<- "CLUSTER_ID"
colnames(mine.clusters)[2]<- "first.mine.year"
colnames(mine.clusters)[3]<- "buffer.forest.prop"
country.mining.raster<- rast(paste0("Data/country.mining.data/rasters/", country.n, ".mining.loss.year.tif"))
mine.centroids.list<- list()
mine.1km.buffer.list<- list()
mine.5km.buffer.list<- list()
mine.2.5km.buffer.list<- list()
mine.5km.buffer.list<- list()
mine.10km.buffer.list<- list()
mine.20km.buffer.list<- list()
mine.5km.buffer.master.list<- list()
mine.10km.buffer.master.list<- list()

for (i in 1:nrow(mine.clusters)){ # this loop goes through each mine cluster, creates all the buffers, then saves them into the list
  mine<- mine.clusters[i,]
  print(paste0("-------------MINE ", i, "-----------"))
  mine.area<- sum(as.numeric(st_area(mine)))
  # mask mine loss raster by cluster
  mine.losses<- crop(country.mining.raster, vect(mine))
  mine.losses<- mask(mine.losses, vect(mine))
  mine.year<- mine$first.mine.year
  
  ## get centroid of first area to make rings by ##
  #plot(mine.losses)
  mine.centroids<- subst(mine.losses, mine.year, 1, others = NA)
  #plot(mine.centroids)
  mine.centroids.vect<- as.polygons(mine.centroids)
  mine.centroids.point<- st_as_sf(centroids(mine.centroids.vect, inside=TRUE))
  mine.centroids.point<- st_transform(mine.centroids.point, "ESRI:54009")
  mine.centroids.point<- cbind(mine.centroids.point, st_drop_geometry(mine)) %>% select(CLUSTER_ID, first.mine.year, buffer.forest.prop)

  ## create buffers around centroid of first year
  # 0-1km
  mine.1km.buffer<- st_as_sf(buffer(vect(mine.centroids.point), 1000))
  mine.1km.buffer<- mine.1km.buffer %>% 
    mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)

  
  # 1-2.5km
 mine.2.5km.buffer<- st_as_sf(buffer(vect(mine.1km.buffer), 1500))
#  #plot(mine.2.5km.buffer)
  # remove old buffer
  mine.2.5km.buffer<- st_difference(mine.2.5km.buffer, mine.1km.buffer)
  mine.2.5km.buffer<- mine.2.5km.buffer %>% 
    mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)

  
  # 2.5km-5km
 mine.5km.buffer<- st_as_sf(buffer(vect(mine.2.5km.buffer), 2500))
#  # remove old buffer
  mine.5km.buffer<- st_difference(mine.5km.buffer, mine.1km.buffer)
  mine.5km.buffer<- st_difference(mine.5km.buffer, mine.2.5km.buffer)
  mine.5km.buffer<- mine.5km.buffer %>% 
    mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)#plot(mine.5km.buffer)

  #1-5km
  mine.5km.buffer<- st_as_sf(buffer(vect(mine.1km.buffer), 4000))
  # remove old buffer
  mine.5km.buffer<- st_difference(mine.5km.buffer, mine.1km.buffer)
  mine.5km.buffer<- mine.5km.buffer %>% 
    mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)#plot(mine.5km.buffer)
  
  # 5-10km
 mine.10km.buffer<- st_as_sf(buffer(vect(mine.5km.buffer), 5000))
  #plot(mine.10km.buffer)
  # remove old buffer
  mine.10km.buffer<- st_difference(mine.10km.buffer, mine.1km.buffer)
  mine.10km.buffer<- st_difference(mine.10km.buffer, mine.2.5km.buffer)
  mine.10km.buffer<- st_difference(mine.10km.buffer, mine.5km.buffer)
  mine.10km.buffer<- mine.10km.buffer %>% 
    mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)#plot(mine.10km.buffer)

  # 10-20km
  mine.20km.buffer<- st_as_sf(buffer(vect(mine.10km.buffer), 10000))
  #plot(mine.20km.buffer)
  # remove old buffer
  mine.20km.buffer<- st_difference(mine.20km.buffer, mine.1km.buffer)
  mine.20km.buffer<- st_difference(mine.20km.buffer, mine.2.5km.buffer)
  mine.20km.buffer<- st_difference(mine.20km.buffer, mine.5km.buffer)
  mine.20km.buffer<- st_difference(mine.20km.buffer, mine.10km.buffer)
  mine.20km.buffer<- mine.20km.buffer %>% 
    mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)
#plot(mine.20km.buffer)
  
  ## master 5km ring
  mine.5km.buffer.master<-st_as_sf(buffer(vect(mine.centroids.point), 5000))
  mine.5km.buffer.master<- mine.5km.buffer.master %>% 
    mutate(mine.year.10p = mine.year) %>% 
   select(CLUSTER_ID, first.mine.year, buffer.forest.prop)
  
  ## master 10km ring
   mine.10km.buffer.master<-st_as_sf(buffer(vect(mine.centroids.point), 10000))
    mine.10km.buffer.master<- mine.10km.buffer.master %>% 
     mutate(mine.year.10p = mine.year) %>% 
    select(CLUSTER_ID, first.mine.year, buffer.forest.prop)

  ## convert back to wgs84 and save to list
  # save centroid
  mine.centroids.point<- st_transform(mine.centroids.point, "EPSG:4326")
 mine.centroids.list[[i]]<- mine.centroids.point
  mine.1km.buffer<- st_transform(mine.1km.buffer, "EPSG:4326")
  mine.1km.buffer.list[[i]]<- mine.1km.buffer
  mine.2.5km.buffer<- st_transform(mine.2.5km.buffer, "EPSG:4326")
  mine.2.5km.buffer.list[[i]]<- mine.2.5km.buffer
  mine.5km.buffer<- st_transform(mine.5km.buffer, "EPSG:4326")
  mine.5km.buffer.list[[i]]<- mine.5km.buffer
  mine.10km.buffer<- st_transform(mine.10km.buffer, "EPSG:4326")
  mine.10km.buffer.list[[i]]<- mine.10km.buffer
  mine.20km.buffer<- st_transform(mine.20km.buffer, "EPSG:4326")
 mine.20km.buffer.list[[i]]<- mine.20km.buffer
  mine.5km.buffer.master<- st_transform(mine.5km.buffer.master, "EPSG:4326")
  mine.5km.buffer.master.list[[i]]<- mine.5km.buffer.master
  mine.10km.buffer.master<- st_transform(mine.10km.buffer.master, "EPSG:4326")
  mine.10km.buffer.master.list[[i]]<- mine.10km.buffer.master
}

## put together and export
#centroids
mine.centroids.df<- do.call(rbind, mine.centroids.list)
st_write(mine.centroids.df, paste0("Data/country.mining.data/centroids/", country.n, ".50p.forest.final.centroids.shp"), append = F)
#buffers
mine.1km.buffer.df<- do.call(rbind, mine.1km.buffer.list)
st_write(mine.1km.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.1km.buffer.shp"), append = F)
mine.2.5km.buffer.df<- do.call(rbind, mine.2.5km.buffer.list)
st_write(mine.2.5km.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.2.5km.buffer.shp"), append = F)
 mine.5km.buffer.df<- do.call(rbind, mine.5km.buffer.list)
st_write(mine.5km.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.5km.buffer.shp"), append = F)
mine.10km.buffer.df<- do.call(rbind, mine.10km.buffer.list)
st_write(mine.10km.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.20km.buffer.shp"), append = F)
mine.20km.buffer.df<- do.call(rbind, mine.20km.buffer.list)
st_write(mine.20km.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.20km.buffer.shp"), append = F)

mine.5km.master.buffer.df<- do.call(rbind, mine.5km.buffer.master.list)
st_write(mine.5km.master.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.5km.master.buffer.shp"), append = F)
mine.10km.master.buffer.df<- do.call(rbind, mine.10km.buffer.master.list)
st_write(mine.10km.master.buffer.df, paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.10km.master.buffer.shp"), append = F)

rm(list = ls()[ls() != "country.forest.mine.files"])
gc()
tmpFiles(remove = T)


}


#### Get forest loss in each set of buffer rings for each conutry ####


## get annual raster for each year ##
rm(list = ls())
gc()
years<- c(2001:2023)
hansen.forest.loss<- rast("Data/hansen.forest.loss/sub_saharan_Africa_forest_loss_50p_30m.tif")
for (i in 1:length(years)){
  year.filename<- years[i]
  print(year.filename)
  year<- year.filename - 2000
  # retain only pixels of forest lost in that year #
  year.forest.loss<- subst(hansen.forest.loss, year, 1, others = NA, filename = paste0("Data/hansen.forest.loss/annual/sub_saharan_Africa_forest_loss_50p_30m_", year.filename, ".tif"), overwrite = T, datatype='INT1U')
}


rm(list = ls())
gc()
tmpFiles(remove = T)
library(tidyr)
country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")

country.iso.code.df<- data.frame(
  country = c("Angola", "Burundi", "Cameroon", "Central African Republic", "Comoros", "Côte d'Ivoire","Democratic Republic of the Congo","Equatorial Guinea",
              "Ethiopia", "Gabon", "Ghana", "Guinea-Bissau", "Guinea", "Kenya", "Liberia", "Madagascar", "Malawi", "Mozambique", "Nigeria",
              "Republic of the Congo", "Rwanda", "Sierra Leone", "South Africa", "South Sudan", "Swaziland", "Tanzania", "Uganda", "Zambia", "Zimbabwe"),
  iso.code = c("AGO", "BDI", "CMR", "CAF", "COM", "CIV", "COD", "GNQ", "ETH", "GAB", "GHA", "GNB", "GIN", "KEN", "LBR", "MDG", "MWI", "MOZ", "NGA", "COG",
               "RWA", "SLE", "ZAF", "SSD", "SWZ", "TZA", "UGA", "ZMB", "ZWE")
)

for (l in 1:length(country.forest.mine.files)){
  print(l)
  country.filename<- country.forest.mine.files[l]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  # Skip if country is Comoros
  if (country.n == "Comoros") {
    next
  }## load in plantations too
  # get iso code
  country.df<- filter(country.iso.code.df, country == country.n)
  country.iso<- country.df$iso.code
  layer.n<- paste0(country.iso, "_plant_v2")
  gdb_path <- "Data/sdpt_v2_v11282023_public.gdb"
  country.plantations<- (st_read(dsn = gdb_path, layer = layer.n) %>% select(iso3, country))
  country.plantations<- vect(st_transform(country.plantations, 4326))
  # load in 2000 forest cover
  country.forest.2000<- rast(paste0("Data/hansen.forest.cover/country.forest.cover/", country.n, "_2000_tree_cover_30m.tif"))
  buffers<- c(1,5, 10, 20)
  for (i in 1:length(buffers)){
    buffer.dist<- buffers[i]
    # load in relevant buffers
    mine.buffers<- st_read(paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.", buffer.dist, "km.buffer.shp"))
    colnames(mine.buffers)[1]<- "CLUSTER_ID"
    colnames(mine.buffers)[2]<- "first.mine.year"
    colnames(mine.buffers)[3]<- "buffer.forest.prop"
    mine.buffer.forest.loss.list<- list()
    for (j in 1:nrow(mine.buffers)){
      mine<-  mine.buffers[j,]
      # add in buffer size as a variable
      mine<- mutate(mine, buffer.size = buffer.dist)
      print(paste0("-------------MINE ", j, "-----------"))
      # get number of forested cells in 2000
      mine.forest.cells.2000<- crop(country.forest.2000, vect(mine))
      mine.forest.cells.2000<- mask(mine.forest.cells.2000, vect(mine))
      # mask out plantations
      mine.plantations<- crop(country.plantations, vect(mine))
      mine.forest.cells.2000<- mask(mine.forest.cells.2000, mine.plantations, inverse = T)
      # get number of forest cells in 2000
      mine.forest.cells.number<- as.numeric(global(mine.forest.cells.2000, fun="notNA"))
      # get area of forest cells
      mine.forest.cells.area<- as.numeric(expanse(mine.forest.cells.2000, unit = "ha"))[2]
      mine<- mutate(mine, forest.cells.2000 = mine.forest.cells.number, forest.area.2000 = mine.forest.cells.area)
      # get number of forest cells lost in each year according to the jrc #
      for (k in 2001:2023){
        year<- k
        #print(year)
        # load in that year forest loss
        forest.loss.data<- rast(paste0("Data/hansen.forest.loss/annual/sub_saharan_Africa_forest_loss_50p_30m_", year, ".tif"))
        # crop by buffer
        forest.loss.in.buffer<- crop(forest.loss.data, vect(mine))
        # mask by buffer
        forest.loss.in.buffer<- mask(forest.loss.in.buffer, vect(mine))
        # mask out plantations
        forest.loss.in.buffer<- mask(forest.loss.in.buffer, mine.plantations, inverse = T)
        # get loss area as number of pixels
        # get loss area as number of pixels
        forest.loss.area<- as.numeric(freq(forest.loss.in.buffer))[3]
        # get loss as area in hectares
        forest.loss.area.hectares<- as.numeric(expanse(forest.loss.in.buffer, unit = "ha"))[2]
        forest.loss.area.df<- data.frame("forest.loss.cells" = forest.loss.area, "forest.loss.area" = forest.loss.area.hectares)
        colnames(forest.loss.area.df)[1]<- paste0("forest.loss.cells", year)
        colnames(forest.loss.area.df)[2]<- paste0("forest.loss.area", year)
        # append to mine buffer data# appforest.loss.area.dfend to mine buffer data
        mine<- cbind(mine, forest.loss.area.df)
      }
      # pivot to long format
      # Assuming your dataframe is named 'df'
     
      mine.long <- mine %>%
        pivot_longer(
          cols = starts_with("forest.loss."),
          names_to = c(".value", "loss.year"),
          names_pattern = "(forest\\.loss\\.(?:cells|area))(\\d+)"
        )
      
      
      # Convert 'loss.year' to numeric
      mine.long$loss.year <- as.numeric(mine.long$loss.year)
      # convert NAs to 0
      mine.long$forest.loss.cells[is.na(mine.long$forest.loss.cells)] <- 0
      mine.long$forest.loss.area[is.na(mine.long$forest.loss.area)] <- 0
      # get cumulative loss through the years
      mine.long<- mutate(mine.long, cumulative.forest.loss.cells = cumsum(forest.loss.cells))
      mine.long<- mutate(mine.long, cumulative.forest.loss.area = cumsum(forest.loss.area))
      # get as a proportion of original forest cells
      mine.long<- mutate(mine.long, cumulative.forest.loss.prop.cells = cumulative.forest.loss.cells/forest.cells.2000)
      mine.long<- mutate(mine.long, cumulative.forest.loss.prop.area = cumulative.forest.loss.area/forest.area.2000)
      
      # drop geometry but keep coordinates in case needed
      mine.long<- st_drop_geometry(mine.long)
      mine.buffer.forest.loss.list[[j]]<- mine.long
      tmpFiles(remove = T)
    }
    mine.buffer.forest.loss.df<- do.call(rbind, mine.buffer.forest.loss.list)
    # export as dataframe
    save(mine.buffer.forest.loss.df, file = paste0("Data/country.mining.data/forest.loss.in.buffers/", country.n, ".forest.mines.", buffer.dist, "km.buffer.forest.loss.df.RData"))
  tmpFiles(remove = T)
  }
}



### Get in master 5km ring ###

rm(list = ls())
gc()
tmpFiles(remove = T)
library(tidyr)
country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")
country.iso.code.df<- data.frame(
  country = c("Angola", "Burundi", "Cameroon", "Central African Republic", "Comoros", "Côte d'Ivoire","Democratic Republic of the Congo","Equatorial Guinea",
              "Ethiopia", "Gabon", "Ghana", "Guinea-Bissau", "Guinea", "Kenya", "Liberia", "Madagascar", "Malawi", "Mozambique", "Nigeria",
              "Republic of the Congo", "Rwanda", "Sierra Leone", "South Africa", "South Sudan", "Swaziland", "Tanzania", "Uganda", "Zambia", "Zimbabwe"),
  iso.code = c("AGO", "BDI", "CMR", "CAF", "COM", "CIV", "COD", "GNQ", "ETH", "GAB", "GHA", "GNB", "GIN", "KEN", "LBR", "MDG", "MWI", "MOZ", "NGA", "COG",
               "RWA", "SLE", "ZAF", "SSD", "SWZ", "TZA", "UGA", "ZMB", "ZWE")
)
for (l in 1:length(country.forest.mine.files)){
  print(l)
  country.filename<- country.forest.mine.files[l]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  # Skip if country is Comoros
  if (country.n == "Comoros") {
    next
  }
  ## load in plantations too
  # get iso code
  country.df<- filter(country.iso.code.df, country == country.n)
  country.iso<- country.df$iso.code
  layer.n<- paste0(country.iso, "_plant_v2")
  gdb_path <- "Data/sdpt_v2_v11282023_public.gdb"
  country.plantations<- (st_read(dsn = gdb_path, layer = layer.n) %>% select(iso3, country))
  country.plantations<- vect(st_transform(country.plantations, 4326))
  country.forest.2000<- rast(paste0("Data/hansen.forest.cover/country.forest.cover/", country.n, "_2000_tree_cover_30m.tif"))
  # load in relevant buffers
  mine.buffers<- st_read(paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.5km.master.buffer.shp"))
  colnames(mine.buffers)[1]<- "CLUSTER_ID"
  colnames(mine.buffers)[2]<- "first.mine.year"
  colnames(mine.buffers)[3]<- "buffer.forest.prop"
  mine.buffer.forest.loss.list<- list()
  for (j in 1:nrow(mine.buffers)){
    mine<-  mine.buffers[j,]
    # add in buffer size as a variable
    mine<- mutate(mine, buffer.size = "5km master")
    print(paste0("-------------MINE ", j, "-----------"))
    # get number of forested cells in 2000
    mine.forest.cells.2000<- crop(country.forest.2000, vect(mine))
    mine.forest.cells.2000<- mask(mine.forest.cells.2000, vect(mine))
    # mask out plantations
    mine.plantations<- crop(country.plantations, vect(mine))
    mine.forest.cells.2000<- mask(mine.forest.cells.2000, mine.plantations, inverse = T)
    # get number of forest cells in 2000
    mine.forest.cells.number<- as.numeric(global(mine.forest.cells.2000, fun="notNA"))
    # get area of forest cells
    mine.forest.cells.area<- as.numeric(expanse(mine.forest.cells.2000, unit = "ha"))[2]
    mine<- mutate(mine, forest.cells.2000 = mine.forest.cells.number, forest.area.2000 = mine.forest.cells.area)
    # get number of forest cells lost in each year according to the jrc #
    for (k in 2001:2023){
      year<- k
      #print(year)
      # load in that year forest loss
      forest.loss.data<- rast(paste0("Data/hansen.forest.loss/annual/sub_saharan_Africa_forest_loss_50p_30m_", year, ".tif"))
      # crop by buffer
      forest.loss.in.buffer<- crop(forest.loss.data, vect(mine))
      # mask by buffer
      forest.loss.in.buffer<- mask(forest.loss.in.buffer, vect(mine))
      # mask out plantations
      forest.loss.in.buffer<- mask(forest.loss.in.buffer, mine.plantations, inverse = T)
      # get loss area as number of pixels
      # get loss area as number of pixels
      forest.loss.area<- as.numeric(freq(forest.loss.in.buffer))[3]
      # get loss as area in hectares
      forest.loss.area.hectares<- as.numeric(expanse(forest.loss.in.buffer, unit = "ha"))[2]
      forest.loss.area.df<- data.frame("forest.loss.cells" = forest.loss.area, "forest.loss.area" = forest.loss.area.hectares)
      colnames(forest.loss.area.df)[1]<- paste0("forest.loss.cells", year)
      colnames(forest.loss.area.df)[2]<- paste0("forest.loss.area", year)
      # append to mine buffer data# appforest.loss.area.dfend to mine buffer data
      mine<- cbind(mine, forest.loss.area.df)
    }
    # pivot to long format
    # Assuming your dataframe is named 'df'
    
    mine.long <- mine %>%
      pivot_longer(
        cols = starts_with("forest.loss."),
        names_to = c(".value", "loss.year"),
        names_pattern = "(forest\\.loss\\.(?:cells|area))(\\d+)"
      )
    
    
    # Convert 'loss.year' to numeric
    mine.long$loss.year <- as.numeric(mine.long$loss.year)
    # convert NAs to 0
    mine.long$forest.loss.cells[is.na(mine.long$forest.loss.cells)] <- 0
    mine.long$forest.loss.area[is.na(mine.long$forest.loss.area)] <- 0
    # get cumulative loss through the years
    mine.long<- mutate(mine.long, cumulative.forest.loss.cells = cumsum(forest.loss.cells))
    mine.long<- mutate(mine.long, cumulative.forest.loss.area = cumsum(forest.loss.area))
    # get as a proportion of original forest cells
    mine.long<- mutate(mine.long, cumulative.forest.loss.prop.cells = cumulative.forest.loss.cells/forest.cells.2000)
    mine.long<- mutate(mine.long, cumulative.forest.loss.prop.area = cumulative.forest.loss.area/forest.area.2000)
    
    # drop geometry but keep coordinates in case needed
    mine.long<- st_drop_geometry(mine.long)
    mine.buffer.forest.loss.list[[j]]<- mine.long
    tmpFiles(remove = T)
  }
  mine.buffer.forest.loss.df<- do.call(rbind, mine.buffer.forest.loss.list)
  # export as dataframe
  save(mine.buffer.forest.loss.df, file = paste0("Data/country.mining.data/forest.loss.in.buffers/", country.n, ".forest.mines.5km.master.buffer.forest.loss.df.RData"))
  tmpFiles(remove = T)
}


#### Get mining only forest loss per year ####


## get mine loss year only for 50% tree cover ##
africa.tree.cover<- rast(paste0("Data/hansen.forest.cover/sub_saharan_Africa_2000_tree_cover_30m_resamlped.tif"))

# load in mining data #
mining.loss.year<- rast("Data/masolele.processed/africa.mining.loss.year.tif")
#resample to match hansen forest loss
mining.loss.year.resampled<- project(mining.loss.year, africa.tree.cover, filename = "Data/masolele.processed/africa.mining.loss.year.resampled.tif", overwrite = T)
# mask by 50p forest cover
mining.loss.year.50p<- mask(mining.loss.year.resampled, africa.tree.cover, filename = "Data/masolele.processed/africa.mining.loss.year.50p.tree.cover.tif", overwrite = T)

## split up into annual rasters ##
for (i in 1:20){
  year<- i
  year.filesave <- i + 2000
  print(year.filesave)
  mining.loss.annual<- subst(mining.loss.year.50p, i, 1, others = NA, filename = paste0("Data/masolele.processed/africa.mining.loss.year.50p.tree.cover.", year.filesave, ".tif"))
}


rm(list = ls())
gc()
tmpFiles(remove = T)
library(tidyr)
country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")
country.iso.code.df<- data.frame(
  country = c("Angola", "Burundi", "Cameroon", "Central African Republic", "Comoros", "Côte d'Ivoire","Democratic Republic of the Congo","Equatorial Guinea",
              "Ethiopia", "Gabon", "Ghana", "Guinea-Bissau", "Guinea", "Kenya", "Liberia", "Madagascar", "Malawi", "Mozambique", "Nigeria",
              "Republic of the Congo", "Rwanda", "Sierra Leone", "South Africa", "South Sudan", "Swaziland", "Tanzania", "Uganda", "Zambia", "Zimbabwe"),
  iso.code = c("AGO", "BDI", "CMR", "CAF", "COM", "CIV", "COD", "GNQ", "ETH", "GAB", "GHA", "GNB", "GIN", "KEN", "LBR", "MDG", "MWI", "MOZ", "NGA", "COG",
               "RWA", "SLE", "ZAF", "SSD", "SWZ", "TZA", "UGA", "ZMB", "ZWE")
)
mining.loss.year.50p<- rast("Data/masolele.processed/africa.mining.loss.year.50p.tree.cover.tif")
gc()
for (l in 1:length(country.forest.mine.files)){
  print(l)
  country.filename<- country.forest.mine.files[l]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  # Skip if country is Comoros
  if (country.n == "Comoros") {
    next
  }
  ## load in plantations too
  # get iso code
  country.df<- filter(country.iso.code.df, country == country.n)
  country.iso<- country.df$iso.code
  layer.n<- paste0(country.iso, "_plant_v2")
  gdb_path <- "Data/sdpt_v2_v11282023_public.gdb"
  country.plantations<- (st_read(dsn = gdb_path, layer = layer.n) %>% select(iso3, country))
  country.plantations<- vect(st_transform(country.plantations, 4326))
  country.forest.2000<- rast(paste0("Data/hansen.forest.cover/country.forest.cover/", country.n, "_2000_tree_cover_30m.tif"))
#  buffers<- c(1,5,10)
#  for (i in 1:length(buffers)){
#    buffer.dist<- buffers[i]
    # load in relevant buffers
    mine.buffers<- st_read(paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.5km.master.buffer.shp"))
    colnames(mine.buffers)[1]<- "CLUSTER_ID"
    colnames(mine.buffers)[2]<- "first.mine.year"
    colnames(mine.buffers)[3]<- "buffer.forest.prop"
    mine.buffer.forest.loss.list<- list()
    for (j in 1:nrow(mine.buffers)){
      mine<-  mine.buffers[j,]
      # add in buffer size as a variable
      mine<- mutate(mine, buffer.size = "5km Master")
      print(paste0("-------------MINE ", j, "-----------"))
      # get number of forested cells in 2000
      mine.forest.cells.2000<- crop(country.forest.2000, vect(mine))
      mine.forest.cells.2000<- mask(mine.forest.cells.2000, vect(mine))
      # mask out plantations
      mine.plantations<- crop(country.plantations, vect(mine))
      mine.forest.cells.2000<- mask(mine.forest.cells.2000, mine.plantations, inverse = T)
      # get number of forest cells in 2000
      mine.forest.cells.number<- as.numeric(global(mine.forest.cells.2000, fun="notNA"))
      # get area of forest cells
      mine.forest.cells.area<- as.numeric(expanse(mine.forest.cells.2000, unit = "ha"))[2]
      mine<- mutate(mine, forest.cells.2000 = mine.forest.cells.number, forest.area.2000 = mine.forest.cells.area)
      
      for (k in 2001:2020){
        year<- k
        #print(year)
        # load in that year forest loss
        forest.loss.data<- rast(paste0("Data/masolele.processed/africa.mining.loss.year.50p.tree.cover.", year, ".tif"))
        # crop by buffer
        forest.loss.in.buffer<- crop(forest.loss.data, vect(mine))
        # mask by buffer
        forest.loss.in.buffer<- mask(forest.loss.in.buffer, vect(mine))
    #  plot(forest.loss.in.buffer)
        # mask out plantations
        forest.loss.in.buffer<- mask(forest.loss.in.buffer, mine.plantations, inverse = T)
        # get loss area as number of pixels
        forest.loss.area<- as.numeric(freq(forest.loss.in.buffer))[3]
        # get loss as area in hectares
        forest.loss.area.hectares<- as.numeric(expanse(forest.loss.in.buffer, unit = "ha"))[2]
        forest.loss.area.df<- data.frame("forest.loss.cells" = forest.loss.area, "forest.loss.area" = forest.loss.area.hectares)
        colnames(forest.loss.area.df)[1]<- paste0("forest.loss.cells", year)
        colnames(forest.loss.area.df)[2]<- paste0("forest.loss.area", year)
        
        # append to mine buffer data# appforest.loss.area.dfend to mine buffer data
        mine<- cbind(mine, forest.loss.area.df)
      }
      # pivot to long format
      # Assuming your dataframe is named 'df'
      mine.long <- mine %>%
        pivot_longer(
          cols = starts_with("forest.loss."),
          names_to = c(".value", "loss.year"),
          names_pattern = "(forest\\.loss\\.(?:cells|area))(\\d+)"
        )
      
      
      # Convert 'loss.year' to numeric
      mine.long$loss.year <- as.numeric(mine.long$loss.year)
      # convert NAs to 0
      mine.long$forest.loss.cells[is.na(mine.long$forest.loss.cells)] <- 0
      mine.long$forest.loss.area[is.na(mine.long$forest.loss.area)] <- 0
      # get cumulative loss through the years
      mine.long<- mutate(mine.long, cumulative.forest.loss.cells = cumsum(forest.loss.cells))
      mine.long<- mutate(mine.long, cumulative.forest.loss.area = cumsum(forest.loss.area))
      # get as a proportion of original forest cells
      mine.long<- mutate(mine.long, cumulative.forest.loss.prop.cells = cumulative.forest.loss.cells/forest.cells.2000)
      mine.long<- mutate(mine.long, cumulative.forest.loss.prop.area = cumulative.forest.loss.area/forest.area.2000)
      # drop geometry but keep coordinates in case needed
      mine.long<- st_drop_geometry(mine.long)
      mine.buffer.forest.loss.list[[j]]<- mine.long
      tmpFiles(remove = T)
    }
    mine.buffer.forest.loss.df<- do.call(rbind, mine.buffer.forest.loss.list)
    # export as dataframe
    save(mine.buffer.forest.loss.df, file = paste0("Data/country.mining.data/forest.loss.in.buffers/", country.n, ".forest.mines.5km.master.buffer.forest.loss.mining.only.df.RData"))
    tmpFiles(remove = T)
  }
}




#### Get mine covariates for DiD modelling ####


country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")
# Load the rasters
travel_time <- rast("Data/covariates/travel_time_to_cities_12.tif")
elevation<- rast("Data/covariates/africa.elevation.nasadem.100m.tif")
slope<- rast("Data/covariates/africa.slope.nasadem.100m.tif")
population_density <- rast("Data/covariates/gpw.2000.2020.mean.tif")
biomes<- st_read("Data/covariates/biome.numbered.ecoregions.shp")
biomes<- st_transform(biomes, "ESRI:54009")

for (l in 1:length(country.forest.mine.files)){
  print(l)
  country.filename<- country.forest.mine.files[l]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  
  ## load in mine centroids ##
  mine.cluster.points<- st_read(paste0("Data/country.mining.data/centroids/", country.n, ".50p.forest.final.centroids.shp"))
  colnames(mine.cluster.points)[1]<- "CLUSTER_ID"
  colnames(mine.cluster.points)[2]<- "first.mine.year"
  colnames(mine.cluster.points)[3]<- "buffer.forest.prop"
  # get covariates for each point #
  
  ## Travel time ##
  
  # extract values for sample #
  extracted.travel.time<- terra::extract(travel_time,mine.cluster.points, cells=FALSE, method="simple") # this command takes the values from 'travel_time' raster that correspond to the cells that each point in 'mine cluster points' falls into
  colnames(extracted.travel.time)[2]<- "travel.time"
  # bind with other data
  mine.cluster.points<- cbind(mine.cluster.points, extracted.travel.time[2])
  
  ## Population density ##
  
  # extract values for sample #
  extracted.pop.density<- terra::extract(population_density,mine.cluster.points, cells=FALSE, method="simple")
  colnames(extracted.pop.density)[2]<- "pop.density"
  # bind with other data
  mine.cluster.points<- cbind(mine.cluster.points, extracted.pop.density[2])
  
  ## Elevation ##
  
  # extract values for sample #
  extracted.elevation<- terra::extract(elevation,mine.cluster.points, cells=FALSE, method="simple")
  colnames(extracted.elevation)[2]<- "elevation"
  # bind with other data
  mine.cluster.points<- cbind(mine.cluster.points, extracted.elevation[2])
  
  ## slope ##
  
  # extract values for sample #
  extracted.slope<- terra::extract(slope,mine.cluster.points, cells=FALSE, method="simple")
  colnames(extracted.slope)[2]<- "slope"
  # bind with other data
  mine.cluster.points<- cbind(mine.cluster.points, extracted.slope[2])
  
  
  ## biome ##
  mine.cluster.points.reprojected<- st_transform(mine.cluster.points, "ESRI:54009")
  biomes.country<- st_crop(biomes,mine.cluster.points.reprojected )
  # extract biome for sample
  mine.cluster.points$biome<- as.factor(st_join(mine.cluster.points.reprojected, biomes.country, join = st_within) %>% pull(BIOME_NUM))
  ## country ##
  mine.cluster.points<- mutate(mine.cluster.points, country = country.n)
  
  # remove geometry
  mine.cluster.points.df<- st_drop_geometry(mine.cluster.points)
  save(mine.cluster.points.df, file = paste0("Data/country.mining.data/model.covariates/", country.n, ".50p.forest.final.COVARIATES.RData"))
}

### Get first year of spillover with another mine ###
country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")

## this final part is for the spillover supplementary analysis, it just works out if a mine overlaps with the 5km buffer of another mine, and if so, what year does that mine start in

## 10km master buffer spillover ## 

for (l in 1:length(country.forest.mine.files)){
  print(l)
  country.filename<- country.forest.mine.files[l]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  
  load(file = paste0("Data/country.mining.data/model.covariates/", country.n, ".50p.forest.final.COVARIATES.RData"))
  mine.cluster.points.df<- select(mine.cluster.points.df, -spillover.10km.master)
  mine.buffers.10km<- st_read(paste0("Data/country.mining.data/buffers/", country.n, ".50p.forest.final.10km.master.buffer.shp"))
  mine.buffers.10km<- st_transform(mine.buffers.10km, "ESRI:54009")
  colnames(mine.buffers.10km)[1]<- "CLUSTER_ID"
  colnames(mine.buffers.10km)[2]<- "first.mine.year"
  colnames(mine.buffers.10km)[3]<- "buffer.forest.prop"
  mine.buffers.10km<- left_join(mine.buffers.10km, mine.cluster.points.df)
  spillover.data.list<- list()
  mine.clusters<- st_read(paste0("Data/country.mining.data/country.forest.mine.shapefiles/", country.n, ".forest.mining.clusters.shp"))
  mine.clusters<- st_transform(mine.clusters, "ESRI:54009")
  colnames(mine.clusters)[1]<- "CLUSTER_ID"
  colnames(mine.clusters)[2]<- "first.mine.year"
  colnames(mine.clusters)[3]<- "buffer.forest.prop"
  mine.clusters<- left_join(mine.clusters, mine.cluster.points.df)
  for (i in 1:nrow(mine.buffers.10km)){
    print(i)
    mine.buffer<- mine.buffers.10km[i,]
    mine.buffer.n<- mine.buffer$CLUSTER_ID
    # get mines from other sites
    other.mines<- filter(mine.clusters, CLUSTER_ID != mine.buffer.n)
    # filter for other mines that overlap with this buffer
    overlapped.mine.buffers<-  other.mines %>%
      st_filter(y=mine.buffer, .predicate = st_intersects)
    # if overlapping mines, calculate the earliest instance of overlap, if no overlap, set as NA
    if (nrow(overlapped.mine.buffers) > 0){
      first.overlap.year<- overlapped.mine.buffers %>% 
        st_drop_geometry() %>% 
        arrange(year.10p) %>% 
        head(1)
      first.overlap.year<- first.overlap.year$year.10p
    } else {first.overlap.year<- NA}
    # combine with covariate data
    spillover.data.list[[i]]<- data.frame(CLUSTER_ID = mine.buffer.n, spillover.10km.master = first.overlap.year)
  }
  
  ## put back in with mine data
  spillover.data.df<- do.call(rbind, spillover.data.list)
  mine.cluster.points.df<- left_join(mine.cluster.points.df, spillover.data.df)
  
  
  save(mine.cluster.points.df, file = paste0("Data/country.mining.data/model.covariates/", country.n, ".50p.forest.final.COVARIATES.RData"))
}

  
  for (l in 1:length(country.forest.mine.files)){
    print(l)
    country.filename<- country.forest.mine.files[l]
    # get country name
    country.n<- str_extract(country.filename, "^[^\\.]+")
    print(country.n)
    
    ## load in mine clusters ##
    mine.shapefiles<- st_read(paste0("Data/country.mining.data/country.forest.mine.shapefiles/", country.n, ".forest.mining.clusters.shp"))
    colnames(mine.shapefiles)[1]<- "CLUSTER_ID"
    colnames(mine.shapefiles)[2]<- "first.mine.year"
    colnames(mine.shapefiles)[3]<- "buffer.forest.prop"
    ## load in forest loss points ##
    all.mines.forest.loss.points<- st_read(paste0("Data/country.mining.data/points/", country.n, ".mining.loss.year.points.shp"))
    colnames(all.mines.forest.loss.points)[1]<- "year"
    new.mine.year.list<- list()
     for (i in 1:nrow(mine.shapefiles)){
   print(i)
  mine<- mine.shapefiles[i,]
    # get forest loss points for mine
    mine.forest.loss.points<- st_intersection(mine, all.mines.forest.loss.points)
    mine.forest.loss.points<- st_drop_geometry(mine.forest.loss.points)
    total.points<- nrow(mine.forest.loss.points)
    # get summary
    mine.forest.loss.year.summary<- mine.forest.loss.points %>% 
      group_by(year) %>% 
      summarise(num.points = n()) %>% 
      mutate(cum.points= cumsum(num.points), cum.prop = cum.points / total.points)
    # get different props
    year.10p<- mine.forest.loss.year.summary %>% 
      filter(cum.prop >= 0.1) %>% 
      arrange(year) %>% 
      head(1)
    year.10p<- year.10p$year
    year.20p<- mine.forest.loss.year.summary %>% 
      filter(cum.prop >= 0.2) %>% 
      arrange(year) %>% 
      head(1)
    year.20p<- year.20p$year
    year.25p<- mine.forest.loss.year.summary %>% 
      filter(cum.prop >= 0.25) %>% 
      arrange(year) %>% 
      head(1)
    year.25p<- year.25p$year
    new.mine.year.df<- data.frame(CLUSTER_ID = mine$CLUSTER, "year.10p"=year.10p,"year.20p"=year.20p, "year.25p"=year.25p )
    new.mine.year.list[[i]]<- new.mine.year.df
    }
 new.mine.year.df<- do.call(rbind, new.mine.year.list)
load(file = paste0("Data/country.mining.data/model.covariates/", country.n, ".50p.forest.final.COVARIATES.RData"))
mine.cluster.points.df<- left_join(mine.cluster.points.df, new.mine.year.df)
save(mine.cluster.points.df, file = paste0("Data/country.mining.data/model.covariates/", country.n, ".50p.forest.final.COVARIATES.RData"))
}

#### Get commodities for mines from Maus ####

country.forest.mine.files<- list.files(path = "Data/country.mining.data/country.forest.mine.shapefiles/", pattern = ".forest.mining.clusters.shp")
maus.mine.data<- st_read("Data/maus.mining/mine_polygons.gpkg")
maus.mine.data.moll<- st_transform(maus.mine.data, "ESRI:54009")
for (j in 1:length(country.forest.mine.files)){
  print(j)
  country.filename<- country.forest.mine.files[j]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  
  mine.clusters<- st_read(paste0("Data/country.mining.data/country.forest.mine.shapefiles/", country.filename))
  mine.clusters<- st_transform(mine.clusters, "ESRI:54009")
  print(paste0(country.n, " - ", nrow(mine.clusters), " Mines"))
  colnames(mine.clusters)[1]<- "CLUSTER_ID"
  colnames(mine.clusters)[2]<- "first.mine.year"
  colnames(mine.clusters)[3]<- "buffer.forest.prop"
  country.mining.commodity.data.list<- list()
  for (i in 1:nrow(mine.clusters)){
    mine.cluster<- mine.clusters[i,]
    print(paste0("-------------MINE ", i, "-----------"))
    # crop mines down to cluster
    maus.mines.in.cluster<- st_crop(maus.mine.data.moll, mine.cluster)
    if (nrow(maus.mines.in.cluster) > 0) {
      mine.cluster.intersect<- st_intersection(maus.mines.in.cluster)
      # get mining material
      mine.cluster.material<- mine.cluster.intersect$materials_list
      mine.cluster.material<- ifelse(is.na(mine.cluster.material), "Unknown", as.character(mine.cluster.material))
    } else {mine.cluster.material<- "No mine in Maus"}
    country.mining.commodity.data.list[[i]]<- data.frame(CLUSTER_ID = mine.cluster$CLUSTER_ID, country = country.n, material = mine.cluster.material)
  }
  country.mining.commodity.data.df<- do.call(rbind, country.mining.commodity.data.list)
  country.mining.commodity.data.df<- distinct(country.mining.commodity.data.df)
  print(country.n)
  print(country.mining.commodity.data.df %>% group_by(material) %>% summarise(n.mines = n()))
  save(country.mining.commodity.data.df, file = paste0("Data/country.mining.data/country.mining.commodity.type/", country.n, ".mine.clusters.by.commodity.secondary.RData"))
}  

## Commodities within 5km buffer ##

country.forest.mine.files<- list.files(path = "Data/country.mining.data/buffers/", pattern = ".5km.master.buffer.shp")
maus.mine.data<- st_read("Data/maus.mining/mine_polygons.gpkg")
maus.mine.data.moll<- st_transform(maus.mine.data, "ESRI:54009")
for (j in 1:length(country.forest.mine.files)){
  print(j)
  country.filename<- country.forest.mine.files[j]
  # get country name
  country.n<- str_extract(country.filename, "^[^\\.]+")
  print(country.n)
  
  mine.clusters<- st_read(paste0("Data/country.mining.data/buffers/", country.filename))
  mine.clusters<- st_transform(mine.clusters, "ESRI:54009")
  print(paste0(country.n, " - ", nrow(mine.clusters), " Mines"))
  colnames(mine.clusters)[1]<- "CLUSTER_ID"
  colnames(mine.clusters)[2]<- "first.mine.year"
  colnames(mine.clusters)[3]<- "buffer.forest.prop"
  country.mining.commodity.data.list<- list()
  for (i in 1:nrow(mine.clusters)){
    mine.cluster<- mine.clusters[i,]
    print(paste0("-------------MINE ", i, "-----------"))
    # crop mines down to cluster
    maus.mines.in.cluster<- st_crop(maus.mine.data.moll, mine.cluster)
    if (nrow(maus.mines.in.cluster) > 0) {
      mine.cluster.intersect<- st_intersection(maus.mines.in.cluster)
      # get mining material
      mine.cluster.material<- mine.cluster.intersect$materials_list
      mine.cluster.material<- ifelse(is.na(mine.cluster.material), "Unknown", as.character(mine.cluster.material))
    } else {mine.cluster.material<- "No mine in Maus"}
    country.mining.commodity.data.list[[i]]<- data.frame(CLUSTER_ID = mine.cluster$CLUSTER_ID, country = country.n, material = mine.cluster.material)
  }
  country.mining.commodity.data.df<- do.call(rbind, country.mining.commodity.data.list)
  country.mining.commodity.data.df<- distinct(country.mining.commodity.data.df)
  print(country.n)
  print(country.mining.commodity.data.df %>% group_by(material) %>% summarise(n.mines = n()))
  save(country.mining.commodity.data.df, file = paste0("Data/country.mining.data/country.mining.commodity.type/", country.n, ".mine.clusters.by.commodity.5km.buffer.secondary.RData"))
}  

