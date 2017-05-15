# reprojection
library(raster)
library(rgdal)
library(maptools)
library(BAMMtools)
library(rasterVis)
library(latticeExtra)
library(gridExtra)
library(grid)

setwd("/Volumes/dongmeic-10/fire/revision/maps")
infolder <- "/Volumes/dongmeic-10/fire/output/revision/results/"
outfolder <- "/Volumes/dongmeic-10/fire/revision/maps/"
ch.aea <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=30 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
pro.shpfile <- "/Volumes/dongmeic-10/fire/mapchina/China_province_re.shp"
pro.shp <- readShapePoly(pro.shpfile, proj4string = ch.aea)
ch.pro <- spTransform(pro.shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# color <- 'YlOrRd'

## setting for reprojection
# arrow <- list("SpatialPolygonsRescale", layout.north.arrow(type=1), 
#               offset = c(2100000, 500000), scale = 550000)
# scale <- list("SpatialPolygonsRescale", layout.scale.bar(height = 0.15), 
#                offset = c(-2200000, -1100000), scale = 600000, fill=c("transparent","black"))
# text <- list("sp.text", c(-1900000, -810000), "550 km", cex=1.2)

## setting for WGS84
scale <- list("SpatialPolygonsRescale", layout.scale.bar(height = 0.15), 
              offset = c(84, 22), scale = 5, fill=c("transparent","black"))

text <- list("sp.text", c(86.5, 23.5), "550 km", cex=1.2)

arrow <- list("SpatialPolygonsRescale", layout.north.arrow(type=1),
              offset = c(127, 32), scale = 4)

bound <- list("sp.polygons", ch.pro, col="grey", lwd=0.8, first = FALSE)

filenames <- c('mean_annual_number_of_fires.tif','mean_annual_active_fire.tif', 
               'mean_annual_number_of_fires_cv.tif','mean_annual_active_fire_cv.tif', 
              'fire_season_duration.tif', 'peak_month.tif', 'fire_radiative_power.tif', 
              'Gini_index.tif', 'percentage_forest_affected.tif', 
              'percentage_savanna_affected.tif', 'percentage_grassland_affected.tif', 
              'percentage_cropland_affected.tif')
names0 <- c('Mean Annual Number of Fires','Mean Annual active Fire Density',
            'CV in Number of Fires', 'CV in active Fire Density', 
           'Fire Season Duration', 'Fire Peak Month', 'Fire Radiative Power', 
           'Gini Index', 'Percent Forest Affected', 'Percent Savanna Affected', 
           'Percent Grassland Affected', 'Percent Cropland Affected')
#names1 <- c('MANF', 'MAFD', 'CVNF', 'CVFD', 'FSD', 'FPM', 'FRP', 'GI', 'PFA', 'PSA', 'PGA', 'PCA')
#rgtypes <-c('LHLHHf', 'LLLLLfc', 'LHMHHg', 'HLHLLfs', 'HLMLHc', 'LLLHLg')

rstack <- stack()
for (i in 1:length(filenames)){
  r.file <- paste(infolder, filenames[i], sep="")
  r0 <- raster(r.file)
  ## reprojection
  # if (i==6){
  #   r1 <- projectRaster(r0, crs=ch.aea, method="ngb")
  # }else{
  #   r1 <- projectRaster(r0, crs=ch.aea)
  # }
  rstack <- stack(rstack, r0)
  print(paste(filenames[i], "added!"))
}
names(rstack) <- names0
#mapTheme <- rasterTheme(region=brewer.pal(6, color))
mapTheme <- rasterTheme(region=rev(heat.colors(32)))

## test multiple plots
# fun <- function(){
#   plot(pro.shp, add=TRUE, bord='grey', lwd=0.8)
# }
# plot(subset(rstack, 1:4), addfun=fun, nc=2, nr=2)
# levelplot(subset(rstack, 1:4), layout=c(2,2))

## test single plot
# layer1 <- subset(rstack, 1)
# layer1[layer1 <= 0] <- NA; 
# plot(layer1, col = rev(heat.colors(32)), zlim = c(cellStats(layer1,stat='min'), cellStats(layer1,stat='max')))
# spplot(layer1, draw = T, colorkey=list(space="right", height=0.5), cuts=20,
#        scales=list(draw=T), sp.layout=list(bound,text,scale,arrow))

# map variables
text <- c("a.", "b.", "c.", "d.")
for (i in seq(1, 12, by=4)){
  layer1 <- subset(rstack, i); layer1[layer1 <= 0] <- NA; 
  layer2 <- subset(rstack, i+1); layer2[layer2 <= 0] <- NA
  layer3 <- subset(rstack, i+2); layer3[layer3 <= 0] <- NA
  layer4 <- subset(rstack, i+3); layer4[layer4 <= 0] <- NA
  p1 <- levelplot(layer1, at=getJenksBreaks(na.omit(getValues(layer1)), 6), scales=list(draw=FALSE), margin=F, par.settings=mapTheme, main=paste(text[1],names0[i]))+latticeExtra::layer(sp.polygons(ch.pro, lwd=0.8, col='gray'))
  if (i == 5){
    p2 <- levelplot(layer2, at=seq(1,12,by=1), scales=list(draw=FALSE), margin=F, par.settings=mapTheme, main=paste(text[2],names0[i+1]))+latticeExtra::layer(sp.polygons(ch.pro, lwd=0.8, col='gray'))
  }else{
    p2 <- levelplot(layer2, at=getJenksBreaks(na.omit(getValues(layer2)), 6), scales=list(draw=FALSE), margin=F, par.settings=mapTheme, main=paste(text[2],names0[i+1]))+latticeExtra::layer(sp.polygons(ch.pro, lwd=0.8, col='gray')) 
  }
  p3 <- levelplot(layer3, at=getJenksBreaks(na.omit(getValues(layer3)), 6), scales=list(draw=FALSE), margin=F, par.settings=mapTheme, main=paste(text[3],names0[i+2]))+latticeExtra::layer(sp.polygons(ch.pro, lwd=0.8, col='gray'))
  p4 <- levelplot(layer4, at=getJenksBreaks(na.omit(getValues(layer4)), 6), scales=list(draw=FALSE), margin=F, par.settings=mapTheme, main=paste(text[4],names0[i+3]))+latticeExtra::layer(sp.polygons(ch.pro, lwd=0.8, col='gray'))
  png(paste(outfolder, "variable_map_", i, ".png", sep=""), width=12, height=8, units="in", res=300)
  grid.arrange(p1, p2, p3, p4, ncol=2)
  dev.off()
  print(paste("variable map", i, "done!"))
}

# map fire regimes
inpath <- "/Volumes/dongmeic-10/fire/revision/maps/"
# regimes <- stack(paste(inpath, "regime_14-May-2017.tif", sep=""))
# regimes <- projectRaster(regimes, method="ngb", crs=ch.aea, alignOnly=T)
# plotRGB(regimes, r = 1, g = 2, b = 3, stretch="lin", addfun=fun)

# for single bands
plot(regimes)
# c("Fire regime 5", "Fire regime 6", "Fire regime 3", "Fire regime 4", "Fire regime 1", "Fire regime 2")
mycolors1 <- c('#e41a1c', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#377eb8')
mycolors2 <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33') # red, blue, green, purple, orange, yellow
mykey <- list(text=list(lab=c("Fire regime 1", "Fire regime 2", 
                              "Fire regime 3", "Fire regime 4", "Fire regime 5", 
                              "Fire regime 6"), cex=c(1.2,1.2,1.2,1.2,1.2,1.2,1.2)), 
              rectangles=list(col=mycolors2), space="inside", width=0.2, columns=1)
# gamma = 2.2
# Y = .2126 * R^gamma + .7152 * G^gamma + .0722 * B^gamma

# bands
# rgb2gray <- function(r){
#   band1 <- raster(paste(inpath, r, sep=""), band=1); 
#   band2 <- raster(paste(inpath, r, sep=""), band=2); 
#   band3 <- raster(paste(inpath, r, sep=""), band=3);
#   r.g <- band1*0.2126 + 0.7152*band2 + 0.0722*band3
#   r.g[r.g==1] <- NA
#   r.g <- projectRaster(r.g, crs=ch.aea, over=T)
#   return(r.g)
# }

#regimes.g <- rgb2gray("regime_14-May-2017.tif")
regimes <- raster(paste(inpath,"regime_14-May-2017.tif", sep=""))
regimes[regimes==0] <- NA
table(getValues(regimes))
#regimes <- projectRaster(regimes, method="ngb", crs=ch.aea)
#table(getValues(regimes))
png("fire_regimes.png", width=9, height=6, units="in", res=300)
p <- spplot(regimes, draw = T, col.regions=mycolors1, colorkey=FALSE, key=mykey, cuts=5,
            scales=list(draw=F), sp.layout=list(bound,text,scale,arrow))
names(p$legend) <- "inside"
p$legend$inside$x <- 0.35
p$legend$inside$y <- 0.97
plot(p)
dev.off()

# rgb.palette <- colorRampPalette(mycolors, space = "rgb")
# 
# spplot(regimes, col.regions=rgb.palette,
#        colorkey=list(height=0.3),
#        sp.layout=list(text,scale,arrow,bound))

# reproject rasters; should write in RGB
# for (i in 1:7){
#   filename <- paste(inpath, "regime_", i, ".tif", sep="")
#   r <- stack(filename)
#   r <- projectRaster(r, crs=ch.aea, over=T)
#   writeRaster(r, filename, "GTiff", overwrite=T)
#   print(paste(filename, "done!"))
# }

# individual fire regimes
read.raster <- function(r.file){
  r <- raster(paste(inpath,r.file, sep=""))
  r[r==0] <- NA
  #r <- projectRaster(r, method="ngb", crs=ch.aea)
  return(r)
  #plot(r)
}
regime1 <- read.raster("regime_1.tif")
p1 <- spplot(regime1, draw = T, col.regions=mycolors2[1], colorkey=FALSE, main=paste("Fire regime 1:", rgtypes[1]),
            scales=list(draw=F), sp.layout=list(bound))
plot(p1)

regime2 <- read.raster("regime_2.tif")
p2 <- spplot(regime2, draw = T, col.regions=mycolors2[2], colorkey=FALSE, main=paste("Fire regime 2:", rgtypes[2]),
             scales=list(draw=F), sp.layout=list(bound))
plot(p2)

regime3 <- read.raster("regime_3.tif")
p3 <- spplot(regime3, draw = T, col.regions=mycolors2[3], colorkey=FALSE, main=paste("Fire regime 3:", rgtypes[3]),
             scales=list(draw=F), sp.layout=list(bound))
plot(p3)

regime4 <- read.raster("regime_4.tif")
p4 <- spplot(regime4, draw = T, col.regions=mycolors2[4], colorkey=FALSE, main=paste("Fire regime 4:", rgtypes[4]),
             scales=list(draw=F), sp.layout=list(bound))
plot(p4)

regime5 <- read.raster("regime_5.tif")
p5 <- spplot(regime5, draw = T, col.regions=mycolors2[5], colorkey=FALSE, main=paste("Fire regime 5:", rgtypes[5]),
             scales=list(draw=F), sp.layout=list(bound))
plot(p5)

regime6 <- read.raster("regime_6.tif")
p6 <- spplot(regime6, draw = T, col.regions=mycolors2[6], colorkey=FALSE, main=paste("Fire regime 6:", rgtypes[6]),
             scales=list(draw=F), sp.layout=list(bound))
plot(p6)

rgtypes <-c('LLLLLfc', 'LHMHHg', 'HLHLLfs', 'HLMLHc', 'LHLHHf', 'LLLHLg')
png("fire_regimes_all.png", width=10, height=6, units="in", res=300)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)
dev.off()

# make a map of China - better do it in ArcMap
mtains <- readShapeSpatial("/Volumes/dongmeic-10/fire/revision/mountains.shp")
proj4string(mtains) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
LC.China <- raster("/Volumes/dongmeic-10/fire/output/revision/LC_China/reclass/LC_China_recl_2013.tif")
mask <- raster("/Volumes/dongmeic-10/fire/masknew.tif")
LC.China <- mask(LC.China, mask)
LC.China[LC.China==0] <- NA
mykey.1 <- list(text=list(lab=c("Forests", "Savannas", "Grasslands", "Croplands"), cex=c(1.2,1.2,1.2,1.2,1.2,1.2,1.2)), 
              rectangles=list(col=terrain.colors(4)), space="inside", width=0.2, columns=1)

png("land_cover.png", width=9, height=6, units="in", res=300)
p.1 <- spplot(LC.China, draw = T, col.regions= terrain.colors(4, alpha = 1), colorkey=FALSE, key=mykey.1, cuts=3,
             scales=list(draw=F), sp.layout=list(bound,text,scale,arrow))
names(p.1$legend) <- "inside"
p.1$legend$inside$x <- 0.37
p.1$legend$inside$y <- 0.9
plot(p.1)
dev.off()