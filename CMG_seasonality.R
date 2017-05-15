library(gdalUtils)
library(raster)
library(rgdal)
library(maptools)
library(rasterVis)
library(classInt)

inpath1 <- "/Users/dongmeichen/CMG/output/CloudCorrFirePix/"
inpath0 <- "/Users/dongmeichen/CMG/output/CorrFirePix/"
outpath <- "/Volumes/dongmeic-10/fire/revision/maps/"
mask <- raster("/Volumes/dongmeic-10/fire/mask_new.tif")
mask[mask < 0] = 0

ch.shpfile <- "/Volumes/dongmeic-10/fire/chinabord/boundarynew.shp"
ch.shp <- readShapePoly(ch.shpfile, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
bound <- list("sp.polygons", ch.shp, col="grey", lwd=0.8, first = FALSE)

s <- stack()
for (month in 1:12){
  s.m <- stack()
  for (year in 2003:2016){
    r.file1 <- paste(inpath1, "CMG_Global_", year, "-", month, ".tif", sep="")
    r.file0 <- paste(inpath0, "CMG_Global_", year, "-", month, ".tif", sep="")
    r1 <- raster(r.file1)
    r0 <- raster(r.file0)
    rc1 <- crop(r1, extent(mask))
    rc0 <- crop(r0, extent(mask))
    rc2 <- rc1 - rc0
    rc2 <- mask(rc2, mask)
    rc2[rc2==0] <- NA
    s.m <- stack(s.m, rc2)
    print(paste(r.file1, year, "done!"))
  }
  print(paste("month", month, "processed!"))
  s <- stack(s, s.m)
}

ss <- stack()
for (i in 1:12){
  r <- calc(subset(s, (i*14-13):(i*14)), sum)
  ss <- stack(ss, r)
  print(paste("month", i, "processed!"))
} 
names(ss) <- c("January", "Feburary", "March", "April", "May", "June", "July", "August", 
               "September", "October", "November", "December")
png(paste(outpath,"cmg_seasonality.png",sep=""), width=10, height=6, units="in", res=300)
spplot(ss, sp.layout=list(bound), col.regions=heat.colors(32), 
       par.settings = list(strip.background=list(col="white")))
dev.off()


