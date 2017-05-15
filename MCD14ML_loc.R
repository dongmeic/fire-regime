library(raster)
library(maptools)

path <- "/Volumes/dongmeichen/mcd14ml/"
out <- "/Volumes/dongmeichen/output/AF_China/"
mask <- raster("/Volumes/dongmeichen/masknew.tif")
# chinabd <- readShapePoly("/Volumes/dongmeichen/boundarynew.shp")
# proj4string(chinabd) <- crs

# settings
crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# read MCD14ML collection 6
for (year in 2001:2015){
  df.yr <- data.frame(YYYYMMDD=integer(0), HHMM=integer(0), sat=character(0), lat=numeric(0), 
                        lon=numeric(0), T21=numeric(0), T31=numeric(0), sample=numeric(0),
                        FRP=numeric(0), conf=numeric(0), type=numeric(0), DOY=numeric(0))
  for (month in formatC(1:12, width=2, flag="0")){
    df <- read.table(gzfile(paste(path, "MCD14ML.", toString(year), month, ".006.01.txt.gz", sep="")), header=T)
    df$DOY <- as.POSIXlt(as.character(df$YYYYMMDD), format="%Y%m%d")$yday + 1
    head(df)
    # select confidence higher than zero and hot spot type is vegetation fire
    df.n <- df[df$conf != 0 & df$type == 0,]
    xy <- data.frame(df.n[,c("lon","lat")])
    coordinates(xy) <- c("lon","lat")
    proj4string(xy) <- crs
    spdf <- SpatialPointsDataFrame(coords = xy, data=df.n, proj4string = crs)
    afd.doy <- rasterize(spdf, mask, "DOY", fun='last', na.rm=TRUE, crs=crs)
    out.doy <- paste(out, "afd", toString(year), "_", toString(as.numeric(month)), ".tif", sep="")
    writeRaster(afd.doy, filename=out.doy, format="GTiff", overwrite=TRUE)
    df.yr <- rbind(df.yr, df.n)
    print(paste(toString(year), month, "done!"))
  }
  xy.yr <- data.frame(df.yr[,c("lon","lat")])
  coordinates(xy.yr) <- c("lon","lat")
  proj4string(xy.yr) <- crs
  spdf.yr <- SpatialPointsDataFrame(coords = xy.yr, data=df.yr, proj4string = crs)
  afd.yr.frp <- rasterize(spdf.yr, mask, "FRP", fun=mean, na.rm=TRUE, crs=crs)
  afd.yr.doy <- rasterize(spdf.yr, mask, "DOY", fun='last', na.rm=TRUE, crs=crs)
  out.yr.frp <- paste(out, "afd_", toString(year), "_frp.tif", sep="")
  out.yr.doy <- paste(out, "afd_", toString(year), "_doy.tif", sep="")
  writeRaster(afd.yr.frp, filename=out.yr.frp, format="GTiff", overwrite=TRUE)
  print(paste(out.yr.frp, "done!"))
  writeRaster(afd.yr.doy, filename=out.yr.doy, format="GTiff", overwrite=TRUE)
  print(paste(out.yr.doy, "done!"))
}