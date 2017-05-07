library(raster)
library(maptools)
library(parallel)
library(doParallel)
library(foreach)
registerDoParallel(cores=12)

path <- "/home2/dongmeic/fire/data/mcd14ml/"
out <- "/home2/dongmeic/fire/output/AF_China/"
mask <- raster("/home2/dongmeic/fire/mask.tif")

# settings
crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# read MCD14ML collection 6
foreach(year=2001:2015) %dopar% {
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
    afd.frp <- rasterize(spdf, mask, "FRP", fun=mean, na.rm=TRUE, crs=crs)
    out.frp <- paste(out, "afd", toString(year), "_", toString(as.numeric(month)), "_frp.tif", sep="")
    writeRaster(afd.frp, filename=out.frp, format="GTiff", overwrite=TRUE)
    print(paste(toString(year), month, "done!"))
  }
}

