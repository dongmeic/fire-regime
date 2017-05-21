library(gdalUtils)
library(raster)
library(rgdal)
library(maptools)
library(rasterVis)
library(classInt)
library(animation)

ch.shpfile <- "/Volumes/home/fire/chinabord/boundarynew.shp"
ch.shp <- readShapePoly(ch.shpfile, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
bound <- list("sp.polygons", ch.shp, col="grey", lwd=0.8, first = FALSE)
outpath <-  "/Volumes/home/fire/revision/maps/"
mask <- raster("/Volumes/home/fire/mask_new.tif")

s1 <- stack()
s2 <- stack()
for (month in 1:12){
  s.m.1 <- stack()
  s.m.2 <- stack()
  for (year in 2003:2016){
    r.file2 <- paste(inpath2, "CMG_Global_", year, "-", month, ".tif", sep="")
    r.file1 <- paste(inpath1, "CMG_Global_", year, "-", month, ".tif", sep="")
    r.file0 <- paste(inpath0, "CMG_Global_", year, "-", month, ".tif", sep="")
    r2 <- raster(r.file2)
    #print(range(na.omit(getValues(r2))))
    r1 <- raster(r.file1)
    r0 <- raster(r.file0)
    rc2 <- crop(r2, extent(mask))
    rc1 <- crop(r1, extent(mask))
    rc0 <- crop(r0, extent(mask))
    rc <- rc1 - rc0
    rc.m1 <- mask(rc, mask)
    rc.m2 <- mask(rc2, mask)
    rc.m1[rc.m1==0] <- NA
    rc.m2[rc.m2==0] <- NA
    s.m.1 <- stack(s.m.1, rc.m1)
    s.m.2 <- stack(s.m.2, rc.m2)
    print(paste(r.file1, year, "done!"))
  }
  print(paste("month", month, "processed!"))
  s1 <- stack(s1, s.m.1)
  s2 <- stack(s2, s.m.2)
}

ss1 <- stack()
k1 <- numeric()
ss2 <- stack()
k2 <- numeric()
for (i in 1:12){
  r1 <- calc(subset(s1, (i*14-13):(i*14)), sum)
  c1 <- na.omit(getValues(r1))
  print(range(c1))
  k1 <- c(k1, c1)
  r2 <- calc(subset(s2, (i*14-13):(i*14)), function(x) mean(x)/255)
  c2 <- na.omit(getValues(r2))
  k2 <- c(k2, c2)
  print(range(c2))
  ss1 <- stack(ss1, r1)
  ss2 <- stack(ss2, r2)
  print(paste("month", i, "processed!"))
} 
names(ss1) <- c("January", "Feburary", "March", "April", "May", "June", "July", "August", 
               "September", "October", "November", "December")
names(ss2) <- c("January", "Feburary", "March", "April", "May", "June", "July", "August", 
                "September", "October", "November", "December")

png(paste(outpath,"cmg_seasonality.png",sep=""), width=10, height=6, units="in", res=300)
spplot(ss1, sp.layout=list(bound), col.regions=rev(heat.colors(32)), at=c(0, 100, 300, 500, max(k1)), cuts=4,
       par.settings = list(strip.background=list(col="white")))
dev.off()

png(paste(outpath,"cmg_MeanCloundFraction.png",sep=""), width=10, height=6, units="in", res=300)
spplot(ss2, sp.layout=list(bound), col.regions=rev(heat.colors(32)), 
       par.settings = list(strip.background=list(col="white")))
dev.off()

# annual by month
for (year in 2008:2016){
  s.m <- stack()
  for (month in 1:12){
    r.file <- paste(inpath2, "CMG_Global_", year, "-", month, ".tif", sep="")
    r <- raster(r.file)
    rc <- crop(r, extent(mask))
    rc.m <- mask(rc, mask)
    rc.m <- rc.m/255
    s.m <- stack(s.m, rc.m)
    print(paste("year", year, "month", month, "processed!"))
  }
  names(s.m) <- c("January", "Feburary", "March", "April", "May", "June", "July", "August", 
                  "September", "October", "November", "December")
  png(paste(outpath,"cmg_MeanCloundFraction_", year, ".png",sep=""), width=10, height=6, units="in", res=300)
  print(spplot(s.m, sp.layout=list(bound), col.regions=rev(heat.colors(32)), main=paste("Mean Cloud Fraction in", year),
         par.settings = list(strip.background=list(col="white"))))
  dev.off()
  print(paste(year, "done!"))
}

im.convert(paste(outpath, "cmg_MeanCloundFraction_20*.png", sep = ""), output = paste(outpath, "cmg_MeanCloundFraction.gif", sep=""))

# annual
for (year in 2003:2016){
  s.m <- stack()
  for (month in 1:12){
    r.file <- paste(inpath2, "CMG_Global_", year, "-", month, ".tif", sep="")
    r <- raster(r.file)
    rc <- crop(r, extent(mask))
    rc.m <- mask(rc, mask)
    rc.m <- rc.m/255
    s.m <- stack(s.m, rc.m)
    print(paste("year", year, "month", month, "processed!"))
  }
  s <- calc(s.m, mean)
  names(s.m) <- c("January", "Feburary", "March", "April", "May", "June", "July", "August", 
                  "September", "October", "November", "December")
  png(paste(outpath,"MeanCloundFraction_", year, ".png",sep=""), width=10, height=6, units="in", res=300)
  print(spplot(s, sp.layout=list(bound), col.regions=rev(heat.colors(32)), main=paste("Annual Mean Cloud Fraction in", year)))
  dev.off()
  print(paste(year, "done!"))
}

# monthly active fire pixels
af.s.m <- stack()
path <- "/Volumes/home/fire/output/revision/results/"
r.file <- paste(path, "mean_monthly_active_fire.tif", sep="")
k <- numeric()
for (i in 1:12){
  r <- raster(r.file, band=i)
  c <- na.omit(getValues(r))
  print(range(c))
  k <- c(k, c)
  r.m <- mask(r, mask)
  r.m[r.m==0] <- NA
  af.s.m <- stack(af.s.m, r.m)
  print(paste("month", i, "processed!"))
}
names(af.s.m) <- c("January", "Feburary", "March", "April", "May", "June", "July", "August", 
                "September", "October", "November", "December")
png(paste(outpath,"active_fire_pixels.png",sep=""), width=10, height=6, units="in", res=300)
print(spplot(af.s.m, sp.layout=list(bound), col.regions=rev(heat.colors(32)), at=c(0, 25, 50, 100, max(k)), cuts=4,
             par.settings = list(strip.background=list(col="white"))))
dev.off()

# annual number of fires
r.file <- paste(path, "annual_number_of_fires.tif", sep="")
k <- numeric()
nof <- stack()
years <- 2001:2016
for (i in 1:16){
  r <- raster(r.file, band=i)
  c <- na.omit(getValues(r))
  print(sum(c))
  k <- c(k, sum(c))
  r.m <- mask(r, mask)
  nof <- stack(nof, r.m)
  print(paste("year", years[i], "processed!"))
}
df <- as.data.frame(t(rbind(years, k)))
plot(df$years, df$k)
