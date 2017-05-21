library(maptools)
library(classInt)
library(GISTools)

setwd("/Volumes/dongmeic/fire/revision/figures")
area <-read.csv("landcover_area_table.csv")
head(area)
class(area$X0); class(area$X1)
png("forest_area.png", width=12, height=8, units="in", res=300)
plot(area$X0,area$X1/10000, type="l", xlab="Year", ylab="Forest Area (10,000 ha)", lwd =8, col="darkgreen", font.lab=8, font.axis=6, cex.axis=1.5, cex.lab=1.5)
#text(2005, 48000, "Afforestation in China", cex=3)
#text(2009, 43000, "Increasing Forested Area", cex=1.8, col="darkgreen")
dev.off()

png("forest_area_bar.png", width=12, height=8, units="in", res=300)
barplot(area$X1/10000000, names.arg=area$X0, main=expression('Annual forested area in China (10'^7*'ha)'), cex.names=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=3)
dev.off()

png("forest_area_plot.png", width = 2298, height = 1000, units = "px", res=300)
par(mar=c(4,5,1,1))
plot(area$X0, area$X1/10000000, type="h", lwd=8, xlab ="Year", ylab="Area", ylim=c(40, 50),
     cex.lab=1.3, cex.axis=1.3)
text(2006, 49, expression('Annual forested area in China (10'^7*'ha)'), cex=1.5)
axis(1,at=seq(40,50,by=2),labels=seq(40,50,by=2))
dev.off()

# burned forest area
ffpath <- "/Volumes/dongmeic/fire/forest_statistics/forest_fires/"
ffdf <- read.csv(paste(ffpath,"2002.csv", sep=""))
ffdf <- as.data.frame(as.character(ffdf$Region))
colnames(ffdf) <- "Region"
for (year in 2002:2013){
  df <- read.csv(paste(ffpath,year,".csv", sep=""))
  df$Total.area.burned <- as.numeric(as.character(df$Total.area.burned))
  df <- df[1:32,7]
  ffdf <- cbind(ffdf,df)
  print(year)
}
colnames(ffdf)[2:13] <- c("2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013")
ffdf$Sum <- rowSums(ffdf[,2:13])
head(ffdf)
pro.map <- readShapePoly("/Volumes/dongmeic/fire/China_pro_wgs84_dis.shp")
colnames(ffdf)[1] <- "ENAME"
pro.ffdf <- merge(pro.map, ffdf, by="ENAME")
proj4string(pro.ffdf) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ch.aea <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=30 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
pro.ffdf.tr <- spTransform(pro.ffdf,ch.aea)
pro.ffdf.tr$BAsum <- pro.ffdf.tr$BAsum * 25

add.northarrow <- function(){
  GISTools::north.arrow(xb=-1800000, yb=-600000, len=80000)
}

add.scale <- function(){
  GISTools::map.scale(-1800000, -880000, 600000,"100km",3,2,sfcol='brown')
} 

quick.map <- function(spdf,var,legend.title,main.title,color,outcome,style) {
  x <- spdf@data[,var]
  if (style == "jenks"){
    plotvar <- x
    nclr <- 5
    plotclr <- brewer.pal(nclr,color)
    class <- classIntervals(plotvar, nclr, style="jenks", dataPrecision=1)
    colcode <- findColours(class, plotclr)
    png(paste(outname,".png", sep=""), width=9, height=8, units="in", res=300)
    par(xpd=FALSE,mfrow=c(1,1),mar=c(0.5,0.5,2.5,0.5))
    plot(spdf, col=colcode)
    title(main.title,cex.main=1.5)
    legend(-680000,2480000, legend=names(attr(colcode, "table")),
           fill=attr(colcode, "palette"), title=legend.title, bty="n")
    add.northarrow()
    add.scale()
    dev.off()
  } else if (style == "quantileCuts") {
    shades <- auto.shading(x, cutter = quantileCuts, n=5, cols=brewer.pal(5, color))
    png(paste(outname,".png", sep=""), width=9, height=8, units="in", res=300)
    par(xpd=FALSE,mfrow=c(1,1),mar=c(0.5,0.5,2.5,0.5))
    choropleth(spdf,x,shades)
    title(main.title,cex.main=1.5)
    choro.legend(-680000,2480000, shades, bty='n', title=legend.title)
    add.northarrow()
    add.scale()
    dev.off()
  }else{
    cat("The style is either \"quantileCuts\" (only for spatial polygon) or \"jenks\". Please try again.")
  }
}

quick.map2 <- function(spdf,var,legend.title,main.title,color,style) {
  x <- spdf@data[,var]
  if (style == "jenks"){
    plotvar <- x
    nclr <- 5
    plotclr <- brewer.pal(nclr,color)
    class <- classIntervals(plotvar, nclr, style="jenks", dataPrecision=1)
    colcode <- findColours(class, plotclr)
    plot(spdf, col=colcode)
    title(main.title,cex.main=1.5)
    legend(-730000,2900000, legend=names(attr(colcode, "table")),
           fill=attr(colcode, "palette"), cex=0.95, title=legend.title, bty="n")
    add.northarrow()
    add.scale()
  } else if (style == "quantileCuts") {
    if (length(na.omit(x)) != length(x)){
      x[is.na(x)] <- 0
    }
    shades <- auto.shading(x, cutter = quantileCuts, n=5, cols=brewer.pal(5, color))
    choropleth(spdf,x,shades)
    title(main.title,cex.main=1.5)
    choro.legend(-730000,2900000, shades, bty='n', cex=0.95, title=legend.title)
    add.northarrow()
    add.scale()
  }else{
    cat("The style is either \"quantileCuts\" (only for spatial polygon) or \"jenks\". Please try again.")
  }
}

pro.ffdf.tr$BAsum_ann <- pro.ffdf.tr$BAsum/12
pro.ffdf.tr$Sum_ann <- pro.ffdf.tr$Sum/12

png("burned_area_ann.png", width=12, height=6, units="in", res=300)
par(xpd=FALSE,mfrow=c(1,2),mar=c(0.5,0.5,2.5,0.5))
quick.map2(pro.ffdf.tr,"BAsum_ann",expression('Burned area (ha·yr'^-1*')'), "Mean annual burned area from MODIS", "Reds", "quantileCuts")
quick.map2(pro.ffdf.tr,"Sum_ann",expression('Burned area (ha·yr'^-1*')'), "Mean annual burned area from CFSY", "Reds", "quantileCuts")
dev.off()

png("burned_area_ann.png", width=12, height=6, units="in", res=300)
par(xpd=FALSE,mfrow=c(1,2),mar=c(0.5,0.5,2.5,0.5))
shades <- shading(c(250, 1000, 4000, 8000), cols = brewer.pal(5, "Reds"))
choropleth(pro.ffdf.tr, pro.ffdf.tr$BAsum_ann, shades)
title("Mean annual burned area from MODIS",cex.main=1.5)
choro.legend(-870000,2900000, shades, bty='n', cex=1, title=expression('Burned area (ha·yr'^-1*')'))
add.northarrow()
add.scale()
choropleth(pro.ffdf.tr, pro.ffdf.tr$Sum_ann, shades)
title("Mean annual burned area from CFSY",cex.main=1.5)
choro.legend(-870000,2900000, shades, bty='n', cex=1, title=expression('Burned area (ha·yr'^-1*')'))
add.northarrow()
add.scale()
dev.off()
