setwd("/Volumes/dongmeic/fire/revision/figures")
area <-read.csv("landcover_area_table.csv")
head(area)
class(area$X0); class(area$X1)
png("forest_area.png", width=12, height=8, units="in", res=300)
plot(area$X0,area$X1/10000, type="l", xlab="Year", ylab="Forest Area (10,000 ha)", lwd =16, col="darkgreen", font.lab=4, font.axis=4)
text(2005, 22000, "Afforestation in China", cex=4)
text(2009, 18500, "Increasing Forest Area", cex=3, col="darkgreen")
dev.off()

png("forest_area_bar.png", width=12, height=8, units="in", res=300)
barplot(area$X1/10000, names.arg=area$X0, main = "Annual forested area in China (10,000 ha)", cex.names=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=3)
dev.off()

