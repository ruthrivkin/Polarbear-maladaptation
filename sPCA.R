library(lme4)
library(car)
library(dplyr)
library(tidyverse)
library(getopt)
library(adegenet)
library(ggplot2)
library(dartR)
library(hierfstat)
library(LEA)
library(poppr)
library(maps)
library(sp)
library(rworldmap)
library(mmod)
library(rnaturalearth)
library(plyr)
library(sf)

setwd("/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/Gradient_Forest/Subset/")


myCol <- c("brown","purple","green","orange","red","blue", "pink", "yellow", "darkgreen", "darkblue", "black", "darkgray", "lightgrey")
ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "NA"),
            axis.line = element_line(size=0.5), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=10, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,size=12),
            axis.title.x=element_text(vjust=0.1,size=12),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10))

#import files
infile <- "Subset.ld.hwe.raw"
myData <- read.PLINK(infile)
myData.nona <- gl.impute(myData, method = "HW")
myData.nona@pop

pop.info <- read.delim("sPCA_Subset_IDs.txt", header = F)
names(pop.info) <- c("subpop", "IID")

coord <- read.delim("sPCA_Subset_Coord.txt", header = FALSE)
plot(coord, pch = 19, cex = .5, 
     xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F)

gen <- gl2gi(myData.nona, probar = FALSE, verbose = NULL)


#Spatial pca (https://phylogatr.org/assets/modules/Spatial_Genetic_Structure.html)

#Add coordinates to genind, jitter coordinates to remove duplicate localities

coord.matrix <- as.matrix(coord)
coord.matrix=jitter(coord.matrix, factor=2)
#add jittered coordinates to genind object 
gen$other$xy<-coord.matrix

#Visualize sample density
nInd(gen)

#generate locality IDs and add assignments to genind object 
pop.info.location <- cbind(pop.info, coord)
coord1 <- pop.info.location %>% 
  dplyr::rename("Lon" = "V1",
                "Lat" = "V2")
pop.coords <- coord1  %>% 
  group_by(subpop)  %>% 
  dplyr::summarize(N = n(),
            Mean.Lon = mean(Lon),
            Mean.Lat = mean(Lat)
  )
pop.coords

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-141,-60, 50, 90) 
boundary

# Get map outlines from rworldmap package
library(rworldmap)
library(rworldxtra)
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
library(raster)
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
localsmap = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="white",
               colour="black", size=0.5)+
  geom_point( data=pop.info.location, aes(x=V1, y=V2, color = subpop), cex = 1.5) +
  coord_quickmap(expand=F)+
  ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  xlab("Longitude (°W)")+
  ylab("Latitude (°N)")+
  ng
localsmap


##perform a sPCA choosing connection network

library(spdep)
library(adespatial)


# Test spatial structure
pb.X <- scaleGen(gen)
mtest <- mantel.randtest(dist(pb.X), dist(gen$other$xy))
plot(mtest, nclass=30) #Not a clear signal

#check euclidean distance 
dist <- dist(gen$other$xy)
d <- as.matrix(dist)
range(d)

#create connection network using delaunay triangulation
myCn1 <- chooseCN(gen$other$xy, type=1, plot=TRUE)
myCn1
class(myCn1)


#run spca
#performs a sPCA using the my connection network; first run with all axes retained

mySpca.type1 <- spca(gen, cn = myCn1, ask=FALSE, scannf=FALSE) #so slow

barplot(mySpca.type1$eig, main="sPCA eigenvalues", col=spectral(length(mySpca.type1$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

#Check likelihood of global vs local structure, set permutations higher
myGtest <- global.rtest(gen$tab,mySpca.type1$lw,nperm=999)
myGtest
plot(myGtest)

myLtest <- local.rtest(gen$tab,mySpca.type1$lw,nperm=999)
myLtest #ns
plot(myLtest)

#produce a figure that represents eigenvalues of sPCA
screeplot(mySpca.type1) #Looks like the first 5 positive eigenvectors are significant

#reduced models
mySpca.reduced <- spca(gen,cn = myCn1, scannf=FALSE, nfposi=5, nfnega=0)
mySpca.reduced$as
plot.reduced <- plot(mySpca.reduced)


#plot eigenvalues of the analysis (stored inside the $eig component as a numeric vector) in decreasing order
#The presence of either global or local structure is indicated by a more extreme eigenvalue. 
#As each bar on the plot represents a single PC axis, we want to retain the axes with the most extreme (i.e., noticeably different) values. Global structure most important

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
cols = my.colors(5)



#visualize spatial genetic structure of sPCA axes, similar colors = genetically similar; distinct colors = genetically distinct
pdf("sPCA Colorplot.pdf", width = 6.78, height = 5.3)
colorplot(gen$other$xy, mySpca.reduced$ls, axes=1:3, col = my.colors, transp=T, alpha=0.75, cex=1.5, add=F, main="Colorplot of PC (lagged scores)", xlab = "Longitude", ylab = "Latitude")
maps::map(add = T, interior = F)
dev.off()

#Map the sPCA results using s.value and lagged scores ($ls), first and second axes
plot(coord, pch = 19, cex = 0, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)
s.value(gen$other$xy, mySpca.reduced$ls[,1], add.p=TRUE, csize=0.7)

#Check other axes
plot(coord, pch = 19, cex = 0, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)
s.value(gen$other$xy, mySpca.reduced$ls[,2], add.p=TRUE, csize=0.7)
s.value(gen$other$xy, mySpca.reduced$ls[,3], add.p=TRUE, csize=0.7)

