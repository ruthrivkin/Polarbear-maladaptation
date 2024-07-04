# Gradient Forest analyses on whale data!
setwd("/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/Gradient_Forest/GFAllEnvironments//")
# Lots of help for gradient forest analyses from pgugger's github (https://github.com/pgugger/LandscapeGenomics/blob/master/2019/Exercise4.md)

# For converting vcf to gradient forest file I used beginning of tutorial here (came across this one before the one above, but from same person: https://github.com/pgugger/LandscapeGenomics/blob/master/2018_China/Exercise3.md)
#prep files for gradient forest

#Create compatible file with rows with ID and columns with snps, no na values
library(adegenet)
library(dplyr)
library(tidyverse)
library(LEA)
library(dartR)

#load pop info, remove AB pop and subpop id
pop.info <- read.delim("Subset_IDs.txt", header = F)
names(pop.info) <- c("subpop", "IID")


#Can also do with snmf
vcfin <- "../Subset/forgf.vcf"
vcf2lfmm(vcfin)

genoin <- "../Subset/forgf.geno"

project = snmf(genoin, K = 1:13, ploidy = 2, entropy = T,
                alpha = 100, project = "new", repetitions = 10, seed = 42)
summary(project)

project=load.snmfProject("../Subset/forgf.snmfProject")
ce = cross.entropy(project, K = 5)
lowest.ce = which.min(ce)
impute(project, "../Subset/forgf.geno", method = 'mode', K = 5, run = lowest.ce)

lfmm_data <- "../Subset/forgf.lfmm_imputed.lfmm"
lfmm2geno(lfmm_data)

#LD pruned dataset
imputedgeno <- as.data.frame(read.geno("../Subset/forgf.lfmm_imputed.geno"))

#Read in snp names
snp_info <- read.table("../Subset/forgf.map")
head(snp_info)
# pull out only CHROM and POS
snp_info <- dplyr::select(snp_info, V2)
colnames(snp_info) <- "snp"
#transform snps to columns and set as column names for genotypes
snp_info_t <- as.data.frame(t(snp_info))
colnames(imputedgeno) <-snp_info_t

snps <- cbind(pop.info, imputedgeno) 
gf.snp <- snps %>%
  column_to_rownames(var="IID")   %>% 
  dplyr::select(-"subpop")
#Save
write.table(gf.snp, "snps.ld.gradientforest.forR")




# And help from gradient forest vignette: (https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)



# need to make sure rtools is downloaded (https://cran.rstudio.com/bin/windows/Rtools/rtools42/rtools.html)\
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

# load libraries
library(gradientForest)
library(dplyr)
library(vegan)

# for climate data stuff
library(raster)
library(tidyverse)
library(sdmpredictors)
library(sp)
library(rgdal)

# some more to help with editing map plots
library(rnaturalearth)
library(rworldmap)
library(ggspatial)


# Step 1: Initial run with gradient forest on present data
# Step 2: Create plots for biological space and geographical space
# Step 3: Measuring genetic offset with future data

##########  
##### Step 1: Initial run with gradient forest on present data

## 1.1 prepare data files

# Load snp data - note that data file must not have missing data for GF
snp <- read.table("snps.ld.gradientforest.forR", header=T, row.names=1)
dim(snp)

# use fread() in the data.table package for bowhead with large data set
#library(data.table)
#snp_temp <- fread("bowhead_snps.filter1.miss.biallel.min100kb.autosomes.maf.n21.miss01.imputed.thin1000.forR", header=T)
#snp_temp <- as.data.frame(snp_temp)
#snp <- snp_temp[,-1]
#rownames(snp) <- snp_temp[,1]


# Load environmental data for each site 

site.clim.dat <- read.csv("23.07.20_tempandice_withcoord.csv")

str(site.clim.dat)

# Pull out present data (I also had future data in site.clim.dat)
pres.clim.points <- dplyr::select(site.clim.dat, "IID", "subpop", "Longitude", "Latitude", "BO22_icethickmean_ss", "BO22_icecovermean_ss", "MeanTemp")
colnames(pres.clim.points)
pres.clim.points.thick <- dplyr::select(site.clim.dat, "IID", "subpop", "Longitude", "Latitude", "BO22_icethickmean_ss", "MeanTemp")
pres.clim.points.cov <- dplyr::select(site.clim.dat, "IID", "subpop", "Longitude", "Latitude", "BO22_icecovermean_ss", "MeanTemp")

# can rename here if needed
pres.clim.points <- rename(pres.clim.points, "IceThickness" = "BO22_icethickmean_ss", "IceCover" = "BO22_icecovermean_ss", "Temp" = "MeanTemp")
pres.clim.points.thick <- rename(pres.clim.points.thick, "IceThickness" = "BO22_icethickmean_ss", "Temp" = "MeanTemp")
pres.clim.points.cov <- rename(pres.clim.points.cov, "IceCover" = "BO22_icecovermean_ss", "Temp" = "MeanTemp")

cor.test(pres.clim.points$IceThickness, pres.clim.points$IceCover)
cor.test(pres.clim.points$IceThickness, pres.clim.points$Temp)


## 1.2 Generate PCNM spatial variables 

# Thought about distances between each individual/site. I used distances within water for isolation-by-distance stuff, but here using general lat/long that is not incorporating water/land, just straight distance. Generally OK though for this purpose, so going forward with just using the site coordinates.
library(adespatial)
library(spatialRF)
coord <- pres.clim.points[,c('Longitude','Latitude')]

#create distance matrix
dist <- as.matrix(dist(coord))
class(dist)
mem <- mem(dist)
pcnm <- pcnm(dist(coord)) 

# Keep half of positive ones
keep <- round(length(which(pcnm$value > 0))/2) 
pcnm.keep <- scores(pcnm)[,1:keep]
pcnm.keep

keep1 <- round(length(which(mem$value > 0))/2) 
mem.keep <- scores(mem)[,1:keep1]
mem.keep

# Create file with climate and PCNM variables (no lat/long)
env.gf.pcnm <- cbind(pres.clim.points, pcnm.keep)
str(env.gf.pcnm)
env.gf.pcnm <- env.gf.pcnm[,-c(1:4)]

env.gf.thick<- cbind(pres.clim.points.thick, pcnm.keep)
str(env.gf.thick)
env.gf.thick <- env.gf.thick[,-c(1:4)]

env.gf.cov<- cbind(pres.clim.points.cov, pcnm.keep)
str(env.gf.cov)
env.gf.cov <- env.gf.cov[,-c(1:4)]

# Define max number of splits (what does this mean?)
lev.pcnm <- floor(log2(0.368*nrow(env.gf.pcnm)/2))
lev.pcnm.thick <- floor(log2(0.368*nrow(env.gf.thick)/2))
lev.pcnm.cov <- floor(log2(0.368*nrow(env.gf.cov)/2))

# run GF using the climate, spatial, and snp data inputs
gf <- gradientForest(cbind(env.gf.pcnm, snp), 
                     predictor.vars=colnames(env.gf.pcnm),
                     response.vars=colnames(snp), 
                     ntree=500, 
                     maxLevel=lev.pcnm, 
                     trace=T, 
                     corr.threshold=0.50,
                     nbin=201) 

gf.thick<- gradientForest(cbind(env.gf.thick, snp), 
                     predictor.vars=colnames(env.gf.thick),
                     response.vars=colnames(snp), 
                     ntree=500, 
                     maxLevel=lev.pcnm.thick, 
                     trace=T, 
                     corr.threshold=0.50,
                     nbin=201) 

gf.cov<- gradientForest(cbind(env.gf.cov, snp), 
                          predictor.vars=colnames(env.gf.cov),
                          response.vars=colnames(snp), 
                          ntree=500, 
                          maxLevel=lev.pcnm.cov, 
                          trace=T, 
                          corr.threshold=0.50,
                          nbin=201) 

#301 for narwhal,
# idk if increasing or decreasing nbin is to help with speed or dataset size. 
#experimenting:
#200 bins for narwhal took 1hr 22 min
#301 bins for narwhal took 45 min
# ok so, this means, more bins -> less time (to an extent probably).
# For bowhead, ended up using ntree=100, and nbin=601.

# save here too.
saveRDS(gf.cov, "gf.ld_pcnm_covk.RDS")
gf <- readRDS("gf.ld_pcnm_thick.RDS")

 ## 1.3 Initial gradient forest plots
# tip: use help(par) to see all the plotting options for color, size, etc

# bar graphs of predictor overall variable importance
# chose colors to highlight PCNMs and known environmental variables separately (orange for PCNM, blue for variables)
pdf("GF_predictor_importance.pdf")
plot(gf.thick, plot.type = "O", 
     lwd=1, cex.axis=0.9)
dev.off()

# save as svg

#Calculate total R2 explained for everything and for spatial variables and env variables
importance(gf, type = "Weighted")

imp <- as.data.frame(importance(gf.cov, type = "Weighted"))
imp
summary(imp)

# Organize variables by importance, for other plots, select only the top 4
by.importance <- names(importance(gf.cov))[1:5]
by.importance
par(mgp = c(3,1,0))

# Splits density plot (showing binned split importance and location on each gradient)
pdf("GF_split_density.pdf", width = 6.78, height = 5.3)
plot(gf.thick, plot.type="S", imp.vars=by.importance,
     leg.posn="topright", cex.legend=0.7, cex.axis=0.8,cex.lab = 0.7, line.ylab=0.9, 
     par.args=list(mgp=c(1.5,0.5,0), mar = c(2.5, 1, 0.1, 0.5), omi=c(0.2,0.35,0.2,0.4)))
dev.off()
# Plot turnover functions, showing how allelic composition changes along spatial or environmental gradients. this is a predictor cumulative plot
pdf("GF_turnover_thickness.pdf", width = 6.78, height = 5.3)
plot(gf.cov, plot.type="C", imp.vars=c("IceCover", "Temp"), 
     show.species=F, common.scale=T, cex.axis=1, cex.lab=1.2, line.ylab=1, lwd=2,col="#7393B3",
     par.args=list(mgp=c(1.5,0.5,0),  mar = c(2.5, 1, 0.1, 0.5), omi=c(0.2,0.35,0.2,0.4)))
dev.off()
# can't get it to plot 4x2 instead of 2x4 but can adjust in inkscape

# Plot turnover functions for individual loci. each line represent allelic change at a single SNP
# (can take a long time to load. maybe skip for now)
pdf("GF_turnover_loci.pdf", width = 6.78, height = 5.3)
plot(gf, plot.type="C", imp.vars=by.importance, 
     show.overall=F, legend=T, leg.posn="topleft", leg.nspecies=5, 
     cex.lab=1, cex.legend=0.5, cex.axis=1, line.ylab=0.9, 
     par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0.2,0.35,0.2,0.4)))
dev.off()
 
 # Plot R2 measure of the fit of the random forest model 
plot(gf, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)

##########  
##########  
##### Step 2: Create plots for biological space and geographical space
# Names of present variables of interest

## 2.1 Load climate data rasters
# Note: instead of re-loading the layers from online with sdmpredictors every time I run this during load_layers(), I  downloaded and saved the rasters into a local directory so can load faster from there

# Names of present variables of interest 
#BO22_RCP26_2050_icethickmean_ss","BO22_RCP45_2050_icethickmean_ss","BO22_RCP60_2050_icethickmean_ss","BO22_RCP85_2050_icethickmean_ss", BO22_RCP26_2100_icethickmean_ss","BO22_RCP45_2100_icethickmean_ss","BO22_RCP60_2100_icethickmean_ss","BO22_RCP85_2100_icethickmean_ss"
ice.present.vars <- c("BO22_icethickmean_ss")

ice.2050.raster<- "BO22_RCP85_2050_icethickmean_ss"
ice.2100.raster<- "BO22_RCP85_2100_icethickmean_ss"

ice.present.rasters <- load_layers(ice.present.vars, datadir="./Rasters/")
ice.2050.raster <- load_layers(ice.2050.raster, datadir="./Rasters/")
ice.2100.raster <- load_layers(ice.2100.raster, datadir="./Rasters/")


#load temp data now that everything is at the same extent, and merge with ice rasters
temp.rast <- raster("./Rasters/MeanTemp_Current_StableClim.tif")
temp.2050.rast <- raster("./Rasters/MeanTemp_2050_StableClim.tif")
temp.2100.rast <- raster("./Rasters/MeanTemp_2100_StableClim.tif")

plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice t")
plot(ice.present.rasters, col=my.colors(200), add = TRUE)

#temp.hist.rast <- raster("./Rasters/MeanTemp_historical_1960_1988.tif")

#convert kelvin to celsius 
#temphist.rast.celsius <- temp.hist.rast - 273.15

presenttemp = projectRaster(temp.rast, ice.present.rasters)
present = c(ice.present.rasters, presenttemp)
present.rasters <- raster::brick(present)

temp.2050 = projectRaster(temp.2050.rast, ice.2050.raster)
env.2050 = c(ice.2050.raster, temp.2050)
env.2050.rasters <- raster::brick(env.2050)

temp.2100 = projectRaster(temp.2100.rast, ice.2100.raster)
env.2100 = c(ice.2100.raster, temp.2100)
env.2100.rasters <- raster::brick(env.2100)

#temp.hist = projectRaster(temphist.rast.celsius, ice.2100.raster)
#env.hist = c(ice.2100.raster, temp.hist)
#env.histrasters <- raster::brick(env.hist)

names(present.rasters) <- c("IceThickness", "Temp")
names(env.2050.rasters) <- c("IceThickness","Temp")
names(env.2100.rasters) <- c("IceThickness", "Temp")


## 2.2 Crop to desired area, using species' range here
# Load range shapefile (obtained from IUCN red list database)
species <- readOGR('../redlist_species_data_b48ae38b-6561-49b9-846f-033115599d7e/data_0.shp')

# Initial crop- probably not necessary but it works for this pipeline. can make more efficient later

crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent
crop.map.pres <- crop(present.rasters,crop.extent) #crop rasters
crop.2050 <- crop(env.2050.rasters, crop.extent)
crop.2100 <- crop(env.2100.rasters, crop.extent)



# Crop to species range shape file
present.rasters.range <- mask(crop.map.pres, species)
crop.2050.rasters.range <- mask(crop.2050, species)
crop.2100.rasters.range <- mask(crop.2100, species)

#r.points <- rasterToPoints(data.rasters.crop, spatial=FALSE)

# column names
#colnames(r.points) <- c("Lon", "Lat", "sst","icethick", "salinity", "chloro")

## 2.3 Extract climate data from rasters and transform environmental predictors
# start with present data
clim.land <- raster::extract(present.rasters.range, 1:ncell(present.rasters.range), df = TRUE)
clim.land <- na.omit(clim.land)
row.names(clim.land) <- NULL
clim.land.noid <-(clim.land[,c("IceThickness", "Temp")])
clim.land.noid <-as.data.frame(clim.land.noid)
class(clim.land)
class(clim.land.noid)

# use predict function to transform environmental predictors
pred <- predict(gf, clim.land.noid)  #note the removal of the cell ID column with [,-1])

#ggsave("blank_polar_map_lat45_linesonly.svg", width=5, height=5)

##########  

## 2.3 Extract climate data from rasters and transform environmental predictors
# start with present data
# convert predictions to color scale for map
# use pca on predictions
pca <- prcomp(pred, center=T, scale.=F)
summary(pca)

# assign pcs to colors
r <- pca$x[, 1]
g <- pca$x[, 2]
b <- pca$x[, 2]

#scale colors
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255

#define raster properties with an existing one
mask<-present.rasters.range$IceThickness
mask[]<-as.numeric(mask[]>0)

# assign color to raster
rastR <- rastG <- rastB <- mask

rastR[clim.land$ID] <- r
rastG[clim.land$ID] <- g
rastB[clim.land$ID] <- b

#stack color raster
rgb.rast <- raster::brick(rastR, rastG, rastB)

#pdf("BOW_scaffold0_GF_Map.pdf")
# initial map/plot
pdf("GF_PCA.map.pdf")
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85), col = "lightgrey",
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plotRGB(rgb.rast, bgalpha = 1, add = TRUE, col = "lightgrey")
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()
# pca thing


nvs <- dim(pca$rotation)[1]
vec <- c("IceThickness", "Temp")
lv <- length(vec)
vind <- rownames(pca$rotation) %in% vec
scal <- 1000
xrng <- range(pca$x[, 1], pca$rotation[, 1]/scal) * 1.1
yrng <- range(pca$x[, 2], pca$rotation[, 2]/scal) * 1.1

pdf("GF_PCA_loadings_thickness.pdf")
plot((pca$x[, 1:2]), xlim = xrng, ylim = yrng, 
     pch = ".", cex = 4, col = rgb(r, g, b, max = 255),
     asp = 1)

arrows(rep(0, lv), rep(0, lv), pca$rotation[vec,1]/2000, pca$rotation[vec, 2]/2000, length = 0.05)
jit <- 0.0015

text(pca$rotation[vec, 1]/100 + jit * sign(pca$rotation[vec, 1]), pca$rotation[vec, 2]/1000 + jit * sign(pca$rotation[vec, 2]), labels = vec)
dev.off()

"The colors represent genetic variation (allelic composition) as predicted based on the modeled relationships with environmental and spatial variables. Similar colors are predicted to be more similar genetically."


##########  
##### Step 3: Initial run with gradient forest on present data

## 3.1 prepare data filesm run through each scenario and merge at the end

#RCP8.5 2050
#Run each time for each section
clim.2050 <- raster::extract(crop.2050.rasters.range, 1:ncell(crop.2050.rasters.range), df = TRUE)
clim.2050 <- na.omit(clim.2050)
row.names(clim.2050) <- NULL
clim.2050.noid <-(clim.2050[,c("IceThickness","Temp")])
clim.2050.noid <-as.data.frame(dplyr::select(clim.2050, c("IceThickness","Temp")))

# transform environmental variables
pred.2050 <- predict(gf, clim.2050.noid)

# noting here that changed pred.adaptive (tutorial) to pred (here).
# genetic offset:
genetic.offset.adaptive <- sqrt((pred.2050[,1]-pred[,1])^2) + 
                                  (pred.2050[,2]-pred[,2])^2

# Define raster properties --- the variable doesn't matter here i think
offset2050 <- crop.2050.rasters.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offset2050[clim.2050$ID] <- genetic.offset.adaptive

#Plot with color scale of choice

plot(offset2050, col=rev( rainbow( 99, start=0, end=0.2)))

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

pdf("GeneticOffset_2050.pdf",  width = 6.78, height = 5.3)
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plot(offset2050, col=my.colors(200), add = TRUE)
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()

#RCP8.5 2100
#Run each time for each section
clim.2100 <- raster::extract(crop.2100.rasters.range, 1:ncell(crop.2100.rasters.range), df = TRUE)
clim.2100 <- na.omit(clim.2100)
row.names(clim.2100) <- NULL
clim.2100.noid <-(clim.2100[,c("IceThickness","Temp")])
clim.2100.noid <-as.data.frame(dplyr::select(clim.2100, c("IceThickness","Temp")))


# transform environmental variables
pred.2100 <- predict(gf, clim.2100.noid)

# noting here that changed pred.adaptive (tutorial) to pred (here).
# genetic offset:
genetic.offset.adaptive.2100 <- sqrt((pred.2100[,1]-pred[,1])^2) + 
  (pred.2100[,2]-pred[,2])^2

# Define raster properties --- the variable doesn't matter here i think
offset2100 <- crop.2100.rasters.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offset2100[clim.2100$ID] <- genetic.offset.adaptive.2100

#Plot with color scale of choice

pdf("GeneticOffset_2100.pdf",  width = 6.78, height = 5.3)
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2100 - 8.5 degree increase")
plot(offset2100, col=my.colors(200), add = TRUE)
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()

#Historic backtrack (1960-1988), keeping so I have the info later
#Run each time for each section
clim.hist <- raster::extract(crop.hist.rasters.range, 1:ncell(crop.hist.rasters.range), df = TRUE)
clim.hist <- na.omit(clim.hist)
row.names(clim.hist) <- NULL
clim.hist.noid <-(clim.hist[,c("IceThickness","Temp")]) #Ice thickness is from 2100 but I can't figure out how to do this with only Temp because the extent is different
clim.hist.noid <-as.data.frame(dplyr::select(clim.hist, "Temp"))

# transform environmental variables
pred.hist <- predict(gf, clim.hist.noid)

# noting here that changed pred.adaptive (tutorial) to pred (here).
# genetic offset:
genetic.offset.adaptive.hist <- sqrt((pred.hist[,2]-pred[,2])^2)

# Define raster properties --- the variable doesn't matter here i think
offsethist <- crop.hist.rasters.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offsethist[clim.hist$ID] <- genetic.offset.adaptive.hist



######## extract gradient forest values?
# save the rast.offset for later too
#writeRaster(offset2050, "raster.offset.2050.tif", overwrite = TRUE)
#writeRaster(offset2100, "raster.offset.2100.tif")

# can import raster back in
#rast.offset.26 <- raster("raster.offset.2050.tif")
#rast.offset.45 <- raster("raster.offset.2100.tif")

# treat rast.offset as a raster to extract values by site.
# load subpop locations
colnames(pres.clim.points)
x <- pres.clim.points$Longitude
y <- pres.clim.points$Latitude

pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
projection(pts)==projection(offset2050) # make sure projection matches

#Extract offset values and merge into dataframe
offset.data.2050 <- as.data.frame(raster::extract(offset2050, pts))
names(offset.data.2050) <- "offset.2050"
offset.all <- cbind(pres.clim.points, offset.data.2050)

offset.data.2100 <- as.data.frame(raster::extract(offset2100, pts))
names(offset.data.2100) <- "offset.2100"
offset.all <- cbind(offset.all, offset.data.2100)


#save dataset
write.csv(offset.all, "GF.ld.offset.all.csv")

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "black"),
            axis.line = element_line(size=1), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=15, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,face="bold",size=17),
            axis.title.x=element_text(vjust=0.1,face="bold",size=17),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15))
# legend.position = "top", legend.direction="horizontal", 
#legend.text=element_text(size=12), legend.key = element_rect(fill = "white"), 
#legend.title = element_text(size=12),legend.key.size = unit(1, "cm"))

#load dataset 
#summarize by heterozygosity and missing data

data <- read.csv("GF.offset.all.csv")
colnames(data)

o.2050 <- lm(offset.2050 ~ subpop,
             data = data)
summary(o.2050)
car::Anova(o.2050, type = 2)

o.2100 <- lm(offset.2100 ~ subpop,
             data = data)

summary(o.2100)
car::Anova(o.all, type = 2)


o.hist <- lm(offset.historical ~ subpop,
             data = data)

summary(o.hist)
car::Anova(o.hist, type = 2)

#Check if offset values change significantly with projections
colnames(data)

datalong <- dplyr::select(data, c("subpop", "offset.2050", "offset.2100", "offset.historical"))
datalong = reshape2::melt(datalong, id.vars=c("subpop"),
                          variable.name = "Year", 
                          value.name = "Offset")
str(datalong)

climate.scenarios <- lm(Offset ~ Year,
                        data = datalong)

summary(climate.scenarios)
car::Anova(climate.scenarios, type = 2)


#Plot some stuff

means <- Rmisc::summarySE(datalong, measurevar = "Offset", groupvar = c("subpop", "Year"), na.rm = TRUE)
str(means)

future.scenarios <- means %>%
  filter(!Year=="offset.historical")

hist.scenarios <- means %>%
  filter(Year=="offset.historical")

(offsetplot <- ggplot(means, aes(x=reorder(subpop, Offset, na.rm = T), y=Offset, fill=Year)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=.3) +      # Thinner lines
    geom_errorbar(aes(ymin=Offset-se, ymax=Offset+se),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    scale_y_continuous(limits = c(0, 0.002)) + 
    xlab("Subpopulation") +
    ylab("Genetic Offset")+
    scale_fill_discrete(name="Climate Scenario") +
    ng
)
ggsave("Offset_by_subpop.pdf", width = 6.78, height = 5.3, dpi = 300)

(histoffsetplot <- ggplot(hist.scenarios, aes(x=reorder(subpop, Offset, na.rm = T), y=Offset, fill = Year)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=.3) +      # Thinner lines
    geom_errorbar(aes(ymin=Offset-se, ymax=Offset+se),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Subpopulation") +
    ylab("Genetic Offset")+
    scale_fill_discrete(name="Climate Scenario") +
    ng
)

ggsave("Historical_Offset_by_subpop.pdf", width = 6.78, height = 5.3, dpi = 300)


str(env.gf.pcnm)
#Create map of samples color coded by offset value
library(raster)
library(sf)
# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-141,-60, 50, 90) 
boundary

# Get map outlines from rworldmap package
library(rworldmap)
library(rworldxtra)
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

library(gridExtra)
# Plot basemap

(PCNM2 = ggplot()+
    geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="lightgrey",
                 colour="black")+
    geom_point( data=env.gf.pcnm, aes(x=Longitude, y=Latitude, col = PCNM2), cex = 1) +
    scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
    scale_shape_manual(values = c(3, 1,2,5,6,4)) +
    coord_quickmap(expand=F)+
    ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
    xlab("Longitude (°W)")+
    ylab("Latitude (°N)")+
    
    ng
)
ggsave("Map_points_coded_by_PCNM2.pdf", width = 6.78, height = 5.3, dpi = 300)


# Plot basemap
(PCNM4 = ggplot()+
    geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="lightgrey",
                 colour="black")+
    geom_point( data=env.gf.pcnm, aes(x=Longitude, y=Latitude, col = PCNM4), cex = 1) +
    scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
    scale_shape_manual(values = c(3, 1,2,5,6,4)) +
    coord_quickmap(expand=F)+
    ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
    xlab("Longitude (°W)")+
    ylab("Latitude (°N)")+
    
    ng
)

ggsave("Map_points_coded_by_PCNM4.pdf", width = 6.78, height = 5.3, dpi = 300)

