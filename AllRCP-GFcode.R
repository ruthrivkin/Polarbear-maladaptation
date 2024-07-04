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

#load pop info
pop.info <- read.delim("Subset_IDs.txt", header = F)
names(pop.info) <- c("subpop", "IID")


#Impute genotypes with snmf. We have to run a new snmf analysis because we are not filtering for hwe here... the same number of clusters was selected in both cases
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



#Run the GF analysis (https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)
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
snp <- read.table("Chip12_IceBears_snptable.USETHISONE.forR", header=T, row.names=1)

# Load environmental data for each site 

site.clim.dat <- read.csv("24.04.15_temp_withallcoord.csv")
str(site.clim.dat)

# Pull out present data (I also had future data in site.clim.dat)
pres.clim.points <- dplyr::select(site.clim.dat, "IID", "subpop", "Longitude", "Latitude", "BO22_icethickmean_ss", "BO22_icecovermean_ss", "MeanTemp")
colnames(pres.clim.points)
pres.clim.points.thick <- dplyr::select(site.clim.dat, "IID", "subpop", "Longitude", "Latitude", "BO22_icethickmean_ss", "MeanTemp")

# can rename here if needed
pres.clim.points.thick <- rename(pres.clim.points.thick, "IceThickness" = "BO22_icethickmean_ss", "Temp" = "MeanTemp")
pres.clim.points.cov <- rename(pres.clim.points.cov, "IceCover" = "BO22_icecovermean_ss", "Temp" = "MeanTemp")

cor.test(pres.clim.points$IceThickness, pres.clim.points$IceCover) #These are highly correlated so running them in seperate models
cor.test(pres.clim.points$IceThickness, pres.clim.points$Temp)
cor.test(pres.clim.points.thick$IceCover, pres.clim.points.thick$Temp)

## 1.2 Generate PCNM spatial variables 

library(adespatial)
library(spatialRF)
coord <- pres.clim.points[,c('Longitude','Latitude')]

#create distance matrix
dist <- as.matrix(dist(coord))
class(dist)
pcnm <- pcnm(dist(coord)) 

# Keep half of positive ones
keep <- round(length(which(pcnm$value > 0))/6) 
pcnm.keep <- scores(pcnm)[,1:10]
pcnm.keep

# Create file with climate and PCNM variables (no lat/long)

env.gf.thick<- cbind(pres.clim.points.thick, pcnm.keep)
str(env.gf.thick)
env.gf.thick <- env.gf.thick[,-c(1:4)]
e
nv.gf.cov<- cbind(pres.clim.points.cov, pcnm.keep)
str(env.gf.cov)
env.gf.cov <- env.gf.cov[,-c(1:4)]

# Define max number of splits
lev.pcnm.thick <- floor(log2(0.368*nrow(env.gf.thick)/2))
lev.pcnm.cov <- floor(log2(0.368*nrow(env.gf.cov)/2))

# run GF using the climate, spatial, and snp data inputs
head(env.gf.thick)
snp[1:5,1:5]

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


# save here too.
saveRDS(gf.thick, "gf_pcnm_thick.RDS")
saveRDS(gf.cov, "gf_pcnm_cov.RDS")

#Load thickness, can do the exact same analysis with ice cover just change the file name
gf <- readRDS("gf.ld_pcnm_thick.RDS")

summary(gf)
head(gf$Y)
gf$Y[1:5,1:5]

(gf$imp.rsq)

## 1.3 Initial gradient forest plots
# tip: use help(par) to see all the plotting options for color, size, etc

# bar graphs of predictor overall variable importance
pdf("GF_predictor_importance.pdf")
plot(gf, plot.type = "O", 
     lwd=1, cex.axis=0.9)
dev.off()

# save as svg

#Calculate total R2 explained for everything and for spatial variables and env variables
importance(gf, type = "Weighted")

imp <- as.data.frame(importance(gf, type = "Weighted"))
imp

# Organize variables by importance, for other plots, select only the top 4
by.importance <- names(importance(gf))[1:5]
by.importance
par(mgp = c(3,1,0))

# Splits density plot (showing binned split importance and location on each gradient)
pdf("GF_split_density.pdf", width = 6.78, height = 5.3)
plot(gf, plot.type="S", imp.vars=c("IceThickness", "Temp"),
     leg.posn="topright", cex.legend=0.7, cex.axis=0.8,cex.lab = 0.7, line.ylab=0.9, 
     par.args=list(mgp=c(1.5,0.5,0), mar = c(2.5, 1, 0.1, 0.5), omi=c(0.2,0.35,0.2,0.4)))
dev.off()

# Plot turnover functions, showing how allelic composition changes along spatial or environmental gradients. this is a predictor cumulative plot
pdf("GF_turnover_thick.pdf", width = 6.78, height = 5.3)
plot(gf, plot.type="C", imp.vars=c("IceThickness", "Temp"), 
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

## 2.1 Load climate data rasters.. this step is such a pain
# Note: instead of re-loading the layers from online with sdmpredictors every time I run this during load_layers(), I  downloaded and saved the rasters into a local directory so can load faster from there

# Names of present variables of interest 
#BO22_RCP26_2050_icethickmean_ss","BO22_RCP45_2050_icethickmean_ss","BO22_RCP60_2050_icethickmean_ss","BO22_RCP85_2050_icethickmean_ss", BO22_RCP26_2100_icethickmean_ss","BO22_RCP45_2100_icethickmean_ss","BO22_RCP60_2100_icethickmean_ss","BO22_RCP85_2100_icethickmean_ss"
ice.present.vars <- c("BO22_icethickmean_ss")

ice.rcp26<- "BO22_RCP26_2100_icethickmean_ss"
ice.rcp45<- "BO22_RCP45_2100_icethickmean_ss"
ice.rcp60<- "BO22_RCP60_2100_icethickmean_ss"
ice.rcp85<- "BO22_RCP85_2100_icethickmean_ss"

ice.present.rasters <- load_layers(ice.present.vars, datadir="./Rasters/")
ice.rcp26 <- load_layers(ice.rcp26, datadir="./Rasters/")
ice.rcp45 <- load_layers(ice.rcp45, datadir="./Rasters/")
ice.rcp60 <- load_layers(ice.rcp60, datadir="./Rasters/")
ice.rcp85 <- load_layers(ice.rcp85, datadir="./Rasters/")


#load temp data now that everything is at the same extent, and merge with ice rasters
temp.rast <- raster("./Rasters/MeanTemp_Current_StableClim.tif")
temp.rcP26 <- raster("../Environment//Rasters/rcp26_StableClim.tif")
temp.rcP45 <- raster("../Environment//Rasters/rcp45_StableClim.tif")
temp.rcP60 <- raster("../Environment//Rasters/rcp60_StableClim.tif")
temp.rcP85 <- raster("../Environment//Rasters/rcp85_StableClim.tif")

presenttemp = projectRaster(temp.rast, ice.present.rasters)
present = c(ice.present.rasters, presenttemp)
present.rasters <- raster::brick(present)

temp.rcp26 = projectRaster(temp.rcP26, ice.rcp26)
env.rcp26 = c(ice.rcp26, temp.rcp26)
rcp26 <- raster::brick(env.rcp26)

temp.rcp45 = projectRaster(temp.rcP45, ice.rcp45)
env.rcp45 = c(ice.rcp45, temp.rcp45)
rcp45 <- raster::brick(env.rcp45)

temp.rcp60 = projectRaster(temp.rcP60, ice.rcp60)
env.rcp60 = c(ice.rcp60, temp.rcp60)
rcp60 <- raster::brick(env.rcp60)

temp.rcp85 = projectRaster(temp.rcP85, ice.rcp85)
env.rcp85 = c(ice.rcp85, temp.rcp85)
rcp85 <- raster::brick(env.rcp85)

#Set the names in each raster stack so that they match with the original env input file in the gf analysis
names(present.rasters) <- c("IceThickness", "Temp")
names(rcp26) <- c("IceThickness","Temp")
names(rcp45) <- c("IceThickness","Temp")
names(rcp60) <- c("IceThickness","Temp")
names(rcp85) <- c("IceThickness","Temp")


## 2.2 Crop to desired area, using species' range here
# Load range shapefile (obtained from IUCN red list database)
species <- readOGR('../redlist_species_data_b48ae38b-6561-49b9-846f-033115599d7e/data_0.shp')

# Initial crop- probably not necessary 

crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent
crop.map.pres <- crop(present.rasters,crop.extent) #crop rasters
crop.rcp26 <- crop(rcp26, crop.extent)
crop.rcp45 <- crop(rcp45, crop.extent)
crop.rcp60 <- crop(rcp60, crop.extent)
crop.rcp85 <- crop(rcp85, crop.extent)


# Crop to species range shape file
present.rasters.range <- mask(crop.map.pres, species)
crop.rcp26.range <- mask(crop.rcp26, species)
crop.rcp45.range <- mask(crop.rcp45, species)
crop.rcp60.range <- mask(crop.rcp60, species)
crop.rcp85.range <- mask(crop.rcp85, species)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

#Plot
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plot(present.rasters.range, col=my.colors(200), add = TRUE)
plot(temp.rast)


## 2.3 Extract climate data from rasters and transform environmental predictors
# start with present data
clim.land <- raster::extract(present.rasters.range, 1:ncell(present.rasters.range), df = TRUE)
clim.land <- na.omit(clim.land)
row.names(clim.land) <- NULL
clim.land.noid <-(clim.land[,c("IceThickness", "Temp")])
clim.land.noid <-as.data.frame(clim.land.noid)
class(clim.land.noid) #dataframe

# use predict function to transform environmental predictors
pred <- predict(gf, clim.land.noid) 

# start with present data
# convert predictions to color scale for map
# use pca on predictions
pca <- prcomp(pred, center=T, scale.=F)
summary(pca)

# assign pcs to colors
r <- pca$x[, 1]
g <- pca$x[, 2]
b <- pca$x[, 2] #We are setting this to be the same as g because we have only 2 predictors. Might be a bit hacky but I couldn't find a neater alternative

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

#Plot 
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plotRGB(rgb.rast, bgalpha = 1, add = TRUE, col = "lightgrey")
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)


##########  
##### Step 3: Initial run with gradient forest on present data

#Only using the ice thickness gf model because we don't have forecasted ice cover in the dataset

#RCP26 - 2100
#Run each time for each section
clim.rcp26 <- raster::extract(crop.rcp26.range, 1:ncell(crop.rcp26.range), df = TRUE)
clim.rcp26 <- na.omit(clim.rcp26)
row.names(clim.rcp26) <- NULL
clim.rcp26.noid <-(clim.rcp26[,c("IceThickness","Temp")])
clim.rcp26.noid <-as.data.frame(dplyr::select(clim.rcp26, c("IceThickness","Temp")))

# transform environmental variables
pred.rcp26 <- predict(gf, clim.rcp26.noid)

# genetic offset
genetic.offset.adaptive <- sqrt((pred.rcp26[,1]-pred[,1])^2) + 
  (pred.rcp26[,2]-pred[,2])^2

# Define raster properties --- the variable doesn't matter here i think
offset26<- crop.rcp26.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offset26[clim.rcp26$ID] <- genetic.offset.adaptive

#Plot with color scale of choice

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

pdf("GeneticOffset_rcp26.pdf",  width = 6.78, height = 5.3)
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plot(offset26, col=my.colors(200), add = TRUE)
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()


#RCP45 - 2100
#Run each time for each section
clim.rcp45 <- raster::extract(crop.rcp45.range, 1:ncell(crop.rcp45.range), df = TRUE)
clim.rcp45 <- na.omit(clim.rcp45)
row.names(clim.rcp45) <- NULL
clim.rcp45.noid <-(clim.rcp45[,c("IceThickness","Temp")])
clim.rcp45.noid <-as.data.frame(dplyr::select(clim.rcp45, c("IceThickness","Temp")))

# transform environmental variables
pred.rcp45 <- predict(gf, clim.rcp45.noid)

# genetic offset:
genetic.offset.adaptive <- sqrt((pred.rcp45[,1]-pred[,1])^2) + 
  (pred.rcp45[,2]-pred[,2])^2

# Define raster properties --- the variable doesn't matter here i think
offset45<- crop.rcp45.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offset45[clim.rcp45$ID] <- genetic.offset.adaptive

#Plot with color scale of choice

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

pdf("GeneticOffset_rcp45.pdf",  width = 6.78, height = 5.3)
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plot(offset45, col=my.colors(200), add = TRUE)
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()



#RCP60 - 2100
#Run each time for each section
clim.rcp60 <- raster::extract(crop.rcp60.range, 1:ncell(crop.rcp60.range), df = TRUE)
clim.rcp60 <- na.omit(clim.rcp60)
row.names(clim.rcp60) <- NULL
clim.rcp60.noid <-(clim.rcp60[,c("IceThickness","Temp")])
clim.rcp60.noid <-as.data.frame(dplyr::select(clim.rcp60, c("IceThickness","Temp")))

# transform environmental variables
pred.rcp60 <- predict(gf, clim.rcp60.noid)

# genetic offset:
genetic.offset.adaptive <- sqrt((pred.rcp60[,1]-pred[,1])^2) + 
  (pred.rcp60[,2]-pred[,2])^2

# Define raster properties --- the variable doesn't matter here i think
offset60<- crop.rcp60.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offset60[clim.rcp60$ID] <- genetic.offset.adaptive

#Plot with color scale of choice

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

pdf("GeneticOffset_rcp60.pdf",  width = 6.78, height = 5.3)
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plot(offset60, col=my.colors(200), add = TRUE)
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()



#RCP85 - 2100
#Run each time for each section
clim.rcp85 <- raster::extract(crop.rcp85.range, 1:ncell(crop.rcp85.range), df = TRUE)
clim.rcp85 <- na.omit(clim.rcp85)
row.names(clim.rcp85) <- NULL
clim.rcp85.noid <-(clim.rcp85[,c("IceThickness","Temp")])
clim.rcp85.noid <-as.data.frame(dplyr::select(clim.rcp85, c("IceThickness","Temp")))

# transform environmental variables
pred.rcp85 <- predict(gf, clim.rcp85.noid)

# genetic offset:
genetic.offset.adaptive <- sqrt((pred.rcp85[,1]-pred[,1])^2) + 
  (pred.rcp85[,2]-pred[,2])^2

# Define raster properties --- the variable doesn't matter here i think
offset85<- crop.rcp85.range$IceThickness 

#Assign genetic offset values (difference between future and present predictions) to raster                
offset85[clim.rcp85$ID] <- genetic.offset.adaptive

#Plot with color scale of choice

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

pdf("GeneticOffset_rcp85.pdf",  width = 6.78, height = 5.3)
par(mar = c(5.1, 4.1, 4.1, 5))
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "2050 - 8.5 degree increase")
plot(offset85, col=my.colors(200), add = TRUE)
points(pres.clim.points$Lon, pres.clim.points$Lat, cex = 0.6)
dev.off()


######## extract gradient forest values
# save the rast.offset for later too
writeRaster(offset26, "offset.rcp26.2100.tif", overwrite = TRUE)
writeRaster(offset45, "offset.rcp45.2100.tif", overwrite = TRUE)
writeRaster(offset60, "offset.rcp60.2100.tif", overwrite = TRUE)
writeRaster(offset85, "offset.rcp85.2100.tif", overwrite = TRUE)

# can import raster back in
#rast.offset <- raster("offset.rcp26.2100.tif")

# treat as a raster to extract values by site.
# load subpop locations
colnames(pres.clim.points)
x <- pres.clim.points$Longitude
y <- pres.clim.points$Latitude

pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
projection(pts)==projection(offset26) # make sure projection matches

#Extract offset values and merge into dataframe
offset.data.26 <- as.data.frame(raster::extract(offset26, pts))
names(offset.data.26) <- "offset.rcp26"

offset.data.45 <- as.data.frame(raster::extract(offset45, pts))
names(offset.data.45) <- "offset.rcp45"

offset.data.60 <- as.data.frame(raster::extract(offset60, pts))
names(offset.data.60) <- "offset.rcp60"

offset.data.85 <- as.data.frame(raster::extract(offset85, pts))
names(offset.data.85) <- "offset.rcp85"

offset2 <- cbind(pres.clim.points, offset.data.26, offset.data.45, offset.data.60, offset.data.85)

na.check = map_int(offset1, ~sum(is.na(.)))

#save dataset
write.csv(offset2, "2024.04.16.GF.offset.all.rcp.csv")

