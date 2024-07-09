#Extracting marine data from Bio-ORACLE (https://tomjenkins.netlify.app/tutorials/r-extract-marine-data/ and Evelien code)
library(tidyverse)
library(sdmpredictors)
library(raster)
library(sp)
library(dismo)

setwd("~/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/Gradient_Forest/Environment/")

#Export a csv file containing marine variables of interest
# List marine data sets
datasets = list_datasets(terrestrial = FALSE, marine = TRUE)

# Variables of interest
variables = c("temp","ice")

# Extract present-day data sets
present = list_layers(datasets) %>%
  # select certain columns
  dplyr::select(dataset_code, layer_code, name, units, description, contains("cellsize"), version) %>%
  # keep variables of interest using a regular expression
  dplyr::filter(grepl(paste(variables, collapse = "|"), layer_code))

# Export present-day data sets to csv file
write.csv(present, file = "bio-oracle-present-datasets.csv")

# Future Representative Concentration Pathway (RCP) scenarios of interest
rcp = c("RCP26","RCP45","RCP60","RCP85")

# Extract future data sets
future = list_layers_future(datasets) %>%
  # keep RCP scenarios of interest
  dplyr::filter(grepl(paste(rcp, collapse = "|"), scenario)) %>% 
  # keep data for 2050 and 2100
  dplyr::filter(year == 2050 | year == 2100) %>% 
  # keep variables of interest using a regular expression
  dplyr::filter(grepl(paste(variables, collapse = "|"), layer_code))

# Export future data sets to csv file
write.csv(future, file = "bio-oracle-future-datasets.csv")


#2. Download and import Bio-ORACLE rasters ice data
ice.present = c("BO22_icethickmean_ss", "BO22_icecovermean_ss")

ice.future = c("BO22_RCP26_2050_icethickmean_ss","BO22_RCP45_2050_icethickmean_ss","BO22_RCP60_2050_icethickmean_ss","BO22_RCP85_2050_icethickmean_ss",
                "BO22_RCP26_2100_icethickmean_ss","BO22_RCP45_2100_icethickmean_ss","BO22_RCP60_2100_icethickmean_ss","BO22_RCP85_2100_icethickmean_ss"
)

# Combine present-day and future vectors
ice = c(ice.present, ice.future)

# Download rasters  and import into R
ice.raster = load_layers(ice)
ice.present.raster = load_layers("BO22_icethickmean_ss")


#4. Extract data from raster
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("chip12_ind_coord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(ice.raster)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
ice.data = tibble(ID = 1:nrow(pts@coords),
                     Lon = pts$x,
                     Lat = pts$y
)
ice.data

rasters <- raster::stack(ice.raster.crop)
nlayers(rasters)

#Create raster list with 20 km2 buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
, fun=mean)
}

#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
names(store_data) = names(rasters)
ice.data = bind_cols(ice.data, as_tibble(store_data))
ice.data
# Check each column for NA values
na.check = map_int(ice.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
ice.data.nona = ice.data %>% drop_na


#Visualise the spread of present-day sea temperature values for our points data set.

# Prepare a custom theme for ggplot
theme1 = theme(
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text = element_text(size = 9, face = "bold"),
  axis.title = element_text(size = 12, face = "bold")
)

# Violin plot and raw data
ice.data.nona %>% 
  # select only columns 5-8 (present-day sea temperature variables)
  dplyr::select(4:12) %>% 
  # transform data to long format for plotting
  pivot_longer(names_to = "Variable", values_to = "Values", cols = everything()) %>% 
  # plot data
  ggplot(data = .)+
  geom_violin(aes(x = Variable, y = Values, fill = Variable), show.legend = FALSE)+
  geom_jitter(aes(x = Variable, y = Values), show.legend = FALSE, alpha = 0.30)+
  scale_y_continuous(expand = c(0,0), limits = c(-1,6), breaks = c(seq(0,6,0.5)))+
  scale_fill_manual(values = heat.colors(9))+
  xlab("Ice Thickness (m)")+
  ylab(expression(bold("Ice Thickness (m)")))+
  coord_flip() +
  theme1

#Export data to a csv file and individual id coordinates
sites$ID <- c(1:1476)
head(sites)
ice_thickness_withcoord <- merge(sites, ice.data.nona)
head(ice_thickness_withcoord)
write_csv(ice_thickness_withcoord, path = "23.06.20_ice_thickness_cover_withcoord.csv")


#Create arctic zone map
library(rgdal)
my_spdf <- readOGR("ABA-Boundaries/Arctic_Zones_complete_polygons.shp")

summary(my_spdf)

library(broom)
spdf_fortified <- tidy(my_spdf,region = "Zone")

# Plot it
library(ggplot2)
ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group), fill="#69b3a2", color="white") +
  theme_void() 


#TEMPERATURE
# loads the packages used in this guide
library(ncdf4)

#Download StableClim Dataset (https://adelaide.figshare.com/articles/dataset/StableClim/12197976?file=24142217)
tsannual <- brick("StableClim_AnnMean_rcp85_ts.nc") #251 layers corresponding to 1850-2100 for rcp85 as a test

# plots the variable
plot(b)
plot(tsannual)
plot(tsmonthly)


# plot time points
plot(tsannual[[240:251]])

#Create new rasterbrick for the same time period as bio-oracle (2000-2014) and future environments (2040-2050 annual mean and 2090-2100 annual means). This is the same timeframe as sdmpredictors projections.

#first can we rename layers so they make sense
names(tsannual)
names(tsannual) <- c(1850:2100)
names(tsannual) #okay that worked well

#Subset current (2000-2014) so need to select the right layers (2000-1850 = 150 so start there and go to 164). Off by one so 151-165 it is

current <- tsannual[[151:166]]
plot(current)

#Summarize across years
meantemp <- mean(current)
plot(meantemp)


#Lets do the same for 2040-2050 and 2090-2100

rcp85.2050 <- tsannual[[191:201]]
plot(rcp85.2050)
rcp85.2050 <- mean(rcp85.2050)

rcp85.2100 <- tsannual[[241:251]]
plot(rcp85.2100)
rcp85.2100 <- mean(rcp85.2100)


#convert longitude from 0-360 to -180-180, already done for the StableCLim dataset
#raster_mean180 <- rotate(meantemp)

#save as new raster
writeRaster(meantemp, filename=file.path("MeanTemp_Current_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp85.2050, filename=file.path("MeanTemp_2050_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp85.2100, filename=file.path("MeanTemp_2100_StableClim.tif"), format="GTiff", overwrite=TRUE)

#Convert to stack

rastcurrent <- raster::stack(meantemp)
rast2050 <- raster::stack(rcp85.2050)
rast2100 <- raster::stack(rcp85.2100)

# Initial crop- probably not necessary but it works for this pipeline. can make more efficient later
crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent

#crop rasters
crop.temp <- crop(rastcurrent,crop.extent) 
crop.2050 <- crop(rast2050, crop.extent)
crop.2100 <- crop(rast2100, crop.extent)




#4. Extract data from raster for current conditions
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../Species_Dist_Modelling/Gradient_Forest/Environment/chip12_ind_coord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(crop.temp)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
temp.data = tibble(ID = 1:nrow(pts@coords),
                   Lon = pts$x,
                   Lat = pts$y
)
temp.data


#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
#Create raster list with buffer around each sample
temp = c(crop.temp, crop.2050, crop.2100)

rasters <- raster::stack(temp)
nlayers(rasters)

#Create raster list with 20 km2  buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
                                    , fun=mean)
}

names(store_data) = names(rasters)
temp.data = bind_cols(temp.data, as_tibble(store_data))
temp.data

# Check each column for NA values
na.check = map_int(temp.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
temp.data.nona = temp.data %>% drop_na


#Export data to a csv file and individual id coordinates
sites$ID <- c(1:1476)
head(sites)
temp_withcoord <- merge(sites, temp.data.nona)
head(temp_withcoord)
str(temp_withcoord)

# Rename columns and convert to celsius, not necessary for StableCLim data
library(dplyr)

temp_mean <- temp_withcoord %>% 
  # mutate(MeanTemp = layer.1 - 273.15, Temp2050 = layer.2 - 273.15, Temp2100 = layer.3 - 273.15) %>% 
  dplyr::rename( MeanTemp = layer.1 , Temp2050 = layer.2, Temp2100 = layer.3)
head(temp_mean)
str(temp_mean)

#write to csv
write_csv(temp_mean, "23.07.20_temp_withallcoord.csv")


#I added this in after performing an initial round of gradient forest analyses, and there is some redundancy with the earlier env data extractions

#Download other RCP scenarios
#Download StableClim Dataset (https://adelaide.figshare.com/articles/dataset/StableClim/12197976?file=24142217)
tsannual_26 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp26_ts.nc") #251 layers corresponding to 1850-2100
tsannual_45 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp45_ts.nc") #251 layers corresponding to 1850-2100
tsannual_60 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp60_ts.nc") #251 layers corresponding to 1850-2100
tsannual_85 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp85_ts.nc") #251 layers corresponding to 1850-2100


# plot time points
plot(tsannual[[240:251]])

#Create new rasterbrick for the same time period as bio-oracle (2000-2014) and future environments (2040-2050 annual mean and 2090-2100 annual means)

#first can we rename layers so they make sense
names(tsannual_26)
names(tsannual_26) <- c(1850:2100)
names(tsannual_45) <- c(1850:2100)
names(tsannual_60) <- c(1850:2100)
names(tsannual_85) <- c(1850:2100)


#Get mean value for 2090-2100

rcp26 <- tsannual_26[[241:251]]
plot(rcp26)
rcp26 <- mean(rcp26)

rcp45 <- tsannual_45[[241:251]]
rcp45 <- mean(rcp45)

rcp60 <- mean(tsannual_60[[241:251]])

rcp85 <- mean(tsannual_85[[241:251]])

#save as new raster
writeRaster(rcp26, filename=file.path("Rasters/rcp26_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp45, filename=file.path("Rasters/rcp45_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp60, filename=file.path("Rasters/rcp60_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp85, filename=file.path("Rasters/rcp85_StableClim.tif"), format="GTiff", overwrite=TRUE)

#Convert to stack

rastcurrent <- raster::stack(meantemp)
rast2050 <- raster::stack(rcp85.2050)
rast2100 <- raster::stack(rcp85.2100)

# Initial crop- probably not necessary but it works for this pipeline. can make more efficient later
crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent

#crop rasters
crop.rcp26 <- crop(rcp26,crop.extent) 
crop.rcp45 <- crop(rcp45, crop.extent)
crop.rcp60 <- crop(rcp60, crop.extent)
crop.rcp85 <- crop(rcp85, crop.extent)




#4. Extract data from raster for current conditions
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../GFAllEnvironments/23.11.01_tempandice_withcoord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(crop.temp)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
temp.data = tibble(ID = 1:nrow(pts@coords),
                   Lon = pts$x,
                   Lat = pts$y
)
temp.data


#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
#Create raster list with buffer around each sample
temp = c(crop.rcp26, crop.rcp45, crop.rcp60, crop.rcp85)
names(temp) = c("rcp26", "rcp45", "rcp60", "rcp85")
rasters <- raster::stack(temp)
nlayers(rasters)

#Create raster list with 20 km2 buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
                                    , fun=mean)
}

names(store_data) = names(rasters)
temp.data = bind_cols(temp.data, as_tibble(store_data))
temp.data

# Check each column for NA values
na.check = map_int(temp.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
temp.data.nona = temp.data %>% drop_na


#Export data to a csv file and individual id coordinates
sites$ID <- c(1:411)
head(sites)
temp_withcoord <- merge(sites, temp.data.nona)
head(temp_withcoord)
str(temp_withcoord)


#write to csv
write_csv(temp_withcoord, "../GFAllEnvironments/24.04.15_temp_withallcoord.csv")


#Do it again for  ice thickness (pull other scenarios)

ice.future = c("BO22_RCP26_2100_icethickmean_ss","BO22_RCP45_2100_icethickmean_ss","BO22_RCP60_2100_icethickmean_ss","BO22_RCP85_2100_icethickmean_ss"
)

# Combine present-day and future vectors

# Download rasters  and import into R
ice.raster = load_layers(ice.future)

# Crop rasters to boundary extent
crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent

ice.raster.crop = crop(ice.raster, crop.extent)
names(ice.raster.crop) <- c("ice.rcp26", "ice.rcp45", "ice.rcp60", "ice.rcp85")


#4. Extract data from raster
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../GFAllEnvironments/24.04.15_temp_withallcoord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(ice.raster)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
ice.data = tibble(ID = 1:nrow(pts@coords),
                  Lon = pts$x,
                  Lat = pts$y
)
ice.data

rasters <- raster::stack(ice.raster.crop)
nlayers(rasters)

#Create raster list with 20 km2 buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
                                    , fun=mean)
}

#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
names(store_data) = names(rasters)
ice.data = bind_cols(ice.data, as_tibble(store_data))
ice.data
# Check each column for NA values
na.check = map_int(ice.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
ice.data.nona = ice.data %>% drop_na



#Export data to a csv file and individual id coordinates
sites$ID <- c(1:411)
head(sites)
ice_thickness_withcoord <- merge(sites, ice.data.nona)
head(ice_thickness_withcoord)
write_csv(ice_thickness_withcoord, "../GFAllEnvironments/24.04.15_temp_withallcoord.csv")



#Make some plots
#3. Plot rasters
# Define colour scheme
cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Define a boundary
boundary = extent(c(xmin = -180, xmax = 20, ymin = 50, ymax = 90))


#load rasters
sites <- read.csv("../GFAllEnvironments/23.11.01_tempandice_withcoord.csv")
coord <- dplyr::select(sites, c("Longitude", "Latitude"))

ice.thick <- c("BO22_icethickmean_ss")
ice.thick <- load_layers(ice.thick, datadir="Rasters/")

ice.cov <- c("BO22_icecovermean_ss")
ice.cov <- load_layers(ice.cov, datadir="Rasters/")

temp <- raster("../Environment/MeanTemp_Current.tif")


# Crop rasters to boundary extent
it.crop = crop(ice.thick, boundary)
ic.crop = crop(ice.cov, boundary)
t.crop = crop(temp, boundary)
t.crop = t.crop -273.15

# Plot ice variables
par(mar = c(5.1, 4.1, 4.1, 5))

pdf("Ice Thickness.pdf",  width = 6.78, height = 5.3)
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice Thickness")
plot(it.crop, col=cols(200), add = TRUE)
points(coord, pch = 19, cex = .75, 
     xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F, fill = T,col = "lightgrey" )
dev.off()

pdf("Ice Cover.pdf",  width = 6.78, height = 5.3)
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice Thickness")
plot(ic.crop, col=cols(200), add = TRUE)
points(coord, pch = 19, cex = .75, 
       xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F, fill = T,col = "lightgrey" )
dev.off()

pdf("Temperature.pdf",  width = 6.78, height = 5.3)
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice Thickness")
plot(t.crop, col=cols(200), add = TRUE)
points(coord, pch = 19, cex = .75, 
       xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F, fill = T,col = "lightgrey" )
dev.off()
