setwd("~/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/Gradient_Forest/GFAllEnvironments/")
library(ggplot2)
library(tidyverse)

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "black"),
            axis.line = element_line(linewidth=1), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=12, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,face="plain",size=15),
            axis.title.x=element_text(vjust=0.1,face="plain",size=15),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15))

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# legend.position = "top", legend.direction="horizontal", 
#legend.text=element_text(size=12), legend.key = element_rect(fill = "white"), 
#legend.title = element_text(size=12),legend.key.size = unit(1, "cm"))

#load dataset 
rcp.data <- read.csv("2024.04.16.GF.offset.all.rcp.csv")
colnames(rcp.data)
dim(rcp.data)

#Add  chip number and sampling year to the dataset
chip2 <- read.csv('/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/ConGen_Revision/Revision2024/Files to Submit/Table S1.csv')
head(chip2)
names(chip2)[3] <- "IID"
names(chip2)[4] <- "Sample.Year"
chip2$Chip <- "Chip2"
chip2 <- dplyr::select(chip2, c("IID", "Sample.Year", "Chip"))
dim(chip2)

chip1 <- read.csv("/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/SourceFiles/Chip1SampleData.csv")
head(chip1)
dim(chip1)
names(chip1)[1] <- "IID"
names(chip1)[4] <- "Sample.Year"
chip1$Chip <- "Chip1"
chip1 <- dplyr::select(chip1, c("IID", "Sample.Year", "Chip"))

chips <- rbind(chip1,chip2)
dim(chips)

match <- left_join(rcp.data, chips)
dim(match)
head(match)


#Add population census sizes (https://www.iucn-pbsg.org/wp-content/uploads/2021/11/July-2021-Status-Report-Web.pdf)
##Status Report on the World’s Polar Bear Subpopulations, IUCN/SSC Polar Bear Specialist Group, July 2021

census <- read.csv("CensusSize.csv")
str(census)

#Add cluster membership data (labeled as Assigned, idk why I picked that name instead of cluster but I'm too lazy to fix it)
cluster <- read.delim("../Subset/Subset_ld_cluster_membership.txt", sep ="")
str(cluster)
cluster <- dplyr::select(cluster, "Ind", "Site", "Assigned")
colnames(cluster) <- c("IID", "subpop", "Assigned")

#Add census size to dataset
data.census <- inner_join(match, census)
str(data.census)

data.all <- inner_join(data.census, cluster)
str(data.all)

#Set subpopulation location in high, low, sub Arctic
data.all$location <- case_when(data.all$subpop %in% "SH" ~ "sub",
                               data.all$subpop %in% c("WH", "FB","DS","GB") ~ "low",
                               data.all$subpop %in% c("BB","KB","LS","NB","SB", "MC","VM","NW") ~ "high"
)


write.csv(data.all,"GF.ld.offset.rcp.compare.csv")

#Load data
data.all <- read.csv("GF.ld.offset.rcp.compare.csv")
str(data.all)


#Make a single column for RCP scenario

datalong <- dplyr::select(data.all, c("Sample.Year", "Chip","Latitude", "IID", "location", "Assigned" ,"offset.rcp26", "offset.rcp45", "offset.rcp60", "offset.rcp85", "subpop", "Population.Size"))
datalong = reshape2::melt(datalong, id.vars=c("Sample.Year","Chip","Latitude", "Assigned", "subpop", "IID", "Population.Size", "location"),
                                      variable.name = "RCP", 
                                      value.name = "Offset")
str(datalong)

write.csv(datalong, "PopulationandClusterOffset.ld.rcp.compare.csv")

datalong$location <- as.factor(datalong$location)
datalong$Assigned <- as.factor(datalong$Assigned)

str(datalong)
dim(datalong)


#Run all rcps in a single model.. didn't use this model in the paper
library(car)
library(stats)
climate.scenarios.sq.cluster <- glm(sqrt(Offset) ~RCP + location + Assigned  + Population.Size +  Chip + Sample.Year,
                                    family = quasibinomial,
                                    data = datalong)

summary(climate.scenarios.sq.cluster)
aov.cluster <- Anova(climate.scenarios.sq.cluster, type = 2)
aov.cluster

anova(climate.scenarios.sq.cluster, test = "F")

TukeyHSD(aov(climate.scenarios.sq.cluster), which = "location" )

#Run each scenario separately
rcp26.glm<- glm(offset.rcp26 ~ location + Assigned  + Population.Size+ Chip + Sample.Year,
                                    family = quasibinomial,
                                    data = data.all)

summary(rcp26.glm)
aov.cluster <- Anova(rcp26.glm)
aov.cluster
anova(rcp26.glm, test = "F")

rcp45.glm<- glm(offset.rcp45 ~ location + Assigned  + Population.Size + Chip + Sample.Year,
                family = quasibinomial,
                data = data.all)

summary(rcp45.glm)
aov.cluster <- Anova(rcp45.glm)
aov.cluster
anova(rcp45.glm, test = "F")

rcp60.glm<- glm(offset.rcp60 ~ location + Assigned  + Population.Size+ Chip + Sample.Year,
                family = quasibinomial,
                data = data.all)

summary(rcp60.glm)
aov.cluster <- Anova(rcp60.glm)
aov.cluster
anova(rcp60.glm, test = "F")


rcp85.glm<- glm(offset.rcp85 ~ location + Assigned  + Population.Size+ Chip + Sample.Year,
                family = quasibinomial,
                data = data.all)

summary(rcp85.glm)
aov.cluster <- Anova(rcp85.glm)
aov.cluster

anova(rcp85.glm, test = "F")



#Calculate means for reporting
meancl <- Rmisc::summarySE(datalong, measurevar = "Offset", groupvars = "Assigned", na.rm = TRUE)
meancl

#Assess fit of model (did this for all models, just change the model name)
plot(rcp85.glm)
hist((data.all$offset.rcp85), xlab = "Response", main = "") #looks okay
vif(rcp85.glm) #Less than 10 so no issues

res <- leveneTest(sqrt(Offset) ~ location, data = datalong) #not great but okay
res

#Plot some stuff
cols = my.colors(length(unique(data.all$Assigned)))
cols = c("grey","#EDF0C0" ,"#5E85B8", "#A5BABC"  ,"#D79073" ,"#C13127")

(offsetplot.subpop <- ggplot(data.all, aes(x=Assigned, y=offset.2100))+ 
    geom_boxplot(outlier.size = 1, color = cols ) +
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("Cluster") +
    ylab("Genetic Offset")+
    ng
)

ggsave("Future_Offset2100_by_cluster.pdf", width = 6.78, height = 5.3, dpi = 300)

cols = my.colors(length(unique(data.all$subpop)))
cols = c("#EDF0C0" ,"#5E85B8", "#A5BABC"  ,"#D79073" ,"#C13127")

site.ordered = sort(c("SH","WH", "FB","DS","SB","GB","BB","NB","KB","MC","LS", "VM","NW"))

(offsetplot.subpop <- data.all %>%
    mutate(subpop = fct_relevel(subpop, "SH","WH", "FB","DS","GB","BB","SB","NB","KB","MC","LS", "VM","NW")) %>%
    ggplot(aes(x=subpop, y=offset.rcp85)) + 
    geom_boxplot(outlier.size = 1, fill = cols) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("Subpopulation") +
    ylab("Genetic Offset")+
ng
)

ggsave("Future_Offset2100_by_subpop.pdf", width = 6.78, height = 5.3, dpi = 300)


#Plot by population size
(offsetplot.census <- ggplot(data.all, aes(x=Sample.Year, y=offset.rcp85, col = subpop)) + 
    geom_point(alpha = 0.5) + 
    stat_smooth(method = "glm", col = "black") +
    xlab("Subpopulution Abundance") +
    ylab("Genetic Offset")+
    ng
)

ggsave("Future_Offset2100_by_size.pdf", width = 6.78, height = 5.3, dpi = 300)


#Plot by rcp*location
cols = my.colors(length(unique(datalong$RCP)))

(offsetplot.rcp <- datalong %>%
    mutate(subpop = fct_relevel(location, "sub","low", "high")) %>%
    ggplot(aes(x=subpop, y=Offset, fill = RCP)) + 
    geom_boxplot(outlier.size = 1) + 
#    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    scale_fill_manual(values = cols) +
    geom_jitter(shape=1, position=position_jitter(0.2), alpha = 0.3) + 
    xlab("Subpopulation") +
    ylab("Genetic Offset")+
    ng
)

ggsave("Offset_by_RCP.pdf", width = 6.78, height = 5.3, dpi = 300)

#Plot by rcp*latitude
cols = my.colors(length(unique(datalong$RCP)))

(offsetplot.lat <- ggplot(datalong, aes(x=Latitude, y=Offset, col = RCP)) + 
    geom_point(alpha = 0.3) + 
    geom_smooth(method = "glm") +
   scale_color_manual(values = cols) +
    xlab("Latitude") +
    ylab("Genetic Offset")+
    ng
)

ggsave("Future_Offset_by_lat.pdf", width = 6.78, height = 5.3, dpi = 300)



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


# Plot basemap
(localsmap = ggplot()+
    geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="lightgrey",
                 colour="black")+
    geom_point( data=data.all, aes(x=Longitude, y=Latitude, col = offset.2100, shape = Assigned), cex = 1.5) +
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    scale_shape_manual(values = c(3, 1,2,5,6,4)) +
    coord_quickmap(expand=F)+
    ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
    xlab("Longitude (°W)")+
    ylab("Latitude (°N)")+
    
    ng
)

ggsave("Map_points_coded_by_offset.pdf", width = 6.78, height = 5.3, dpi = 300)




#Correlation between pedigree 

pedigree <- read.csv("../../../Heterozygosity Correlates of Fitness/Source Files/WH_pedigree_cubcounts.csv")
str(pedigree)
ped <- pedigree %>% dplyr::select(SampleName, Cubs, Sex) 
colnames(ped)[1] ="IID"

#Condensed the capture data file so that each line is a single individual, did this as part of another project

#Check overlap between fitness and Offset datasets
offset.fitness <- merge(ped, data.all, by = "IID") #50 nice

cor.test(offset.fitness$offset.2100, offset.fitness$Cubs,  method = "kendall") #Reported this

lh.test <- glm(Cubs ~ offset.2100, data = offset.fitness, family = "poisson") #Same result
summary(lh.test)
anova(lh.test, test = "Chisq")

#plot
plot(offset.fitness$offset.2100, offset.fitness$Cubs)
abline(lm(offset.fitness$offset.2100 ~ offset.fitness$Cub))

par()

(offsetplot.lh <- ggplot(offset.fitness, aes(x=offset.2100, y=Cubs,)) + 
    geom_point(alpha = 1) + 
    stat_smooth(method = "glm", col = "black") +
    xlab("Genetic Offset Score") +
    ylab("Number of Cubs per Lifetime")+
    ng
)

ggsave("LifetimeReprodSuccess.pdf", width = 6.78, height = 5.3, dpi = 300)


#Check relatedness
ped2 <- read.csv("../../../Heterozygosity Correlates of Fitness/Source Files/bearPED.csv")
str(ped2)
colnames(ped2)[1] ="IID"

relatedness <- inner_join(ped2, data.all)

#Found 6 parent-offspring pairs, remove from analysis and see if theres is a difference, no difference

nopar <- datalong %>%
  filter(IID != "X00520",
         IID != "X04188", 
         IID != "X09579",
         IID != "X10997",
         IID != "X13018",
         IID != "X12219")

nopar.analysis <- glm(Offset ~ Year + Assigned  + location + Population.Size,
                                    family = quasibinomial,
                                    data = nopar)

summary(nopar.analysis)
aov.cluster <- aov(nopar.analysis)
summary(aov.cluster)
tukey.cluster <- TukeyHSD(aov.cluster)


