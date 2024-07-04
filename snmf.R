library(adegenet)
library(poppr)
library(LEA)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rworldmap)
library(rworldxtra)
library(ggsn)
library(sf)
library(raster)
library(rgeos)
library(maps)
library(maptools)
library(grid)
library(miscTools)
library(stringr)
library(ggpubr)
library(tidyverse)


setwd("/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/Gradient_Forest/Subset/")
ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "black"),
            axis.line = element_line(linewidth=1), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=15, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,face="plain",size=17),
            axis.title.x=element_text(vjust=0.1,face="plain",size=17),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15))

#snmf

# Create geno file
vcfin <- "Subset.ld.hwe.vcf"
vcf2lfmm(vcfin)

pop.info <- read.delim("Subset_IDs.txt", header = F)
names(pop.info) <- c("subpop", "IID")

genoin <- read.geno("Subset.ld.hwe.geno")
dim(genoin)
table(genoin)


##identify best clusters

obj.snmf = snmf(genoin, K = 1:13, ploidy = 2, entropy = T,
                alpha = 100, project = "new", repetitions = 10, seed = 42)
summary(obj.snmf)

obj.snmf1 = load.snmfProject("Subset.ld.hwe.snmfProject")

K <- summary(obj.snmf1)$crossEntropy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("temp") %>%
  mutate(K = as.numeric(str_extract(temp, "(\\d)+"))) %>%
  dplyr::select(-temp)

# choose omptimal (K = 3)
ggplot(K, aes(x = K, y = mean)) +
  geom_line(color = "black", size = 0.25 ) +
  geom_segment(aes(x = K, y = min, xend = K, yend = max)) +
  geom_point(shape = 21, size = 4, color = "black", fill = "blue") +
  scale_x_continuous(breaks = seq(0, 16, by = 1)) +
  labs(x = "Number of Clusters", y = "Cross-entropy criterion") +
  ng
ggsave("Cross-entropy.pdf", width = 6.78, height = 5.3, dpi=300)


#plot ancestry matrix
ce = cross.entropy(obj.snmf1, K = 5)
ce

lowest.ce = which.min(ce)
lowest.ce

qmatrix = as.data.frame(Q(obj.snmf1, K = 5, run = lowest.ce))
head(qmatrix)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

cols = my.colors(5)
pdf("IndividualsAncestry.pdf",  width = 6.78, height = 5.3)
barchart(obj.snmf1, K = 5, run = lowest.ce,
         border = NA, space = 0,
         col = cols,
         xlab = "Individuals",
         ylab = "Ancestry proportions K = 5",
         main = "Ancestry matrix") -> bp

#add column for order
pop.info$order <- c(1:411)
bp <- as.data.frame(bp)
bp <- left_join(bp, pop.info)

#create axis
axis(1, at = 1:length(bp$order), 
     labels = bp$subpop, las = 3, 
     cex.axis = .4)
dev.off()

# Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# Add individual IDs
qmatrix$Ind = pop.info$IID

# Add site IDs
qmatrix$Site = pop.info$subpop
head(qmatrix)

# Convert dataframe to long format

qlong = reshape2::melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

# Change order of sites by using the factor function
#site.order = c("BB","DS", "FB","GB","KB","LS","MC","NB","NW","SB","SH", "VM","WH")
qlong$Site= as.factor(qlong$Site)
levels(qlong$Site)

# Define colour palette
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

cols = my.colors(length(unique(qlong$variable)))

#Sort by cluster

admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
ng
admix.bar
ggsave("Ind_admixture_barplot_for merging.pdf", width = 6.78, height = 5.3, dpi=300)

# Plot admixture barplot 
admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Site, scales = "free", ncol = 2)+
  scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar

ggsave("Pop_admixture_barplot.pdf", width = 6.78, height = 5.3, dpi=300)


# ----------------- #
#
# Prepare pie charts
#
# ----------------- #

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = reshape2::melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

# Test function on one site
pie_charts(avg_admix, site = "FB", cols = cols)

# Apply function to all sites using for loop
site.ordered = sort(c("BB","DS", "FB","GB","KB","LS","MC","NB","NW","SB","SH", "VM","WH"))

pies = list()
for (i in site.ordered){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
}


# ----------------- #
#
# Prepare basemap
#
# ----------------- #

# Import txt file containing coordinates
coord = read.delim("Subset_Coord.txt", header = FALSE)
coord = cbind(coord, pop.info)
coord <- coord %>% 
  dplyr::rename("Lon" = "V1",
         "Lat" = "V2")
pop.coords <- coord  %>% 
  group_by(subpop)  %>% 
  dplyr::summarize(N = n(),
         Mean.Lon = mean(Lon),
         Mean.Lat = mean(Lat)
    )
pop.coords


# Order alphabetically by site]
coords = pop.coords[order(pop.coords$subpop), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$subpop)

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-141,-60, 50, 90) 
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
basemap = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="lightgrey",
               colour="black", size=0.5)+
  coord_quickmap(expand=F)+
  ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outline, dist = 500, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", 
                 location = "bottomleft", anchor = c(x = -78, y = 50.4),
                 st.bottom = FALSE, st.size = 3, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  ng
basemap



# ----------------- #
#
# Add pie charts to basemap
#
# ----------------- #

# Extract coordinates for each site
coord.list = list()


for (i in site.ordered){
  coord.list[[i]] = c(subset(coords, subpop == i)$Mean.Lon, subset(coords, subpop == i)$Mean.Lat)
}
coord.list

# Define pie chart sizes
radius = 3.5

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(site.ordered)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
pie.map = basemap + pies.ac
pie.map
ggsave("Admixture_pie_charts_map.pdf", width = 6.78, height = 5.3, dpi = 300)



#Assign Individual cluster based on 75% ancestry proportion
names(qmatrix) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Ind", "Site")
cluster.assign <- qmatrix %>% 
  mutate(
    Assigned = case_when(
      Cluster1 > .75 ~ "Cluster 1",
      Cluster2 > .75 ~ "Cluster 2",
      Cluster3 > .75 ~ "Cluster 3",
      Cluster4 > .75 ~ "Cluster 4",
      Cluster5 > .75 ~ "Cluster 5",
      .default = "Admixed"
    )
  )

write_delim(cluster.assign, "Subset_ld_cluster_membership.txt")
