#!/usr/bin/env Rscript
#SBATCH -J conStructxval2
#SBATCH -o CSxval2.out
#SBATCH -e CSxval2.err
#SBATCH -n 10
#SBATCH -t 9000
#SBATCH --mem=300G
#SBATCH --partition=standard
#SBATCH --account=grdi_genarcc
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04,tmpfs_size=20G"

.libPaths(new = "/gpfs/fs7/grdi/genarcc/wp3/PolarBears/env/R")
setwd("/gpfs/fs7/grdi/genarcc/wp3/PolarBears/Offsets/")

#Load genetic data
library(conStruct)
library(parallel)
library(foreach)
library(doParallel)

load("ld.myConStructData.RData")

#Load coordinates
coord = read.delim("Subset_Coord.txt", header = FALSE)
names(coord) <- c("Lon", "Lat")
coordinates <- as.matrix(coord)

#Geographic distance
#install.packages("geosphere")
library(geosphere)

dist <- distm(coordinates)

cl <- makeCluster(10,type="FORK")
registerDoParallel(cl)

#Run cross validation for 1-10 layers with 8 reps per layer
my.xvals <- x.validation(train.prop = 0.85,
                         n.reps = 8,
                         K = 1:10,
                         freqs = freqs,
                         data.partitions = NULL,
                         geoDist = dist,
                         coords = coordinates,
                         prefix = "xval2",
                         n.iter = 1e3,
                         make.figs = FALSE,
                         save.files = FALSE,
                         parallel = TRUE,
                         n.nodes = 10)


stopCluster(cl)