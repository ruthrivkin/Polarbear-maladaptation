#!/usr/bin/env Rscript
#SBATCH -J conStruct10
#SBATCH -o CSk10.out
#SBATCH -e CSk10.err
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=200G
#SBATCH --partition=standard
#SBATCH --account=grdi_genarcc
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04,tmpfs_size=20G"

.libPaths(new = "/gpfs/fs7/grdi/genarcc/wp3/PolarBears/env/R")
setwd("/gpfs/fs7/grdi/genarcc/wp3/PolarBears/Offsets/")

#Load genetic data
library(conStruct)
load("ld.myConStructData.RData")

#Load coordinates
coord = read.delim("Subset_Coord.txt", header = FALSE)
names(coord) <- c("Lon", "Lat")
coordinates <- as.matrix(coord)

#Geographic distance
#install.packages("geosphere")
library(geosphere)

dist <- distm(coordinates)

# run a conStruct analysis

#   you have to specify:
#       the number of layers (K)
#       the allele frequency data (freqs)
#       the geographic distance matrix (geoDist)
#       the sampling coordinates (coords)

#Change value of K (1-10) to run analysis on other layers
my.run.k10 <- conStruct(spatial = TRUE, 
                    K = 10, 
                    freqs = freqs,
                    geoDist = dist, 
                    coords = coordinates,
                    prefix = "pb.spatial.k10",
                    make.figs = TRUE)

