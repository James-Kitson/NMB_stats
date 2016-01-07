###########################################################################
############# Script for summarising MySeq reads  #################
###########################################################################

### Clear the workspace
rm(list=ls())

### load the libraries you need
library(akima)
library(plyr)
library(dplyr)

########################################################################################
################### summarisation of clustering parameters #############################
########################################################################################

### read in the data
my.clust<-read.csv("DATA/combined_reads_stats.csv", stringsAsFactors=FALSE)

### set a grouping variable for the data type
my.clust$type<-ifelse(grepl("DNA",my.clust$sample),"DNApositive",
                      ifelse(grepl("PCR",my.clust$sample),"PCRpositive",
                             ifelse(grepl("Negative",my.clust$sample),"negative","Moth sample")))

### subset the data to only include sample wells
my.clust.subs<-subset(my.clust, my.clust$type=="Moth sample",
                      select=c("cluster_thres", "clusters_min_cov", "cluster_above_thres"))


########################################################################################
### instructions for the next bit of code taken from
### http://stackoverflow.com/questions/19339296/plotting-contours-on-an-irregular-grid
########################################################################################

# interpolate data into a regulae grid of values for plotting - currently cluster numbers are logged to spread the colour ramp
fld <- with(my.clust.subs, interp(x = cluster_thres , y = clusters_min_cov, z = cluster_above_thres, duplicate="mean",
                                  xo=seq(min(cluster_thres), max(cluster_thres), length = 100),
                                  yo=seq(min(clusters_min_cov), max(clusters_min_cov), length = 100)))

### create a nice colour scale
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
### write this to an svg object
svg(filename="Clustering_contours_filtered_data.svg")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color=jet.colors,
               xlab = "Similarity threshold for clustering",
               ylab = "Minimum coverage of retained clusters",
               key.title = title(main = "Clusters\nretained", cex.main = 0.8))

dev.off()

