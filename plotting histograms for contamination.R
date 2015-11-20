###########################################################################
###### Script for plotting the distributions of read depths and % reads for each assignment
##### this is to look for potential cut-off values to clean my data
###########################################################################

### Clear the workspace
rm(list=ls())

### load the libraries you need
library(reshape2)
library(ggplot2)
library(akima)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scales)

### read in the data
my.reads<-read.csv(file=paste("DATA/metaBEAT.tsv",sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

### transpose the read data
my.reads.trans<-recast(my.reads, variable~OTU_ID)
colnames(my.reads.trans)[1]<-"sample"
my.reads.trans$sample<-as.character(my.reads.trans$sample)

### melt the data into long format
my.reads.trans.melt<-melt(my.reads.trans, id.vars="sample",variable.name="species", value.name="reads")

###########################################################################
### Plot histograms of each assignment read depth panelled by assignment
###########################################################################

### create a nice colour scale using colourRampPalette and RColorBrewer
colors <-brewer.pal(length(unique(my.reads.trans.melt$species)), "Set3")

############## plotting order highest to lowest threshold - hashed out to only select preferred order
### make the ggplot object
hl<-ggplot(my.reads.trans.melt, aes(x=reads))
### add the kernal density plots
hl + geom_density(aes(fill=factor(species)), alpha=0.5) +
  ### give the legend a sensible name
  scale_fill_discrete(name="BLAST assignment") +
  ### separate by clustering similarity threshold
  facet_wrap(~species, scales="free", ncol=2) +
  ### set nice axis labels
  labs(x = "Read Depth", y = "Kernel density")

############## plotting order lowest to highest threshold
### make the ggplot object
lh<-ggplot(my.clust2, aes(x=cluster_above_thres))
### add the kernal density plots
lh + geom_density(aes(fill=rev(factor(clusters_min_cov))), alpha=0.5) +
  ### colour with our colour scale
  scale_fill_manual(name="Minimum number of\nsequences in cluster",
                    values = jet.colors2(length(unique(my.clust2$clusters_min_cov)))) +
  ### separate by clustering similarity threshold
  facet_wrap(~cluster_thres, scales="free", ncol=2) +
  ### set nice axis labels
  labs(x = "Clusters retained", y = "Kernel density")

### save the graph to an svg plot
ggsave(filename="clustering_kernel_density.svg")