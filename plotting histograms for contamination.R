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

### Plot