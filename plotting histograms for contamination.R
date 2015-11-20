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

###########################################################################
### read in all the data and sort out the samples for labelleing and subsetting
###########################################################################
### read in the data
my.reads<-read.csv(file=paste("DATA/metaBEAT.tsv",sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

### read in the sample by plate data
my.plates<-read.csv(file="DATA/Samples_and_MIDS_corrected.txt",
                    sep="\t", stringsAsFactors=FALSE, header=TRUE)

### trim the plate data to the necessary columns
my.plates<-my.plates[,1:3]

### greedy regex to split the sample string by the underscore and leave us with a column for nest and a column for indentifier
my.plates<-cbind(my.plates, do.call(rbind, strsplit(as.character(my.plates$sample), "_|_.*_")))

### trim the plate data to the necessary columns (i.e. drop the identifier column)
my.plates<-my.plates[,1:4]
### name the columns
colnames(my.plates)<-c("sample","plate","plate.numeric","nest")

my.plates$type<-ifelse(my.plates$nest=="Negative","Negative",
                            ifelse(my.plates$nest=="DNApositive","DNApositive",
                                   ifelse(my.plates$nest=="PCRpositive","PCRpositive","Sample")))

###########################################################################
### process the data as read depths
###########################################################################

### transpose the read data
my.reads.trans<-recast(my.reads, variable~OTU_ID)
colnames(my.reads.trans)[1]<-"sample"
my.reads.trans$sample<-as.character(my.reads.trans$sample)

### use match to add the type data to the read data
my.reads.trans$type<-my.plates$type[match(my.reads.trans$sample,my.plates$sample)]
my.reads.trans$type.numeric<-my.plates$type.numeric[match(my.reads.trans$sample,my.plates$sample)]

### use match to add the plate.numeric data to the read data
my.reads.trans$plate.numeric<-my.plates$plate.numeric[match(my.reads.trans$sample,my.plates$sample)]
my.reads.trans$plate.numeric.numeric<-my.plates$plate.numeric.numeric[match(my.reads.trans$sample,my.plates$sample)]

### subset the data to only contain the sample data and drop the type variable
my.reads.trans.subs<-subset(my.reads.trans, my.reads.trans$type=="Sample", select=-type)

### calculate a total read depth for each sample
my.reads.trans.subs$total<-rowSums(my.reads.trans.subs[c(2:7)])

### calculate the % of each well that each assignment makes
my.reads.trans.subs[,10:15]<-my.reads.trans.subs[,2:7]/my.reads.trans.subs$total*100

### label the columns appropriately
colnames(my.reads.trans.subs)<-gsub(".1", "_percentage", colnames(my.reads.trans.subs))

### drop the total variable as we don't need it in the stacked data
my.reads.trans.subs<-my.reads.trans.subs[,c(1:8,10:15)]

### calculate the maximum % of each sample that any +ve assignment makes
my.reads.trans.subs$cutoff<-apply(my.reads.trans.subs[, c(2,4,6)], 1, max)

### replace all the species that are less than 20% of the total reads with zero
#my.reads.trans[,2:16][my.reads.trans[,2:16] < my.reads.trans$cutoff] <- 0

### melt the data twice on each of reads and percentages
read.depth<-melt(my.reads.trans.subs[,1:8], id.vars=c("sample","plate.numeric"),variable.name="species", value.name="reads")
percentages<-melt(my.reads.trans.subs[,c(1,8:14)], id.vars=c("sample","plate.numeric"),variable.name="species", value.name="percentage")

### bind all this together to give the full long format data
my.reads.trans.melt<-cbind(read.depth,percentages$percentage)
colnames(my.reads.trans.melt)<-c("sample","plate.numeric","species","reads","percentage")

### create a nice colour scale using colourRampPalette and RColorBrewer
colors <-brewer.pal(length(unique(my.reads.trans.melt$species)), "Set3")

###########################################################################
### Plot histograms of each assignment read depth panelled by assignment
###########################################################################

### make the ggplot object
reads.plot<-ggplot(my.reads.trans.melt, aes(x=reads))
### add the kernal density plots
reads.plot + geom_histogram(aes(fill=factor(species)), alpha=1, binwidth=30, origin=-15) +
  ### give the legend a sensible name
  scale_fill_discrete(name="BLAST assignment") +
  ### separate by clustering similarity threshold
  facet_wrap(~species, scales="free", ncol=2) +
  ### set nice axis labels
  labs(x = "Read Depth", y = "Frequency")

### save the graph to an svg plot
ggsave(filename="Diagrams/Read_depth_distribution_histogram.svg")

###########################################################################
### Plot the data as % of each well for sample wells only
###########################################################################

### make the ggplot object
percentage.plot<-ggplot(my.reads.trans.melt, aes(x=percentage))
### add the kernal density plots
percentage.plot + geom_histogram(aes(fill=factor(species)), alpha=1,binwidth=5, origin=-2.5) +
  ### give the legend a sensible name
  scale_fill_discrete(name="BLAST assignment") +
  ### separate by clustering similarity threshold
  facet_wrap(~species, scales="free", ncol=2) +
  ### set nice axis labels
  labs(x = "Percentage of reads in well", y = "Frequency")

### save the graph to an svg plot
ggsave(filename="Diagrams/Percentage_depth_distribution_histogram.svg")
  
  