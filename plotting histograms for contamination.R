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
### process the data
###########################################################################

### transpose the read data
my.reads.trans<-recast(my.reads, variable~OTU_ID)
colnames(my.reads.trans)[1]<-"sample"
my.reads.trans$sample<-as.character(my.reads.trans$sample)

### use match to add the type data to the read data
my.reads.trans$type<-my.plates$type[match(my.reads.trans$sample,my.plates$sample)]
my.reads.trans$type.numeric<-my.plates$type.numeric[match(my.reads.trans$sample,my.plates$sample)]

### subset the data to only contain the sampel data and dropt the type variable
my.reads.trans.subs<-subset(my.reads.trans, my.reads.trans$type=="Sample", select=-type)

### melt the data into long format
my.reads.trans.melt<-melt(my.reads.trans.subs, id.vars="sample",variable.name="species", value.name="reads")

###########################################################################
### Plot histograms of each assignment read depth panelled by assignment
###########################################################################

### create a nice colour scale using colourRampPalette and RColorBrewer
colors <-brewer.pal(length(unique(my.reads.trans.melt$species)), "Set3")

############## plotting order highest to lowest threshold - hashed out to only select preferred order
### make the ggplot object
hl<-ggplot(my.reads.trans.melt, aes(x=reads))
### add the kernal density plots
hl + geom_histogram(aes(fill=factor(species)), alpha=1) +
  ### give the legend a sensible name
  scale_fill_discrete(name="BLAST assignment") +
  ### separate by clustering similarity threshold
  facet_wrap(~species, scales="free", ncol=2) +
  ### set nice axis labels
  labs(x = "Read Depth", y = "Frequency")

### save the graph to an svg plot
ggsave(filename="Read_depth_distribution_histogram.svg")