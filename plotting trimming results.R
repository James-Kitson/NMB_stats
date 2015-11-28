###########################################################################
###### Script for summarising MySeq reads into per well composition  #####
###########################################################################

### Clear the workspace
rm(list=ls())

### load the libraries you need
library(reshape2)
library(ggplot2)
library(akima)
library(plyr)
library(dplyr)
library(ggplot2bdc)
library(ggExtra)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scales)

########################################################################################
################### box plots of read depth   ###########################
########################################################################################

### read in the data
my.reads<-read.csv(file="DATA/reads_stats.csv", stringsAsFactors=FALSE, header=TRUE)
#colnames(my.reads)

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

### use match to add the plate data to the read data
my.reads$plate<-my.plates$plate[match(my.reads$sample,my.plates$sample)]
my.reads$plate.numeric<-my.plates$plate.numeric[match(my.reads$sample,my.plates$sample)]

### set the plotting order of the plates
my.reads$plate<-factor(reorder(my.reads$plate, my.reads$plate.numeric))

### make sample and plate factors for faceting and ordering
my.reads$sample<-as.factor(my.reads$sample)
my.reads$plate<-as.factor(my.reads$plate)

### set a grouping variable for the data type
my.reads$type<-ifelse(grepl("DNA",my.reads$sample),"DNApositive",
                      ifelse(grepl("PCR",my.reads$sample),"PCRpositive",
                             ifelse(grepl("Negative",my.reads$sample),"negative","Moth sample")))

### subset the data to only keep the numbers of reads at each stage
my.reads.subs<-subset(my.reads, select=c("sample","total","trimmed.total","plate","type"))

### melt the data into long format
my.reads.subs.melt<-melt(my.reads.subs, id.vars=c("sample","type", "plate"))

########################################################################################
################### Plot the data including the +ves and -ves  ###########################
########################################################################################

### make the ggplot object + add the jittered dots + add the boxplots and colour them by trim level and make them a bit transparent
ggplot(aes(y = value, x = plate, fill = variable), data = my.reads.subs.melt) +
  ### make the boxplot and suppress the outliers as we are plotting the points anyway
  geom_boxplot(aes(fill=variable), alpha=0.5, position = position_dodge(width = 0.85), outlier.shape = NA) +
  ### plot the points
  geom_point(pch = 21, position = position_jitterdodge()) +
  ### fix the axes titles
  labs(y = "Reads per PCR well", x="PCR plate") +
  scale_fill_discrete(name="Trim level",
                      breaks=c("total", "trimmed.total"),
                      labels=c("Raw reads", "Trimmed reads")) +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(size = rel(1.9)),
        axis.text.y = element_text(size = rel(1.9)),
        axis.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.9)),
        legend.title = element_text(size = rel(1.9)),
        legend.position = "right")

### save the graph to an svg plot
ggsave(filename="filtered_trimming_summary_all_samples.svg")

########################################################################################
################### Plot the data excluding the +ves and -ves  ###########################
########################################################################################

### make the ggplot object + add the jittered dots + add the boxplots and colour them by trim level and make them a bit transparent
ggplot(aes(y = value, x = plate, fill = variable), data = subset(my.reads.subs.melt, my.reads.subs.melt$type=="Moth sample")) +
  ### make the boxplot and suppress the outliers as we are plotting the points anyway
  geom_boxplot(aes(fill=variable), alpha=0.5, position = position_dodge(width = 0.85), outlier.shape = NA) +
  ### plot the points
  geom_point(pch = 21, position = position_jitterdodge()) +
  ### fix the axes titles
  labs(y = "Reads per PCR well", x="PCR plate") +
  scale_fill_discrete(name="Trim level",
                      breaks=c("total", "trimmed.total"),
                      labels=c("Raw reads", "Trimmed reads")) +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(size = rel(1.9)),
        axis.text.y = element_text(size = rel(1.9)),
        axis.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.9)),
        legend.title = element_text(size = rel(1.9)),
        legend.position = "right")

### save the graph to an svg plot
ggsave(filename="filtered_trimming_summary_moths_only.svg")

####################################################################################################################################
################################## counting reads etc ##############################################################################
####################################################################################################################################

### totals and counts for text
counting<-subset(my.reads, my.reads$type=="Moth sample")

max(counting$total)
min(counting$total)

mean(counting$total)
sd(counting$total)

max(counting$trimmed.total)
min(counting$trimmed.total)

mean(counting$trimmed.total)
sd(counting$trimmed.total)

mean(counting$cluster_above_thres)
sd(counting$cluster_above_thres)

sum(my.reads$total)
sum(my.reads$trimmed.total)
