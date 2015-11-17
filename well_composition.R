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
################### stacked bar charts of well composition   ###########################
########################################################################################

### read in the data
my.reads<-read.csv(file=paste("DATA/metaBEAT.tsv",sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
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

### set a minimum occurance for an assignment to be trusted
occurance=1

### subset the data frame to drop all assignments occuring fewer than the frequency specified above
my.reads.subs<-subset(my.reads, subset = rowSums(my.reads[2:ncol(my.reads)] > 0) > occurance)

### transpose the read data
my.reads.trans<-recast(my.reads.subs, variable~OTU_ID)
colnames(my.reads.trans)[1]<-"sample"
my.reads.trans$sample<-as.character(my.reads.trans$sample)

### use match to add the plate data to the read data
my.reads.trans$plate<-my.plates$plate[match(my.reads.trans$sample,my.plates$sample)]
my.reads.trans$plate.numeric<-my.plates$plate.numeric[match(my.reads.trans$sample,my.plates$sample)]

### use match to add the nest data to the read data
my.reads.trans$nest<-my.plates$nest[match(my.reads.trans$sample,my.plates$sample)]
my.reads.trans$nest.numeric<-my.plates$nest.numeric[match(my.reads.trans$sample,my.plates$sample)]

### Pull out the aissignments for further manual examination
my.assignments<-colnames(my.reads.trans)
write.csv(my.assignments, file="DATA/assignments_out.csv")

### Pull in the assigned colours as a .csv
Taxa.col<-read.csv("DATA/assignments_in.csv", stringsAsFactors=FALSE)

### drop the OPM NuMts and the higher order Thaumetopoea hits
#my.reads.trans<-my.reads.trans[,c(1:13,16:21)]

### combine the higer order hits
#my.reads.trans$higher.order.hits<-apply(my.reads.trans[,c(2,3,5,6,9:13)], MARGIN=1, FUN=sum)
### combine the positive hits
#my.reads.trans$positive.hits<-apply(my.reads.trans[,c(4,8,15)], MARGIN=1, FUN=sum)

### Drop all the hits we've just combined
#my.reads.trans<-my.reads.trans[,c(1,7,14,16:21)]

my.reads.trans$type<-ifelse(my.reads.trans$nest=="Negative","Negative",
                            ifelse(my.reads.trans$nest=="DNApositive","DNApositive",
                                   ifelse(my.reads.trans$nest=="PCRpositive","PCRpositive","Sample")))

### to remove the positive and negative samples for OTU counting in a barplot, I need to iteratively subset them out using the nest column.
### I've not found a more elegant way of doing this yet
my.reads.trans.samps.only<-my.reads.trans[my.reads.trans$nest != "Negative", ]
my.reads.trans.samps.only<-my.reads.trans.samps.only[my.reads.trans.samps.only$nest != "DNApositive", ]
my.reads.trans.samps.only<-my.reads.trans.samps.only[my.reads.trans.samps.only$nest != "PCRpositive", ]

### make sample and plate factors for faceting and ordering
my.reads.trans$sample<-as.factor(my.reads.trans$sample)
my.reads.trans$plate<-as.factor(my.reads.trans$plate)

### total all the reads
my.reads.trans$total<-rowSums(my.reads.trans[c(2:7)])
### calculate the percentage of reads in each well that are OPM
my.reads.trans$Percent.Thau<-(my.reads.trans$Thaumetopoea_processionea/my.reads.trans$total)*100
### work out a 20% of reads inclusion cutoff
#my.reads.trans$cutoff<-my.reads.trans$total*0.2

### replace all the species that are less than 20% of the total reads with zero
#my.reads.trans[,2:16][my.reads.trans[,2:16] < my.reads.trans$cutoff] <- 0

### subset the data frame to drop all coloumns containing only zeros
###my.reads.trans.drop<-cbind(my.reads.trans[,c(1, 34:39)], subset(my.reads.trans[,c(2:33)], select = colSums(my.reads.trans[,c(2:33)]) !=0))

### remove the total variable
my.reads.trans.drop<-my.reads.trans[,c(1:11,13)]

### reorder my.reads.melt by both percentage Thau and plate
#my.reads.melt<-my.reads.melt[order(my.reads.melt$Plate.numeric,-my.reads.melt$Percent.Thau),]
my.reads.trans.drop<-my.reads.trans.drop[order(-my.reads.trans.drop$Percent.Thau),]

### make a panel factor after setting the decreasing OPM order
my.reads.trans.drop$panel<-as.factor(c(rep(1,times=960/3),
                                       rep(2,times=960/3),
                                       rep(3,times=960/3)))

### melt the data into long format
my.reads.melt<-melt(my.reads.trans.drop, id.vars=c("sample","plate", "plate.numeric","Percent.Thau","nest","type","panel"))
colnames(my.reads.melt)<-c("Sample","Plate","Plate.numeric","Percent.Thau","Nest","Type","Panel","Species","Reads")

### create an ID variable for the order
my.reads.melt$level.order<-seq(1, length(my.reads.melt$Percent.Thau),1)

### order the species for plotting the legend
#my.reads.melt$Species <- factor(my.reads.melt$Species,
#                      levels=c("Thaumetopoea.processionea",
#                               "Triops.cancriformis",
#                               "Comaster.audax",
#                               "Astatotilapia.calliptera",
#                               "Bilateria",
#                               "Pancrustacea",
#                               "nohit"))

### order the samples from highest OPM % to lowest and in increasing plate number
my.reads.melt$Sample <- reorder(my.reads.melt$Sample, my.reads.melt$level.order)
### order the plates from lowest to highest
#my.reads.melt$Plate <- reorder(my.reads.melt$Plate, my.reads.melt$Plate.numeric)

my.reads.melt$Species <- factor(my.reads.melt$Species, Taxa.col$Assignment)

levels(my.reads.melt$Species)

### create a nice colour scale using colourRampPalette and RColorBrewer
### jet.colors3 <- colorRampPalette(brewer.pal(length(unique(my.reads.melt$Species)), "Spectral"))

### create a nice colour scale manually
colours<-c("#FFFF75",              # moth colour
           "#FFA375",              # Carcelia colour
           "#33AD5C",              # Bad hit colours
           "#3366FF",              # Positive colours
           "#000000")              # unassigned

vivid.colours<-c("#FFFF00",              # moth colour
                 "#FF3300",              # Carcelia colour
                 "#009933",              # Bad hit colours
                 "#0000FF",              # Positive colours
                 "#000000")              # unassigned

vivid.colours2<-as.character(Taxa.col$Colour)

### count the columns greater than zero and write to a new data frame
hit.hist<-data.frame(OTUs = rowSums(my.reads.trans.samps.only[c(2:7)] != 0), Type=my.reads.trans.samps.only$type)


########################################################################################
################### make a plot of % composition by PCR plate  ######################################
########################################################################################
svg(file="well_composition_reads_by_plate.svg", width=10, height=8)
### set up the ggplot
well.composition<-ggplot(data=my.reads.melt, aes(x=Sample, y=Reads, fill=Species)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(position="fill", stat="identity") +
  ### wrap the plot by plate  
  facet_wrap(~Panel, scales="free_x", nrow=3, ncol=1) +
  ### give it a percentage scale
  scale_y_continuous(labels = percent_format()) +
  ### set the colours
  # scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.reads.melt$Species)))) +
  scale_fill_manual(name=as.character(Taxa.col$Assignment), values = vivid.colours2) +
  ### add a sensible y axis label
  labs(y = "% of reads per well", x="PCR wells") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1), face="italic"),
        legend.title = element_text(size = rel(1)),
        strip.text.x = element_blank(),
        legend.position = "right")


### Make a ggplot object of our OTU counts
c<- ggplot(hit.hist, aes(factor(reorder(OTUs, -OTUs)))) +
  geom_bar(alpha=0.5) + labs(y = "Frequency", x="OTUs per well") +
  coord_flip() +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)))

grid.arrange(well.composition, c, heights=c(3/4, 1/4), ncol=1)

dev.off()