###########################################################################
############# Script for summarising MySeq reads  #################
###########################################################################

### Clear the workspace
rm(list=ls())

### set working directory
setwd("/Volumes/JK GEL/OPM_run1_filtered_cluster_test/")

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


### read in your data
my.data<-read.csv("reads_stats.csv")

### set a grouping variable for the data type
my.data$type<-ifelse(grepl("DNA",my.data$sample),"DNApositive",
                     ifelse(grepl("PCR",my.data$sample),"PCRpositive",
                            ifelse(grepl("Negative",my.data$sample),"negative","Moth sample")))
  
### subset the data to only keep the numbers of reads at each stage
my.data.subs<-subset(my.data, select=c("sample","total","trimmed.total","trimmed.pe","merged","type"))

### melt the data into long format
my.data.subs.melt<-melt(my.data.subs, id.vars=c("sample","type"))

### make the ggplot object + add the jittered dots + add the boxplots and colour them by type and make them a bit transparent
p<-ggplot(my.data.subs.melt, aes(x=type, y=value), fill=type) + geom_jitter() + geom_boxplot(aes(fill=type), alpha=0.5)
### tell ggplot that you want your data subsetted by a second variable
p + facet_wrap(~variable, scales="free_y")

### save the graph to an svg plot
ggsave(filename="read_summary.svg")

### calculating the proportional dropoff with each processing step
my.data$total.to.trimmed<-my.data$trimmed.total/my.data$total*100
my.data$trimmed.to.pe<-my.data$trimmed.pe/my.data$total*100

### subset the data to only keep the proportions of reads at each stage
my.propdata.subs<-subset(my.data, select=c("sample","total.to.trimmed","trimmed.to.pe","type"))
my.propdata.subs$total<-rep(100, length(my.propdata.subs))
my.propdata.subs<-my.propdata.subs[,c(1,5,2,3,4)]

### melt the data into long format
my.propdata.subs.melt<-melt(my.propdata.subs, id.vars=c("sample","type"))

qplot(variable,value, data=my.propdata.subs.melt, geom="boxplot",
      fill=variable,
      xlab="Level of processing",
      ylab="Proportion of total reads remaining")

### save the graph to an svg plot
ggsave(filename="proportion_summary.svg")

########################################################################################
################### summarisation of clustering parameters #############################
########################################################################################

### read in the data
my.clust<-read.csv("combined_reads_stats.csv", stringsAsFactors=FALSE)

### read in the intermediate data
###my.clust<-read.csv("combined_reads_stats_filter.csv", stringsAsFactors=FALSE)

### set a grouping variable for the data type
my.clust$type<-ifelse(grepl("DNA",my.clust$sample),"DNApositive",
                      ifelse(grepl("PCR",my.clust$sample),"PCRpositive",
                             ifelse(grepl("Negative",my.clust$sample),"negative","Moth sample")))

### subset the data appropriately

my.clust.subs<-subset(my.clust, my.clust$type=="Moth sample",
                      select=c("cluster_thres", "clusters_min_cov", "cluster_above_thres"))

#my.clust.subs<-subset(my.clust,
#                      grepl(paste("-", c("positive", "Negative"), sep="", collapse= "|"), my.clust$sample),
#                      select=c("cluster_thres", "clusters_min_cov", "cluster_above_thres"))

#test<-c("DNApositive","PCRpositive","Negative","Win")
#grepl(paste("-", c("positive", "Negative"), sep="", collapse= "|"), test)

#newdata <- test[paste(c("positive", "Negative"), sep="", collapse= "|")]

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

########################################################################################
################### kernel plots of clusters retained  #############
########################################################################################

### subset the data appropriately
my.clust2<-subset(my.clust, select=c("cluster_thres",
                                      "clusters_min_cov",
                                     "cluster_above_thres"))

### make the factors into factor format
my.clust2$cluster_thres<-as.factor(my.clust2$cluster_thres)
my.clust2$clusters_min_cov<-as.factor(my.clust2$clusters_min_cov)

### give the facets nice names by correcxtly naming the  factor levels
for(i in 1:length(levels(my.clust2$cluster_thres))){
levels(my.clust2$cluster_thres)<-paste(as.numeric(my.clust2$cluster_thres[1]*100))
}

### create a nice colour scale using colourRampPalette and RColorBrewer
jet.colors2 <-
  colorRampPalette(brewer.pal(length(unique(my.clust2$clusters_min_cov)), "Set3"))

############## plotting order highest to lowest threshold - hashed out to only select preferred order
### make the ggplot object
hl<-ggplot(my.clust2, aes(x=cluster_above_thres))
### add the kernal density plots
hl + geom_density(aes(fill=factor(clusters_min_cov)), alpha=0.5) +
  ### colour with our colour scale
  scale_fill_manual(name="Minimum number of\nsequences in cluster",
                    values = jet.colors2(length(unique(my.clust2$clusters_min_cov)))) +
  ### separate by clustering similarity threshold
  facet_wrap(~cluster_thres, scales="free", ncol=2) +
  ### set nice axis labels
  labs(x = "Clusters retained", y = "Kernel density")

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

########################################################################################
################### summarisation of positive read contamination position  #############
########################################################################################

### read in the data
my.contamination<-read.csv("Positive_contamination.csv")

### subset the data to include only the bits we need
my.contamination.subs<-subset(my.contamination, select=c("WELL","Total"))

### We will be logging the colour scale later so I have to add 1 to all my zero values
#my.contamination.subs$Total<-ifelse(my.contamination.subs$Total==0, my.contamination.subs$Total+1, my.contamination.subs$Total)
colnames(my.contamination.subs)<-c("WELL","Total")

### use dplyr to generate row and column numbers from the well references
my.contamination.mut<- mutate(my.contamination.subs,
                        Row=as.numeric(match(toupper(substr(WELL, 1, 1)), LETTERS)),
                        Column=as.numeric(substr(WELL, 2, 5)))

### plot this onto a plate with a colopur scale for positive sequences using the instructions
### at http://www.r-bloggers.com/plotting-microtiter-plate-maps/

### set up the ggplot, tell it that you want points for the wells and that these should be overlaid with points for read intensity based on values in "Total"
ggplot(data=my.contamination.mut, aes(x=Column, y=Row, z=Total)) +
  geom_point(size=25) +
  geom_point(aes(colour=Total), size=20) +
  
### Colour the points above how you like and set the position and style of the legend - unfortunately this can't be set out nicely in Rstudio
  scale_color_gradient("Positive DNA read intensity", low="#e5f6ff", high="red", na.value="black", trans="log", guide = guide_colourbar(barwidth=40, barheight=2, label=FALSE, ticks=FALSE, direction="horizontal", title.position="top", title.hjust=0.5)) +
  
  ### set up the axis labels and spacing of points in the style of a PCR plate
  coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
  scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
  scale_x_continuous(breaks=seq(1, 12)) +
  
  ### Add a title
  labs(title="Spatial contamination of positive DNA") +
  
  ### Add some last theme touches to make the diagram look like a plate (from the package ggplot2bdc)
  ### and boost the size of the labels a bit
  theme_bdc_microtiter() +
  theme(legend.position = "bottom", text = element_text(size=20))

### save the graph to an svg plot
ggsave(filename="contamination_distribution.svg")

########################################################################################
################### stacked bar charts of well composition   ###########################
########################################################################################

### read in the data
my.reads<-read.csv(file="metaBEAT.tsv", sep="\t", stringsAsFactors=FALSE, header=TRUE)
#colnames(my.reads)

### read in the sample by plate data
my.plates<-read.csv(file="/Volumes/JK BACKUPS/Illumina Data/OPM_MiSeq1/Samples_and_MIDS_corrected.txt",
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
occurance=5

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
write.csv(my.assignments, file="assignments_out.csv")

### Pull in the assigned colours as a .csv
Taxa.col<-read.csv("assignments_in.csv")

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
my.reads.trans$total<-rowSums(my.reads.trans[c(2:15)])
### calculate the percentage of reads in each well that are OPM
my.reads.trans$Percent.Thau<-(my.reads.trans$Thaumetopoea_processionea/my.reads.trans$total)*100
### work out a 20% of reads inclusion cutoff
#my.reads.trans$cutoff<-my.reads.trans$total*0.2

### replace all the species that are less than 20% of the total reads with zero
#my.reads.trans[,2:16][my.reads.trans[,2:16] < my.reads.trans$cutoff] <- 0

### subset the data frame to drop all coloumns containing only zeros
###my.reads.trans.drop<-cbind(my.reads.trans[,c(1, 34:39)], subset(my.reads.trans[,c(2:33)], select = colSums(my.reads.trans[,c(2:33)]) !=0))

### remove the total and variables
my.reads.trans.drop<-my.reads.trans[,c(1:19,21)]

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
hit.hist<-data.frame(OTUs = rowSums(my.reads.trans.samps.only[c(2:4,8,9)] != 0), Type=my.reads.trans.samps.only$type)


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

########################################################################################
# plotting a kernel density of the number of retained species per well after filtering #
########################################################################################



ggsave_golden("OTU_count_frequency.svg")

########################################################################################
################### make a plot of % composition by nest  ######################################
########################################################################################
svg(file="well_composition_reads_by_nest.svg", width=50, height=30)
### set up the ggplot
x<-ggplot(data=my.reads.melt, aes(x=Sample, y=Reads, fill=Species)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(position="fill", stat="identity") +
  ### wrap the plot by plate  
  facet_wrap(~Nest, scales="free_x") +
  ### give it a percentage scale
  scale_y_continuous(labels = percent_format()) +
  ### set the colours
  #  scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.reads.melt$Species)))) +
  scale_fill_manual(name="Species", values = colours) +
  ### add a sensible y axis label
  labs(y = "% of reads per well") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(2)),
        strip.text= element_text(size = rel(2)))

########################################################################################
################### make a plot of absolute read counts by nest ###############################
########################################################################################

### set up the ggplot
y<-ggplot(data=my.reads.melt, aes(x=Sample, y=Reads, fill=Species)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(stat="identity") +
  ### wrap the plot by plate  
  facet_wrap(~Nest, scales="free_x") +
  ### set the colours
  #  scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.reads.melt$Species)))) +
  scale_fill_manual(name="Species", values = colours) +
  ### add a sensible y axis label
  labs(y = "Absolute numbers of reads per well") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(2)),
        strip.text= element_text(size = rel(2)))

### arrange the plots into a grid - this is just a visualisation, this does not end up in the svg
grid.arrange(x, y, ncol=2)

### arrange the plots into a grid and close the svg dev
grid.arrange(x, y, ncol=2)
dev.off()

########################################################################################
################### stacked bar charts of well composition using clusters  ###########################
########################################################################################

### read in the data
my.clusters<-read.csv(file="inter_test_clusters.tsv", sep="\t", stringsAsFactors=FALSE, header=TRUE)
#colnames(my.clusters)

### read in the sample by plate data
my.plates<-read.csv(file="/Volumes/JK BACKUPS/Illumina Data/OPM_MiSeq1/Samples_and_MIDS_corrected.txt",
                    sep="\t", stringsAsFactors=FALSE, header=TRUE)
### trim the plate data to the necessary columns
my.plates<-my.plates[,1:3]

### transpose the cluster data
my.clusters.trans<-recast(my.clusters, variable~OTU_ID)
colnames(my.clusters.trans)[1]<-"sample"
my.clusters.trans$sample<-as.character(my.clusters.trans$sample)

### use match to add the plate data to the read data
my.clusters.trans$plate<-my.plates$plate[match(my.clusters.trans$sample,my.plates$sample)]
my.clusters.trans$plate.numeric<-my.plates$plate.numeric[match(my.clusters.trans$sample,my.plates$sample)]

### total all the clusters
my.clusters.trans$total<-rowSums(my.clusters.trans[c(2:18)])
### calculate the percentage of clusters in each well that are OPM
my.clusters.trans$Percent.Thau<-(my.clusters.trans$Thaumetopoea_processionea/my.clusters.trans$total)*100
### remove the total variable
my.clusters.trans<-my.clusters.trans[c(1:20,22)]

### make sample and plate factors for faceting and ordering
my.clusters.trans$sample<-as.factor(my.clusters.trans$sample)
my.clusters.trans$plate<-as.factor(my.clusters.trans$plate)

### melt the data into long format
my.clusters.melt<-melt(my.clusters.trans, id.vars=c("sample","plate", "plate.numeric","Percent.Thau"))
colnames(my.clusters.melt)<-c("Sample","Plate","Plate.numeric","Percent.Thau","Species","clusters")

### reorder my.clusters.melt2 by both percentage Thau and plate
my.clusters.melt<-my.clusters.melt[order(-my.clusters.melt$Percent.Thau,my.clusters.melt$Plate.numeric),]

### create an ID variable for the order
my.clusters.melt$level.order<-seq(1, length(my.clusters.melt$Percent.Thau),1)

### order the species for plotting the legend
#my.clusters.melt$Species <- factor(my.clusters.melt$Species,
#                      levels=c("Thaumetopoea.processionea",
#                               "Triops.cancriformis",
#                               "Comaster.audax",
#                               "Astatotilapia.calliptera",
#                               "Bilateria",
#                               "Pancrustacea",
#                               "nohit"))

### order the samples from highest OPM % to lowest and in increasing plate number
my.clusters.melt$Sample <- reorder(my.clusters.melt$Sample, my.clusters.melt$level.order)
### order the plates from lowest to highest
my.clusters.melt$Plate <- reorder(my.clusters.melt$Plate, my.clusters.melt$Plate.numeric)

### create a nice colour scale using colourRampPalette and RColorBrewer
#jet.colors3 <-
#  colorRampPalette(brewer.pal(length(unique(my.clusters.melt$Species)), "Spectral"))

### for specifying specific colours - you need to know how many levels you have
colours <- c("#a6cee3", "#1f78b4", "#b2df8a", "#FF0066", "#33a02c", "#fb9a99", 
             "#e31a1c", "#669999", "#fdbf6f", "#009933", "#ff7f00", "#cab2d6",
             "#FFFF00", "#6a3d9a", "#ffff99", "#b15928", "#000000")

########################################################################################
################### make a plot of % composition  ######################################
########################################################################################
svg(file="well_composition_clusters.svg", width=50, height=30)
### set up the ggplot
w<-ggplot(data=my.clusters.melt, aes(x=Sample, y=clusters, fill=Species)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(position="fill", stat="identity") +
  ### wrap the plot by plate  
  facet_wrap(~Plate, scales="free_x") +
  ### give it a percentage scale
  scale_y_continuous(labels = percent_format()) +
  ### set the colours
  #  scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.clusters.melt$Species)))) +
  scale_fill_manual(name="Species", values = colours) +
  ### add a sensible y axis label
  labs(y = "% of clusters per well", x="Individual Samples") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(2)),
        strip.text= element_text(size = rel(2)))

########################################################################################
################### make a plot of absolute read counts  ###############################
########################################################################################

### set up the ggplot
z<-ggplot(data=my.clusters.melt, aes(x=Sample, y=clusters, fill=Species)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(stat="identity") +
  ### wrap the plot by plate  
  facet_wrap(~Plate, scales="free_x") +
  ### set the colours
  #  scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.clusters.melt$Species)))) +
  scale_fill_manual(name="Species", values = colours) +
  ### add a sensibleaxis labels
  labs(y = "Absolute numbers of clusters per well", x="Individual Samples") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(2)),
        strip.text= element_text(size = rel(2)))

### arrange the plots into a grid and close the svg dev
grid.arrange(w, z, ncol=2)
dev.off()


########################################################################################
# work out percentage parasitism with carcelia #
########################################################################################

### divide the number of well containing Carcelia by all well that are not +ve or -ve to get percentage
percent.cacelia<-(colSums(my.reads.trans[2] > 0)/919)*100

### I may want to set a specific Carcelia cutoff at any value of carcelia lower than the highest +ve or -ve Carcelia read count should be deleted

### count the columns greater than zero and write to a new data frame
hit.hist2<-data.frame(OTUs = rowSums(my.reads.trans[c(2:17)] != 0), Type=my.reads.trans$type)

table(hit.hist)
table(hit.hist2)

