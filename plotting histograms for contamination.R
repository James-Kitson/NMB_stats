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