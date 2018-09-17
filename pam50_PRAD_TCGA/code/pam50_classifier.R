##############################################################################
# making PAM50 calls  
##############################################################################
# setting up enviroment, load libraries

# .libPaths(c("~/Documents/Rpackagesv2",.libPaths()))
options(stringsAsFactors=FALSE)

library(readr)
library(here)
library(tidyverse)

source("code/medianCenteringFUNCTION.R")

###################################################################################################
# step 1: load and clean data
###################################################################################################

genes <- read_delim("data/tcgaPAM50.txt", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
genesT <- as.data.frame(t(genes))
sample_names <- colnames(genes)[2:ncol(genes)]
gene_names <- genes$Hugo_Symbol

genesT <- as.data.frame(cbind(sample_names,genesT[-1,]))
colnames(genesT)[2:ncol(genesT)] <- gene_names

for(i in 2:ncol(genesT)){
  genesT[,i] <- as.numeric(genesT[,i])
}



###################################################################################################
# step 2: perform median gene centering 
###################################################################################################

# need tidyverse library to run median_centering fucntion

testing <- median_centering(input.data=genesT,
                            subsetN=100,
                            subsetSize=100,
                            percentERpos=NA,
                            idVar="sample_names")

###################################################################################################
# step 3: save median centered values for parker algo
###################################################################################################

write.table(testing, file="data/prostate_medCentering/tcgaPAM50.txt", sep="\t",
            row.names = FALSE, quote = FALSE)

###################################################################################################
# step 4: run parker algo
###################################################################################################

library(ctc)
library(heatmap.plus)

paramDir<- "fromUNC/PAM50/bioclassifier_R" # the location of unchanging files such as the function library and main program
inputDir<- "data/prostate_medCentering"  # the location of the data matrix, and where output will be located

inputFile<- "tcgaPAM50.txt" # the input data matrix as a tab delimited text file
short<-"tcgaPAM50" # short name that will be used for output files

calibrationParameters<- NA 	#the column of the "mediansPerDataset.txt" file to use for calibration; 
#NA will force centering within the test set & -1 will not do any 
#adjustment (when adjustment performed by used)

hasClinical<-FALSE 	#may include tumor size as second row, with 'T' as the gene name, 
#and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
#set this variable to FALSE if tumor size is not available

collapseMethod<-"mean" # can be mean or iqr (probe with max iqr is selected)
# typically, mean is preferred for long oligo and
# iqr is preferred for short oligo platforms


####
# run the assignment algorithm
####

source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))
