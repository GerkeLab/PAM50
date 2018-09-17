# .libPaths(c("~/Documents/Rpackagesv2",.libPaths()))
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(here)

###################################################################################################
# step 1: import data
###################################################################################################

# import pam50 results from parker algo
pam50 <- read_delim("data/prostate_medCentering/tcgaPAM50_pam50scores.txt", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

# all clinical
clinical <- read_delim("data/data_bcr_clinical_data_patient.txt", 
                       "\t", escape_double = FALSE, comment = "#", 
                       trim_ws = TRUE)

# tcga with high quality rna samples 
clinical_pub <- read_delim("data/data_clinical.txt", 
                           "\t", escape_double = FALSE, comment = "#", 
                           trim_ws = TRUE)

###################################################################################################
# Step 2: data cleaning 
###################################################################################################

# shorten barcode to patient ID
pam50$patient_id <- substr(pam50$X1,1,12)

# reassign Calls to only normal, basal, luminal A/B
pam50$Call2 <- colnames(pam50[,c(2,4:6)])[apply(pam50[,c(2,4:6)],1,which.max)]

# merge pam50 calls with clinical info 
dat <- merge(clinical,pam50,by.x="PATIENT_ID",by.y="patient_id",all.x=TRUE)
dat$high_quality <- ifelse(dat$PATIENT_ID %in% clinical_pub$PATIENT_ID,1,0)

# create gleason values
dat$gleason <- ifelse(dat$GLEASON_SCORE==7,
                      paste0(dat$GLEASON_PATTERN_PRIMARY,"+",dat$GLEASON_PATTERN_SECONDARY),
                      dat$GLEASON_SCORE)
dat$gleason <- ifelse(dat$gleason %in% c("9","10"),"9-10",dat$gleason)

# clean psa and categorize 
dat$psa <- as.numeric(dat$PSA_MOST_RECENT_RESULTS)
dat$psa_grp <- ifelse(dat$psa<4,"0-4",
                      ifelse(dat$psa<10,"4-10",
                             ifelse(dat$psa<20,"10-20",">20")))


###################################################################################################
# step 3: summary statistics (feel free to skip)
###################################################################################################

summary(dat[dat$high_quality==1,]$AGE)

table(dat[dat$high_quality==1,]$BIOCHEMICAL_RECURRENCE_INDICATOR)

summary(as.numeric(dat[dat$high_quality==1,]$DAYS_TO_LAST_FOLLOWUP))

xtabs(~dat[dat$high_quality==1,]$gleason)

xtabs(~dat[dat$high_quality==1,]$PATH_T_STAGE)

xtabs(~dat[dat$high_quality==1,]$psa_grp)

xtabs(~dat[dat$high_quality==1,]$Call)

xtabs(~dat[dat$high_quality==1,]$Call2)

xtabs(~dat[dat$high_quality==1,]$gleason + dat[dat$high_quality==1,]$Call)

xtabs(~dat[dat$high_quality==1,]$gleason + dat[dat$high_quality==1,]$Call2)

###################################################################################################
# step 4: plot calls by confidence 
###################################################################################################

# PRIOR TO DROPPING NORMAL/HER2

ggplot(dat[dat$high_quality==1,],aes(x=Call,y=Confidence,fill=Call)) + 
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0,dodge.width=0.75),
              aes(fill=Call,col=Call),alpha=0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  guides(fill=FALSE,colour=FALSE)

# AFTER DROPPING HER2

ggplot(dat[dat$high_quality==1,],aes(x=Call2,y=Confidence,fill=Call2)) + 
  scale_fill_manual(values=rep("white",5)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5,jitter.height = 0,dodge.width=0.75),
              aes(fill=Call2,col=Call2),alpha=0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  guides(fill=FALSE,colour=FALSE)

###################################################################################################
# step 5: create heatmap
###################################################################################################
# load expression data 
exp <- read_delim("data/prostate_medCentering/tcgaPAM50.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
exp <- as.data.frame(exp)
colnames(exp)[2:ncol(exp)] <- substr(colnames(exp)[2:ncol(exp)],1,12)

row.names(exp) <- exp$gene

exp <- exp[,colnames(exp) %in% clinical_pub$PATIENT_ID]

name_order <- c(dat[dat$high_quality==1 & dat$Call2=="Basal",]$PATIENT_ID,
                dat[dat$high_quality==1 & dat$Call2=="LumA",]$PATIENT_ID,
                dat[dat$high_quality==1 & dat$Call2=="LumB",]$PATIENT_ID,
                dat[dat$high_quality==1 & dat$Call2=="Normal",]$PATIENT_ID)

exp <- exp[,name_order]

pam50calls <- dat[match(colnames(exp), dat$PATIENT_ID),]$Call2

dfs <- dat[match(colnames(exp), dat$PATIENT_ID),]$DFS_STATUS

glea <- as.character(dat[match(colnames(exp), dat$PATIENT_ID),]$gleason)

anno1 <- HeatmapAnnotation(df=as.data.frame(cbind(PAM50 = pam50calls,
                                                  Status=dfs,
                                                  Gleason=glea)),
                           col = list(PAM50 = c("Basal"="hotpink",
                                                "LumA"="cyan",
                                                "LumB"="orange",
                                                "Normal"="purple"))
)

Heatmap(exp,cluster_columns = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 7), top_annotation = anno1, show_heatmap_legend = FALSE)