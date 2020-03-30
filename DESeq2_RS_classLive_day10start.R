## Set your working directory
# setwd("~/github/2020_Ecological_Genomics")

## Import the libraries that we're likely to need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

## Import the counts matrix
countsTable <- read.table("RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)

## Import the samples description table - links each sample to factors of the experimental design.
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)

############ Try with only Day 10 data

# grep("10", names(countsTableRound), value = TRUE)
# day10countstable <- subset(countsTableRound, grep("10", names(countsTableRound), value = TRUE)) #doesn't work has to be logical

day10countstable <- countsTableRound %>% select(contains("10"))
dim(day10countstable)

conds10<- subset(conds, day=="10")
dim(conds10)
head(conds10)

## Let's see how many reads we have from each sample:
colSums(day10countstable)
mean(colSums(day10countstable))
barplot(colSums(day10countstable), las=3, cex.names=0.5,names.arg = substring(colnames(day10countstable),1,13))
abline(h=mean(colSums(day10countstable)), col="blue", lwd =2)

# What's the average number of counts per gene
rowSums(day10countstable)
mean(rowSums(day10countstable))
median(rowSums(day10countstable))
# wow! This shows dispersion across genes - differences in magnitude of expression

# What's the average number of counts per gene per sample
apply(day10countstable,2,mean)

## Create a DESeq object and define the experimental design here with the tilde

dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)
dim(dds)
# [1] 66408    30

# Filter out genes with few reads 

dds <- dds[rowSums(counts(dds)) > 30]
dim(dds)
# 24300    30
## Run the DESeq model to test for differential gene expression: 
# 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 
# 3) run negative binomial glm
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# Running the model: design = ~ climate + treatment + climate:treatment

# [1] "Intercept"            "climate_HD_vs_CW"     "treatment_D_vs_C"    
# [4] "treatment_H_vs_C"     "climateHD.treatmentD" "climateHD.treatmentH"


# Order and list and summarize results from specific contrasts
# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
resClimHDtrtD <- results(dds, name="climateHD.treatmentD", alpha = 0.05)
resClimHDtrtD <- resClimHDtrtD[order(resClimHDtrtD$padj),]
head(resClimHDtrtD)

summary(resClimHDtrtD)
# out of 24300 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 118, 0.49%
# low counts [2]     : 0, 0%

### Same with climateHD.treatmentH
resClimHDtrtH <- results(dds, name="climateHD.treatmentH", alpha = 0.05)
resClimHDtrtH <- resClimHDtrtH[order(resClimHDtrtH$padj),]
head(resClimHDtrtH)

summary(resClimHDtrtH)
# out of 24300 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 118, 0.49%
# low counts [2]     : 0, 0%
##

# climate_HD_vs_CW
resClimHDvCW <- results(dds, name="climate_HD_vs_CW", alpha = 0.05)
resClimHDvCW <- resClimHDvCW[order(resClimHDvCW$padj),]
head(resClimHDvCW)

summary(resClimHDvCW)
# out of 24300 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.0041%
# outliers [1]       : 118, 0.49%
# low counts [2]     : 0, 0%

# treatment_D_vs_C
resTrtDvC <- results(dds, name="treatment_D_vs_C", alpha = 0.05)
resTrtDvC <- resTrtDvC[order(resTrtDvC$padj),]
head(resTrtDvC)

summary(resTrtDvC)
# out of 24300 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 257, 1.1%
# LFC < 0 (down)     : 327, 1.3%
# outliers [1]       : 118, 0.49%
# low counts [2]     : 10730, 44%

# treatment_H_vs_C
resTrtHvC <- results(dds, name="treatment_H_vs_C", alpha = 0.05)
resTrtHvC <- resTrtHvC[order(resTrtHvC$padj),]
head(resTrtHvC)

summary(resTrtHvC)
# out of 24300 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0041%
# LFC < 0 (down)     : 2, 0.0082%
# outliers [1]       : 118, 0.49%
# low counts [2]     : 0, 0%


res_interClimTreat <- results(dds, name="climateHD.treatmentD", alpha=0.05)
res_interClimTreat <- res_interClimTreat[order(res_interClimTreat$padj),]
head(res_interClimTreat)


summary(res_interClimTreat)
# out of 23887 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 678, 2.8%
# LFC < 0 (down)     : 424, 1.8%
# outliers [1]       : 61, 0.26%
# low counts [2]     : 7367, 31%

##############################
# contrasts:
# resClimHDtrtD
# resClimHDtrtH

# resClimHDvCW
# resTrtDvC
# resTrtHvC

# MA Plot
# Data visualization -----------------------------------------------------------
# MA plot
ggmaplot(resClimHDvCW, 
         main = expression("Source climate hot/dry" %->% "Source climate cool/wet"),
         fdr = 0.01, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         legend = "top", 
         top = 0,
         # genenames = as.vector(res_treatCD$gene),
         # select.top.method = "padj",
         # font.label = c("bold", 11), 
         # label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  theme(legend.text = element_text(size = 12))

# PCA
vsd <- vst(dds, blind = FALSE, nsub = 10000)
data <- plotPCA(vsd, 
                ntop = 10000,
                intgroup = c("climate", "treatment"), 
                returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

# reordering the factors
data$treatment <- factor(data$treatment, 
                         levels = c("C", "H", "D"), 
                         labels = c("C", "H", "D"))
data$day <- factor(data$day, 
                   levels = c("0", "5", "10"), 
                   labels = c("0","5","10"))

ggplot(data, aes(PC1, PC2, color = climate, shape=treatment)) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()


