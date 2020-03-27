library(tximportData)
library(tximport)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(readr)
library(wesanderson)
library(vsn)

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


## Let's see how many reads we have from each sample:
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)

# what's the average num of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))
# shows dispersion across genes, differences in magnitude of expression

apply(countsTableRound,2,mean)
## Create a DESeq object and define the experimental design here with the tilde
dds <- DESeqDataSetFromMatrix(
  countData = countsTableRound, 
  colData = conds,
  design = ~ climate + day + treatment)
dim(dds)
# [1] 66408    76

# Filter out genes with few reads
dds <-dds[rowSums(counts(dds)) > 76]
dim(dds)
# [1] 23887    76
# filtering to a sum of 76 read across all samples

## Run the DESeq model to test for differential gene expression: 
# 1) estimate size factors (per sample), 
# 2) estimate dispersion (per gene), 
# 3) run negative binomial glm
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

# running model: design ~ pop+day+treatment
# [1] "Intercept"            "pop_BRU_05_vs_ASC_06"
# [3] "pop_CAM_02_vs_ASC_06" "pop_ESC_01_vs_ASC_06"
# [5] "pop_JAY_02_vs_ASC_06" "pop_KAN_04_vs_ASC_06"
# [7] "pop_LOL_02_vs_ASC_06" "pop_MMF_13_vs_ASC_06"
# [9] "pop_NOR_02_vs_ASC_06" "pop_XBM_07_vs_ASC_06"
# [11] "day_10_vs_0"          "day_5_vs_0"          
# [13] "treatment_D_vs_C"     "treatment_H_vs_C" 

# running model: design ~ climate+day+treatment
# [1] "Intercept"        "climate_HD_vs_CW"
# [3] "day_10_vs_0"      "day_5_vs_0"      
# [5] "treatment_D_vs_C" "treatment_H_vs_C"

# Order and list and summarize results from specific contrasts
res <- results(dds, alpha = 0.05)
res <- res[order(res$padj),]
head(res)

summary(res)
# out of 23887 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 16, 0.067% (higher in hot, lower in control)
# LFC < 0 (down)     : 3, 0.013% (higher in control, lower in hot)
# outliers [1]       : 61, 0.26%
# low counts [2]     : 14300, 60%

res_treatCD <- results(dds, name="treatment_D_vs_C", alpha = 0.05)
res_treatCD <- res_treatCD[order(res_treatCD$padj),]
head(res_treatCD)

summary(res_treatCD)
# out of 23887 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 678, 2.8%
# LFC < 0 (down)     : 424, 1.8%
# outliers [1]       : 61, 0.26%
# low counts [2]     : 7367, 31%

# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast


##### Data visualization #####
# MA plot
plotMA(res_treatCD, ylim = c(-3,3))
# red is transcripts that are sig expressed
# range of expression is pretty low, median of ~24 makes sense

# PCA
vsd <- vst(dds,blind = F)

data <- plotPCA(vsd, intgroup = c("climate","treatment", "day"), returnData = T)
percentVar <- round(100 * attr(data, "percentVar"))


data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))

data$day <- factor(data$day, levels=c("0","5","10"), labels = c("0","5","10"))

ggplot(data, aes(PC1, PC2, color=day, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Counts of specific top gene! (important validatition that the normalization, model is working)

d <-plotCounts(dds, gene="MA_10426407g0030", intgroup = (c("treatment","climate")), returnData=TRUE)
d



p <-ggplot(d, aes(x=climate, y=count, shape=climate, colour = treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("CW","HD"))
p


d <-plotCounts(dds, gene="MA_10257300g0010", intgroup = (c("treatment","climate","day")), returnData=TRUE)
d



p <-ggplot(d, aes(x=treatment, y=count, color=day, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p
# Heatmap of top 20 genes sorted by pvalue
library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate")])
pheatmap(mat, annotation_col=df)
