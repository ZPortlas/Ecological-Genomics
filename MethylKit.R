library(methylKit)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# first, we want to read in the raw methylation calls with methylkit

# set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)
dir <- "C:/Users/Zoe/Desktop/Ecological-Genomics"
# read in the sample ids
samples <- read.table("sample_id.txt", header = F)
# now point to coverage files
files <- file.path(dir, samples$V1)
all(file.exists(files))
# convert to list
file.list <- as.list(files)
# get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))
# use methRead to read in the coverage files
myobj <- methRead(location= file.list, sample.id =   nmlist,
                  assembly = "atonsa", # this is just a string. no actual database
                  dbtype = "tabix", context = "CpG",
                  resolution = "base", mincov = 20,
                  treatment = 
                    c(0,0,0,0,
                      1,1,1,1,
                      2,2,2,2,
                      3,3,3,3,
                      4,4,4,4),
                  pipeline = "bismarkCoverage",
                  dbdir = dir)

######
# visualize coverage and filter
######

# We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats(myobj[[1]], plot = T)

# and can plot all of our samples at once to compare.

# filter samples by depth with filterByCoverage()
filtered.myobj <- filterByCoverage(myobj,
                                   lo.count = 20,
                                   lo.perc = NULL,
                                   hi.count = NULL,
                                   hi.perc = 97.5,
                                   dbdir = dir)

######
# merge samples
######

#Note! This takes a while and we're skipping it

# use unite() to merge all the samples. We will require sites to be present in each sample or else will drop it

meth <- unite()


############ 
meth <- methylKit:::readMethylBaseDB(
  dbpath = "methylBase_united.txt.bgz",
  dbtype = "tabix",
  sample.id =   unlist(nmlist),
  assembly = "atonsa", # this is just a string. no actual database
  context = "CpG",
  resolution = "base",
  treatment = c(0,0,0,0,
                1,1,1,1,
                2,2,2,2,
                3,3,3,3,
                4,4,4,4),
  destrand = FALSE)

# percMethylation() calculates the percent methylation for each site and sample
pm <- percMethylation(meth)

#plot methylation histograms
ggplot(gather(as.data.frame(pm)), aes(value)) + 
  geom_histogram(bins = 10, color="black", fill="grey") + 
  facet_wrap(~key, ncol = 4)

# calculate and plot mean methylation
sp.means <- colMeans(pm)

p.df <- data.frame(sample=names(sp.means),
              group = substr(names(sp.means), 1,6),
                   methylation = sp.means)
ggplot(p.df, aes(x=group, y=methylation, color=group)) + 
  stat_summary(color="black") + geom_jitter(width=0.1, size=3)

# is the fact methylation is higher in AA_F00 real or is it because the mapping rate is so much higher?

# sample clustering
clusterSamples(meth, dist = "correlation", method = "ward.D", plot = T)
# things aren't clustering by treatment

# PCA
PCASamples(meth, screeplot = F)

# subset with reorganize()

meth_sub <- reorganize(meth,sample.ids =
                         c("AA_F00_1","AA_F00_2",
                           "AA_F00_3", "AA_F00_4",
                           "HH_F25_1","HH_F25_2",
                           "HH_F25_3","HH_F25_4"),
                       treatment = c(0,0,0,0,1,1,1,1),
                       save.db=FALSE)

# calculate differential methylation

myDiff <- calculateDiffMeth(meth_sub,
                            overdispersion = "MN",
                            mc.cores = 1,
                            suffix = "AA_HH",
                            adjust = "qvalue",
                            test = "Chisq")

# get all differentially methylated bases
myDiff <- getMethylDiff(myDiff,qvalue = 0.05,
                         difference = 10)
# we can visualize the changes in methylation frequencies quickly.
hist(getData(myDiff)$meth.diff)
# get hyper methylated bases
hyper <- getMethylDiff(myDiff, difference = 10, 
                       qvalue = 0.05, type = "hyper")

# get hypo methylated bases
hypo <- getMethylDiff(myDiff, difference = 10, 
                      qvalue = 0.05, type = "hypo")


#heatmaps first

# get percent methylation matrix
pm <- percMethylation(meth_sub)
# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]



# add snp, chr, start, stop

din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, 
                      getData(myDiff)$start, 
                      sep=":"), 
                din, pm.sig)

colnames(df.out) <- c("snp", colnames(din), 
                      colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))
# make a dataframe with snp id's, methylation, etc.

# add snp, chr, start, stop


####
# heatmap
####

pheatmap(pm.sig, show_rownames = FALSE)

# we can also normalize 
ctrmean <- rowMeans(pm.sig[,1:4])
h.norm <- (pm.sig-ctrmean)
pheatmap(h.norm, show_rownames = FALSE)


#####
#let's look at methylation of specific snps
####

# convert data frame to long form
df.plot <- df.out[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")

df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)

# Looking at snp LS049205.1:248
# if you choose a different snp, you can create different plots.

df.plot %>% filter(snp=="LS049205.1:248") %>% 
  ggplot(., aes(x=group, y=methylation, 
                color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")

## write bed file for intersection with genome annotation
setwd(dir)
write.table(file = "diffmeth.bed",
            data.frame(chr = df.out$chr, 
                       start = df.out$start, 
                       end = df.out$end),
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\t")
##################################################
# write up exp
# AA_F25 and AH_F25
#################################################
# subset with reorganize()
meth_sub <- reorganize(meth, sample.ids =
                         c("AA_F25_1","AA_F25_2",
                           "AA_F25_3", "AA_F25_4",
                           "AH_F25_1","AH_F25_2",
                           "AH_F25_3","AH_F25_4"),
                       treatment = c(0,0,0,0,1,1,1,1),
                       save.db=FALSE) 

# diff methylation
myDiff <- calculateDiffMeth(meth_sub,
                            overdispersion = "MN",
                            mc.cores = 1,
                            suffix = "AA_AH",
                            adjust = "qvalue", 
                            test = "Chisq") 

myDiff <- getMethylDiff(myDiff,
                        qvalue = 0.05, 
                        difference = 10) 

summary(getData(myDiff)$meth.diff)
#  Min.     1st Qu. Median  Mean   3rd Qu.  Max. 
# -17.194 -14.451 -11.866  -7.447 -10.417  17.477 
hist(getData(myDiff)$meth.diff) # methylation differences 

# hist below 0 are hypomethylated and those above 0 are hypermethylated
# we see a lower methylation % in AH_F25 samples compared to AA_F25 samples

# get hyper methylated bases
hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper")

# get hypo methylated bases
hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo")

hist(getData(hyper)$meth.diff)
hist(getData(hypo)$meth.diff)

#heatmaps
pm <- percMethylation(meth_sub)

# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]

# add snp, chr, start, stop
din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, 
                      getData(myDiff)$start, sep=":"), 
                din, pm.sig)
colnames(df.out) <- c("snp", 
                      colnames(din), 
                      colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))

####
# heatmap
####

my_heatmap <- pheatmap(pm.sig,show_rownames = F) # red is higher methylation

# we can also normalize 
ctrmean <- rowMeans(pm.sig[,1:4])
h.norm <- (pm.sig-ctrmean)
my_heatmap <- pheatmap(h.norm,show_rownames = FALSE)

#####
#let's look at methylation of specific snps
####

# convert data frame to long form
head(df.out,3)
df.plot <- df.out[,c(1,5:12)] %>% 
  pivot_longer(-snp, values_to = "methylation")

df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)
df.plot[9:16,]
# looking at snp LS051659.1:1214

df.plot %>% filter(snp=="LS051659.1:1214") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black") +
  ggtitle("Methylation at SNP LS051659.1:1214")

# write bed file for genome annotation
write.table(file = "AA_AH_diffmeth.bed",
            data.frame(chr = df.out$chr, 
                       start = df.out$start, 
                       end = df.out$end),
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

#################################################
# AA_F25 and HH_F25
################################################
# subset with reorganize()
meth_sub <- reorganize(meth, sample.ids =
                         c("AA_F25_1","AA_F25_2",
                           "AA_F25_3", "AA_F25_4",
                           "HH_F25_1","HH_F25_2",
                           "HH_F25_3","HH_F25_4"),
                       treatment = c(0,0,0,0,1,1,1,1),
                       save.db=FALSE) 

# diff methylation
myDiff <- calculateDiffMeth(meth_sub,
                            overdispersion = "MN",
                            mc.cores = 1,
                            suffix = "AA_HH",
                            adjust = "qvalue", 
                            test = "Chisq") 

myDiff <- getMethylDiff(myDiff,
                        qvalue = 0.05, 
                        difference = 10) 

summary(getData(myDiff)$meth.diff)
# Min.     1st Qu.  Median  Mean   3rd Qu.  Max. 
# -32.956 -16.331  11.160   1.032  15.490  32.097 
hist(getData(myDiff)$meth.diff) # methylation differences 

# hist below 0 are hypomethylated and those above 0 are hypermethylated
# we see a slightly higher methylation % in HH_F25 samples compared to AA_F25 samples

# get hyper methylated bases
hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper")

# get hypo methylated bases
hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo")

hist(getData(hyper)$meth.diff)
hist(getData(hypo)$meth.diff)

#heatmaps
pm <- percMethylation(meth_sub)

# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]

# add snp, chr, start, stop
din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, 
                      getData(myDiff)$start, sep=":"), 
                din, pm.sig)
colnames(df.out) <- c("snp", 
                      colnames(din), 
                      colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))

####
# heatmap
####

my_heatmap <- pheatmap(pm.sig,show_rownames = F) # red is higher methylation

# we can also normalize 
ctrmean <- rowMeans(pm.sig[,1:4])
h.norm <- (pm.sig-ctrmean)
my_heatmap <- pheatmap(h.norm,show_rownames = FALSE)

#####
#let's look at methylation of specific snps
####

# convert data frame to long form
head(df.out,3)
df.plot <- df.out[,c(1,5:12)] %>% 
  pivot_longer(-snp, values_to = "methylation")

df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)

# looking at snp LS042748.1:7742

df.plot %>% filter(snp=="LS042748.1:7742") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black") +
  ggtitle("Methylation at SNP LS042748.1:7742")

# write bed file for genome annotation
write.table(file = "AA_HH_diffmeth.bed",
            data.frame(chr = df.out$chr, 
                       start = df.out$start, 
                       end = df.out$end),
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\t")
