library(tidyverse)

bs_Gwa <- read.table(gzfile("budset.lrt0.gz"),
                     header = T)
bs_maf <- read.table(gzfile("budset.mafs.gz"),
                     header = T)

bs_Gwa <- bs_Gwa %>%
  filter(LRTscore > -998)
  

summary(bs_Gwa)

hist(bs_Gwa$beta)
# minor allele -> earlier budset (neg. direction)
# look at the subset of alleles with bio effect to look at which pops have sig. minor alleles?

# bedtools intersect

ggplot(bs_Gwa, aes(beta, LRTscore)) +
  geom_point()
# pchisq log.p = T lower.tail = F?

bs_Gwa$p <- -10*(pchisq(bs_Gwa$LRTscore, df = 1, log.p = T, lower.tail = F))
head(p)

bs_Gwa$position2 <- seq(1,nrow(bs_Gwa))

head(bs_Gwa)

ggplot(bs_Gwa, aes(position2, p)) +
  geom_point()

0.05/506042
# 70.05

quantile(bs_Gwa$p,0.9999)


# find outliers for the GWASs
bs_outliers <- bs_Gwa %>%
  filter(p > 184.8012)

write.table(bs_outliers, "bs_outliers.txt", 
            row.names = F)

# find intersection with annotations
# venn diagram
# temp and precipitation

# UTR group? 
###################################
# height
library(tidyverse)

ht_Gwa <- read.table(gzfile("height.lrt0.gz"),
                     header = T)
ht_maf <- read.table(gzfile("height.mafs.gz"),
                     header = T)

ht_Gwa <- ht_Gwa %>%
  filter(LRTscore > -998)


summary(ht_Gwa)

hist(ht_Gwa$beta)
# minor allele -> positive

# bedtools intersect

ggplot(ht_Gwa, aes(beta, LRTscore)) +
  geom_point()
# pchisq log.p = T lower.tail = F?

ht_Gwa$p <- -10*(pchisq(ht_Gwa$LRTscore, df = 1, log.p = T, lower.tail = F))

ht_Gwa$position2 <- seq(1,nrow(ht_Gwa))

head(ht_Gwa)

ggplot(ht_Gwa, aes(position2, p)) +
  geom_point()

0.05/506042
# 70.05

quantile(ht_Gwa$p,0.9999)


# find outliers for the GWASs
ht_outliers <- ht_Gwa %>%
  filter(p > 142.6555)

write.table(ht_outliers, "ht_outliers.txt", 
            row.names = F)


