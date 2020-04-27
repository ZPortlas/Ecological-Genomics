# install.packages('RcppCNPy')
# install.packages('label.switching')

library(RcppCNPy)
library(tidyverse)
library(ggfortify)

# reads results from admixure
edge <-npyLoad('edge_PC.admix.Q.npy')
edge <- as.data.frame(edge)
edge

bp = barplot(t(as.matrix(edge)), 
             space = c(0.2),
             col=rainbow(2),
             xlab="Individual #", 
             ylab="Ancestry",
             border=NA)

# reads estimated covariance matrix
cov <- as.matrix(read.table("edge_PC.cov"))
# compute eigenvalues and eigenvectors
cov.eigen <- eigen(cov)

plot(cov.eigen$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(cov.eigen$values)

# Perform PCA on covariance matrix
edge.pca <- prcomp(cov)
edge.pca
summary(edge.pca)

# plot unscaled PCs
pca.plot <- autoplot(edge.pca, data = edge)
pca.plot

# scaling
scaling <- edge.pca$sdev[1:2] * sqrt(nrow(edge))

pc1 <- rowSums(t(t(sweep(edge, 2 ,colMeans(edge))) *cov.eigen$vectors[,1] * -1) / scaling[1])
pc2 <- rowSums(t(t(sweep(edge, 2, colMeans(edge))) *cov.eigen$vectors[,2]) / scaling[2])

# collect PCs in data.frame and plot 
df <- data.frame(pc1, pc2)

ggplot(df, aes(x=pc1, y=pc2)) + 
  geom_point()

# save df without row or column names for GWAS
write.table(df, "edge_PC.txt", 
            row.names = F, col.names = F)
