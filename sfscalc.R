SFS <- scan("DG_folded_allsites.sfs")
SFS

sumSFS <- sum(SFS)
sumSFS

pctPoly = 100*(1-(SFS[1]/sumSFS))
pctPoly

plotSFS <- SFS[-c(1,length(SFS))]
barplot(plotSFS)
# appear to have lost the rare alleles because of drift
# expect Tajima's D to be above 0 in a bottlenecked population

div <- read.table("DG_folded_allsites.thetas.idx.pestPG")
colnames(div) <- c("window", "chrname","wincenter","tW","tP","tF",
                   "tH","tL","tajD","fulif","fuliD","fayH",
                   "zengsE","numSites")
head(div)

div$tWpersite <- div$tW/div$numSites
div$tPpersite <- div$tP/div$numSites
str(div)


par(mfrow=c(2,2))
barplot(plotSFS)
hist(div$tWpersite, col = "gray", xlab = "Theta-W", main = "")
hist(div$tPpersite, col = "gray", xlab = "Theta-Pi", main = "")
hist(div$tajD, col = "gray", xlab = "Tajima's D", main = "")

summary(div)
