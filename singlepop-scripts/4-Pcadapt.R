#------------------------------- 
# Validate Allele Frequencies
# Author: Gitanshu Munjal
# Affiliation: Brummer Lab, UC Davis 
# Date: 7/29/2016
# Requirements: should already have run sub-CalculateAlleleFrequencies.R
#------------------------------- 




#----------------------Initialize

# #Clear workspace
# dev.off()
# rm(list=ls())                 



#Set working directory
setwd(paste("../",
            populationtitle,
            "/data/", sep = ""))



#Read frequency tables produced using CalculateAlleleFrequencies.R

#All frequencies 
all.freqs <- read.table(
  paste("../results/sub-bulk.freqs-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")

#Mapped frequencies
mapped.freqs <- read.table(
  paste("../results/sub-mapped.freqs-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")



#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")






#----------------------Analysis


#-----------Using all markers
#Matrix where every column is a marker and every row is a sample (data are frequency)
all.M <- t(as.matrix(subset(all.freqs, select = 
                              -c(id, chr, pos, cluster, cadl.scaf, cadl.scaf.pos,base.1,high.1,low.1))))



#pcadapt
library(pcadapt)

#Read frequency data
pcadapt.data <- read.pcadapt(all.M,
                         type="pool",
                         ploidy = 4,
                         pop.sizes = rep(250,nrow(all.M)))



#Perform PCA...retain 20 PCs
pcadapt.results <- pcadapt(pcadapt.data,K=20)



#Adjust depending on population
poplist.names <- c(rep("O", bulk.census[[populationtitle]]$base*250/24),
                   rep("H", bulk.census[[populationtitle]]$high*250/24),
                   rep("L", bulk.census[[populationtitle]]$low*250/24))



#View summary plots
plot(pcadapt.results, option="screeplot")
plot(pcadapt.results, option="scores", pop = poplist.names)



#Write plots to disk
pdf(file = paste("../results/4-screescoreplot-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
plot(pcadapt.results, option="screeplot")
plot(pcadapt.results, option="scores", pop = poplist.names)
dev.off()






#-----------Calculate test statistic (D^2)

#Pull scores for PC1
pc.scores <- pcadapt.results$scores[,1]



#Numerator for Z-score calculation assuming 1 PC (equals 1 with no missing data)
z.num <- sum((pc.scores^2))



#Calculate Z-scores for every marker 
#(page 6 on http://biorxiv.org/content/biorxiv/early/2016/05/30/056135.full.pdf)
g.mat <- t(pcadapt.data)
x.mat <- as.matrix(pc.scores)



#Regression coefficients
betas <- t(g.mat) %*% x.mat



#Residual matrix
resid.mat <- g.mat - (x.mat %*% t(x.mat) %*% g.mat)



#Calculate Z-score for marker
pc.zscores <- c()
for(marker in 1:nrow(pcadapt.data)){
  pcr.resid.var <- sum((resid.mat[,marker])^2)
  
  pc.zscore <- betas[marker]*sqrt(z.num/pcr.resid.var)
  
  pc.zscores <- c(pc.zscores, pc.zscore)
}
all.freqs$zscores <- pc.zscores



library(MASS)
#Estimate of variance
sigma.estimate <- cov.rob(all.freqs$zscores)



#Calculate test statistic
all.freqs$dstat <- (all.freqs$zscores - mean(pc.zscores))^2 / sigma.estimate$cov



#Calculate inflation factor 
gif.chi <- rchisq(length(pc.zscores),1)
gif     <- median(all.freqs$dstat)/median(gif.chi)



#Correct test statistic
all.freqs$dstat.corrected <- all.freqs$dstat/gif



#Write results to disk
write.table(all.freqs, 
            paste("../results/4-pcadaptresults-", 
                  populationtitle, ".txt", sep = ""), 
            quote = FALSE, row.names = FALSE, sep = "\t")




#Threshold
sig.level = 0.999
thresh.dst  <- quantile(all.freqs$dstat.corrected, probs=sig.level, na.rm = T)



#Find candidate markers at threshold and write to disk
candidates <- subset(all.freqs, all.freqs$dstat.corrected >= thresh.dst)
dim(candidates)
write.table(candidates, 
            paste("../results/4-candidates-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, row.names = FALSE, sep = "\t")






#---------------Manhattan plot

#Data frame for plot
man.data.pca <- data.frame(chr=as.character(all.freqs$chr),
                           pos=as.numeric(all.freqs$pos),
                           var=sqrt(all.freqs$dstat.corrected),
                           special=NA)


man.pca <- manhattan.plot.special(man.data.pca,
                                   paste(populationtitle),
                                  sqrt(thresh.dst))
#View plot
man.pca + ylim(0,30) + labs(y = "d")


#Write plot to disk
pdf(file = paste("../results/4-pcadaptplot-",
                 populationtitle,
                 ".pdf", sep = ""), width = 8, height = 6)
man.pca + ylim(0,30) + labs(y = "d")
dev.off()


