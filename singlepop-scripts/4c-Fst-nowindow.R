#------------------------------- 
# Fst Scan
# Author: Gitanshu Munjal
# Affiliation: Brummer Lab, UC Davis 
# Date: 7/29/2016
# Requirements: should already have run CalculateAlleleFrequencies.R
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
  paste("../results/bulk.freqs-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")

#Mapped frequencies
mapped.freqs <- read.table(
  paste("../results/mapped.freqs-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")



#Read count tables produced using CalculateAlleleFrequencies.R

#All counts
all.counts <- read.table(
  paste("../results/bulk.counts-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")

#Mapped counts
mapped.counts <- read.table(
  paste("../results/mapped.counts-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")

# #Neutral POD
# neutral.counts <- read.table("../results/Gpool.POD-NEU-1",  
#   header = F, sep = " ")


#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")






#----------------------Pairwise Bulk Segregant Analysis----------------------#


#----------- Data frame to pull from -----------#
bsa.freqs  <- all.freqs
bsa.counts <- all.counts



#Order both data frames 
bsa.freqs  <- bsa.freqs[order(bsa.freqs$id),]
bsa.counts <- bsa.counts[order(bsa.counts$id),]






#----------- Define outlier/significance level -----------#
sig.level         <- 0.99
bsa.freqs$special <- NA






#----------- Calculate Fst -----------#


#-----------High / Base

#Parameters for fst
x1 = bsa.counts$high.1
m1 = bsa.counts$high.1 + bsa.counts$high.2

x2 = bsa.counts$base.1
m2 = bsa.counts$base.1 + bsa.counts$base.2


#Calculate Fst
bsa.freqs$fbh <- vectorfst.tet(x1, x2, m1, m2)
hist(bsa.freqs$fbh, main = "Base / High Fst", xlab = "", breaks = seq(0,1,0.005))

thresh.bh <- quantile(bsa.freqs$fbh, probs=sig.level, na.rm=T)




#-----------High / Low

#Parameters for Fst
x1 = bsa.counts$low.1
m1 = bsa.counts$low.1  + bsa.counts$low.2

x2 = bsa.counts$high.1
m2 = bsa.counts$high.1 + bsa.counts$high.2


#Calculate Fst
bsa.freqs$fhl <- vectorfst.tet(x1, x2, m1, m2)
hist(bsa.freqs$fhl, main = "High / Low Fst", xlab = "", breaks = seq(0,1,0.005))





#-----------Low / Base

#Parameters for fst
x1 = bsa.counts$base.1
m1 = bsa.counts$base.1+bsa.counts$base.2

x2 = bsa.counts$low.1
m2 = bsa.counts$low.1+bsa.counts$low.2


#Calculate Fst
bsa.freqs$fbl <- vectorfst.tet(x1, x2, m1, m2)
hist(bsa.freqs$fbl, main = "Base / Low Fst", xlab = "", breaks = seq(0,1,0.005), col=viridis(1, alpha = 0.4))





#----Plots
pdf(file = paste("../results/4-fstplots-",
                 populationtitle,
                 ".pdf", sep = ""), width = 11, height = 8)

par(mfrow=c(1,2))
calibration.data <-
  data.frame(Fst=c(bsa.freqs$fhl,
                   bsa.freqs$fbh,
                   bsa.freqs$fbl),
             Population=rep(c("H-L", "O-H","O-L"),
                            each=nrow(bsa.freqs)))




#Plot 1
hist(bsa.freqs$fbh, freq = FALSE, breaks = seq(0,1,0.05), 
     ylim = c(0,30), col = viridis(8, alpha = 0.6)[2], 
     xlab=expression('F'["ST"]), main = populationtitle)
hist(bsa.freqs$fbl, freq = FALSE, breaks = seq(0,1,0.05),
     ylim = c(0,25), add = T, col = viridis(8, alpha = 0.6)[7], lty=3)
abline(v=thresh.bh)
legend(x = "right", y = "top", 
       legend = c(levels(calibration.data$Population)[c(2,3)]),
       fill = viridis(8, alpha=0.7)[c(2,7)], bty = "n", cex=1.2)



#Plot 2
hist(bsa.freqs$fbh, freq = FALSE, breaks = seq(0,1,0.05), 
     ylim = c(0,0.5), col = viridis(8, alpha = 0.6)[2], 
     xlab=expression('F'["ST"]), main = populationtitle)
hist(bsa.freqs$fbl, freq = FALSE, breaks = seq(0,1,0.05),
     ylim = c(0,25), add = T, col = viridis(8, alpha = 0.6)[7], lty=3)
abline(v=thresh.bh)
legend(x = "right", y = "top", 
       legend = c(levels(calibration.data$Population)[c(2,3)]),
       fill = viridis(8, alpha=0.7)[c(2,7)], bty = "n", cex=1.2)

dev.off()
