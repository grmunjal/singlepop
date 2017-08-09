#------------------------------- 
# Effective PopSize Analysis
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



#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")





#----------------------Ne from Nest (https://github.com/ThomasTaus/Nest)

#Load library
library(Nest)


#------------High-Base Analysis

#Parameters (see ?estimateNe help for details...)
p0       <- all.freqs$base.1          
pt       <- all.freqs$high.1 
cov0     <- all.counts$base.1 + all.counts$base.2
covt     <- all.counts$high.1 + all.counts$high.2
t        <- 3
Ncensus  <- 500
poolSize <- c(bulk.census[[populationtitle]]$base,
              bulk.census[[populationtitle]]$high)
ploidy   <-  4
method   <- c("P.planI", "JR.planI", "W.planI",
              "P.planII", "JR.planII", "W.planII")

estNe.bh <- estimateNe(p0=p0, pt=pt, cov0=cov0, covt=covt, 
                       t=t, method=method, 
                       Ncensus=Ncensus, poolSize=poolSize, ploidy = ploidy)



#------------Low-Base Analysis

#Parameters (see ?estimateNe help for details...)
p0       <- all.freqs$base.1          
pt       <- all.freqs$low.1 
cov0     <- all.counts$base.1 + all.counts$base.2
covt     <- all.counts$low.1 + all.counts$low.2 
t        <- 3
Ncensus  <- 500
poolSize <- c(bulk.census[[populationtitle]]$base,
              bulk.census[[populationtitle]]$low)
ploidy   <-  4
method   <- c("P.planI", "JR.planI", "W.planI",
              "P.planII", "JR.planII", "W.planII")

estNe.bl <- estimateNe(p0=p0, pt=pt, cov0=cov0, covt=covt, 
                       t=t, method=method, 
                       Ncensus=Ncensus, poolSize=poolSize, ploidy = ploidy)



#Tabulate results
ne.nest <- data.frame(BaseHighEstimate=estNe.bh, BaseLowEstimate=estNe.bl)
ne.nest


  
#Write results
write.table(ne.nest, 
            paste("../results/5-EffectivePopSize-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, row.names = TRUE, sep = "\t")
