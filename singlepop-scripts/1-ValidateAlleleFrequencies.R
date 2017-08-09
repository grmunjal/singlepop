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

rownames(all.M) <- c("O-1", "O-2","O-3", "O-4",
                     "H-1", "H-2","H-3", "H-4",
                     "L-1", "L-2","L-3", "L-4")


#Make G-matrix and plot
all.G.matrix <- make.G.matrix(all.M)
prettyHeatmap(all.G.matrix)


#Do a PCA and plot
pca.alldata <- marker.pca(all.M)
#biplot(pca.alldata)



#Write plot to disk
pdf(file = paste("../results/1-gpca-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
prettyHeatmap(all.G.matrix)
prettyBiplot(all.M)
dev.off()


#-----------Using mapped markers
#Matrix where every column is a marker and every row is a sample (data are frequency)
mapped.M <- t(as.matrix(subset(mapped.freqs, select = 
                              -c(id, chr, pos, cluster, cadl.scaf, cadl.scaf.pos,base.1,high.1,low.1))))


#Make G-matrix and plot
mapped.G.matrix <- make.G.matrix(mapped.M)
heatmap(mapped.G.matrix,revC = TRUE)


#Do a PCA and plot
pca.mappeddata <- marker.pca(mapped.M)
#biplot(pca.mappeddata)
