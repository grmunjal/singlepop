#------------------------------- 
# AFS
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
  paste("../results/sub-bulk.freqs-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")

#Mapped frequencies
mapped.freqs <- read.table(
  paste("../results/sub-mapped.freqs-", populationtitle, ".txt", sep = ""), 
  header = T, sep = "\t")



#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")






#----------------------Allele frequency distributions

#-----------Using all markers

#Pull frequencies for each pop into a dataframe
afs.data = data.frame(rbind(
  cbind(freq = all.freqs$high.sub1.1, Population = paste(populationtitle, "-H-1", sep = "")),
  cbind(all.freqs$high.sub2.1, paste(populationtitle, "-H-2", sep = "")),
  cbind(all.freqs$high.sub3.1, paste(populationtitle, "-H-3", sep = "")),
  cbind(all.freqs$high.sub4.1, paste(populationtitle, "-H-4", sep = "")),
  
  cbind(all.freqs$base.sub1.1, paste(populationtitle, "-C0-1", sep = "")),
  cbind(all.freqs$base.sub2.1, paste(populationtitle, "-C0-2", sep = "")),
  cbind(all.freqs$base.sub3.1, paste(populationtitle, "-C0-3", sep = "")),
  cbind(all.freqs$base.sub4.1, paste(populationtitle, "-C0-4", sep = "")),
  
  cbind(all.freqs$low.sub1.1, paste(populationtitle, "-L-1", sep = "")),
  cbind(all.freqs$low.sub2.1, paste(populationtitle, "-L-2", sep = "")),
  cbind(all.freqs$low.sub3.1, paste(populationtitle, "-L-3", sep = "")),
  cbind(all.freqs$low.sub4.1, paste(populationtitle, "-L-4", sep = ""))))



#Make sure data frame has appropriate structure and use maf
afs.data$freq <- 1 - as.numeric(as.character(afs.data$freq))


#Plot AFS
#Purple is high, Green is low, Blue is base (colors set in the rtools script)

afs.plot <- ggplot(afs.data, aes(freq, color = Population)) +
  
#   geom_freqpoly(data = afs.data[afs.data$Cycle=="H3",],
#                 binwidth = 0.001, size = 1.35, colour = alpha("black", 0.4)) +
#   
#   geom_freqpoly(data = afs.data[afs.data$Cycle=="C0",],
#                 binwidth = 0.001, size = 1.35, colour = alpha("black", 0.4)) +
#   
#   geom_freqpoly(data = afs.data[afs.data$Cycle=="L3",],
#                 binwidth = 0.001, size = 1.35, colour = alpha("black", 0.4)) +
  
  geom_freqpoly(binwidth = 0.0125, size = 0.65, position = "identity",na.rm = TRUE, pad = TRUE) + 
  
  scale_color_manual(values = c(viridis(3, alpha = 0.4)[c(2,2,2,2,3,3,3,3)],
                     viridis(3, alpha = 0.4)[c(1,1,1,1)])) +
  
  theme_classic.adjust +
  
  theme(legend.title = element_text(size=16),
        legend.text=element_text(size=14),
        legend.position = c(0.85,0.85)) +
  labs(x = "Allele Frequency", y = "Count") +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  theme(axis.title.x = element_text(size = rel(1.5))) +
  theme(legend.title = element_text(size=16))



#View plot
print(afs.plot)



#Write plot to disk
pdf(file = paste("../results/2-afsplot-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
print(afs.plot)
dev.off()



