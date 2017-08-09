#------------------------------- 
# Candidates-Comparison
# Author: Gitanshu Munjal
# Affiliation: Brummer Lab, UC Davis 
# Date: 10/06/2016
# Requirements: should already have run Ztest.R
#------------------------------- 




#----------------------Initialize

# #Clear workspace
# dev.off()
# rm(list=ls())                 



#Set working directory
setwd(paste("../",
            populationtitle,
            "/data/", sep = ""))



#Read candidate frequency tables prduced using F- tests

#PCadapt
can.pcadapt <- read.table(paste("../results/4-candidates-", 
                           populationtitle, ".txt", sep = ""), header = T, sep = "\t")



#Effective pop size estimates
pop.size <- read.table(paste("../results/5-EffectivePopSize-", 
                             populationtitle, ".txt", sep = ""), header = T, sep = "\t")


#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")


#Color scheme for plots
plot.colors <- c(alpha("grey", alpha = 0.5),
                 viridis(9, alpha = 0.5)[c(8)], 
                 viridis(9, alpha = 0.5)[c(3,2)])






#------------------ PCadapt


#Index by minor base allele
id           <- can.pcadapt$id
allele.index <- ifelse(can.pcadapt$base.1 > 0.5, "alt", "ref")
base.index   <- ifelse(allele.index=="alt", 1-can.pcadapt$base.1, can.pcadapt$base.1)
high.index   <- ifelse(allele.index=="alt", 1-can.pcadapt$high.1, can.pcadapt$high.1)
low.index    <- ifelse(allele.index=="alt", 1-can.pcadapt$low.1, can.pcadapt$low.1)
bh.ne        <- pop.size$BaseHighEstimate[5]

#Bin by minor base allele
bin          <- rep(cut(base.index, breaks = seq(0,0.5,0.05), right = TRUE),3)



frequency    <- round(c(base.index, high.index, low.index),2)
cycle        <- c(rep("3-base", nrow(can.pcadapt)), 
                  rep("2-less", nrow(can.pcadapt)), 
                  rep("1-more", nrow(can.pcadapt))) 



#Boxplot data
boxplot.data <- data.frame(id,
                           frequency, 
                           cycle, 
                           bin=factor(as.numeric(bin)))
  




#----- Calculate drift threshold for each marker

#Empty vector to store threhold for each marker
bh.limits <- c()

#loop over all markers
for (i in 1:length(base.index)){
  
  #Parameters for simulations 
  initial.freq <- base.index[i]
  N = round(bh.ne * 4,0)
  t = 4
  j = initial.freq*N
  t0 = 1
  nb_time = 4
  N_sample = 384
  plot <- FALSE
  sims <- 1000
  sleep <- 0
  
  #Run simulation
  sim.results <- WF_trajectory(N, t, j, t0, nb_time, N_sample, plot, sims, sleep)
  
  #max expected change due to drift  
  max.change.driftsim <- 
    quantile(abs(sim.results-initial.freq), probs = 0.99)
  
  obs.change.bh <- ifelse(abs(high.index[i] - base.index[i]) > max.change.driftsim, "s","ns")
  obs.change.bl <- ifelse(abs(low.index[i] - base.index[i]) > max.change.driftsim, "s","ns")
  
  bh.limits <- data.frame(rbind(bh.limits, cbind(obs.change.bh, obs.change.bl)))
  
}



#Define colors 
#If more than upper limit on drift, label it by population
results <- c(rep("3-base", nrow(can.pcadapt)), 
             as.character(bh.limits$obs.change.bh), 
             as.character(bh.limits$obs.change.bl))

boxplot.data$color <- ifelse(results=="s", boxplot.data$cycle, "0-drift")
boxplot.data$color <- ifelse(results=="3-base", boxplot.data$cycle, boxplot.data$color)


#----- Make drift plot
bin.vals <- levels(bin)


bh.plot <- ggplot() + 
  geom_jitter(data = boxplot.data, 
              aes(x = bin, y = frequency, color = color), 
              size = 2.5, height = 0, width = 0.2) +
  scale_color_manual(values = plot.colors) +
  scale_x_discrete(breaks=c(levels(boxplot.data$bin)), 
                   labels=bin.vals[as.numeric(levels(boxplot.data$bin))]) +
  #scale_shape_manual(values = c(22,24,21)) +
  theme_classic.bold + theme(legend.position="none") + 
  labs(x = paste(populationtitle, "-O bin"), y = "Frequency", 
       title = populationtitle)



#Write passed candidates to disk
can.pcadapt.passed <- boxplot.data[boxplot.data$color!="3-base", ]
can.pcadapt.passed <- unique(can.pcadapt.passed$id[can.pcadapt.passed$color!="0-drift"])
can.pcadapt.passed <- can.pcadapt[can.pcadapt$id %in% can.pcadapt.passed,]
write.table(can.pcadapt.passed, 
            paste("../results/6-pca-driftpassed-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, row.names = FALSE, sep = "\t")






#View both plots
drift.plot <- grid.arrange(bh.plot, ncol = 1)



#Write plot to disk
pdf(file = paste("../results/6-driftplot-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
grid.arrange(drift.plot)
dev.off()


pdf(file = paste("../results/6-directionplot-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
par(mfrow=c(1,3))
hist(can.pcadapt$high.1, 
     breaks = seq(0,1,0.1), ylim = c(0,70), 
     main = "H", xlab="Allele Frequency",
     col=viridis(8, alpha = 0.6)[6])
hist(can.pcadapt$base.1, 
     breaks = seq(0,1,0.1), ylim = c(0,70), 
     main = "O", xlab="Allele Frequency",
     col=viridis(8, alpha = 0.6)[8])
hist(can.pcadapt$low.1, 
     breaks = seq(0,1,0.1), ylim = c(0,70), 
     main = "L", xlab="Allele Frequency",
     col=viridis(8, alpha = 0.6)[1])
dev.off()

