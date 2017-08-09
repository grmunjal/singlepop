#------------------------------- 
# Truncatula genes in candidate regions
# Author: Gitanshu Munjal
# Affiliation: Brummer Lab, UC Davis 
# Date: 7/29/2016
# Requirements: should already have run Fst scan and have gff file
#------------------------------- 




#----------------------Initialize

# #Clear workspace
# dev.off()
# rm(list=ls())                 



#Set working directory
setwd(paste("../",
            populationtitle,
            "/data/", sep = ""))



#Read GFF file
genes <- read.table("./genesv2.txt",
                    header=F, sep="\t", quote = "")



#Read candidate frequency tables prduced using F- tests

#PCadapt
can.pcadapt <- read.table(paste("../results/4-candidates-", 
                            populationtitle, ".txt", sep = ""), header = T, sep = "\t")


#High/Base
can.drift.passed <- read.table(paste("../results/6-pca-driftpassed-", 
                                     populationtitle, ".txt", sep = ""), header = T, sep = "\t")



#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")
system("rm ../results/*genes*.txt")






#----------------------Truncatula based annotations for candidates

#Parameters
regionsize <- 1000*20 #(look for genes +/- this many bp from candidate)



#------------PCadapt candidates

#Remove unmapped markers
can.pcadapt.mapped <- subset(can.pcadapt, can.pcadapt$chr != "chr9")

#Upper bound of significant region = SNP located regionsize above significant SNP 
can.pcadapt.mapped$upperbound = can.pcadapt.mapped$pos - regionsize

#Lower bound of significant region = SNP located regionsize below significant SNP 
can.pcadapt.mapped$lowerbound = can.pcadapt.mapped$pos + regionsize



#----Get genes in candidate regions from gff file
#Loop over all candidates
for(i in 1:nrow(can.pcadapt.mapped)) {
  
  #Subset GFF file by chr for candidate i
  selected.genes <- subset(genes, genes[,1] == as.character(can.pcadapt.mapped$chr[i])) 
  
  #Subset above subset for positions above upperbound
  selected.genes <- subset(selected.genes, selected.genes[,4] > can.pcadapt.mapped$upperbound[i])   
  #Subset above subset for positions below lowerbound
  selected.genes <- subset(selected.genes, selected.genes[,5] < can.pcadapt.mapped$lowerbound[i])   
  
  #Write significant genes in region to table (table will have lots of duplicate entries)
  write.table(selected.genes, 
              paste("../results/7-can.pcadapt-genes-",
                    populationtitle, ".txt", sep = ""), 
              quote=F, sep="\t", append=T, col.names=F , row.names = F)
}



#Read above table, remove duplicate entries, and order by chromosome and position
selected.genes <- read.table(paste("../results/7-can.pcadapt-genes-",
                                   populationtitle, ".txt", sep = ""), 
                             header = F, sep="\t" ,quote = "")
selected.genes <- selected.genes[!duplicated(selected.genes), ]        
selected.genes <- selected.genes[order(selected.genes[,1], selected.genes[,4]), ]



#Write to disk
write.table(selected.genes, paste("../results/7-can.pcadapt-genes-",
                                  populationtitle, ".txt", sep = ""), 
            quote=F, sep="\t", col.names=F, row.names=F, append = FALSE)





