#------------------------------- 
# Calculate allele frequencies
# Author: Gitanshu Munjal
# Affiliation: Brummer Lab, UC Davis 
# Date: 7/29/2016
#------------------------------- 




#----------------------Initialize

# #Clear workspace
# dev.off()
# rm(list=ls())                 



#Set working directory
setwd(paste("../",
            populationtitle,
            "/data/", sep = ""))



#Read master matrix from SNP-calling using GBS-SNP-CROP
print(getwd())
unfiltered.snps <- read.table("parsed/demultiplexed/SNPsCalled/SNPS_counts-R.txt",
                              header = FALSE,
                              row.names = 1,
                              sep = "\t")



#Number of SNPs and 2*Number of individuals
dim(unfiltered.snps)



#Source packages and functions
source("../../singlepop-scripts/singlepop-rtools.R")






#----------------------Make matrices of bulk counts by cycles


#Allele 1 counts
allele1.counts <- unfiltered.snps[, seq(1, ncol(unfiltered.snps), 2)]
ncol(allele1.counts)



#Allele 2 counts
allele2.counts <- unfiltered.snps[, seq(2, ncol(unfiltered.snps), 2)]
ncol(allele2.counts)



#Define bulks
bulk.descriptions[[populationtitle]]$base     -> base 
bulk.descriptions[[populationtitle]]$high     -> high 
bulk.descriptions[[populationtitle]]$low      -> low 




#Cycle specific bulk matrices
base.bulk <- cbind(base.1 = rowSums(allele1.counts[, base]),
                   base.2 = rowSums(allele2.counts[, base]))

high.bulk <- cbind(high.1 = rowSums(allele1.counts[, high]),
                   high.2 = rowSums(allele2.counts[, high]))

low.bulk <- cbind(low.1 = rowSums(allele1.counts[, low]),
                  low.2 = rowSums(allele2.counts[, low]))



#All bulks combined into a dataframe
bulk.counts <- data.frame(cbind(base.bulk, high.bulk, low.bulk))






#----------------------Depth Filters


#-----1-Filter lower limit on marker depth within each population

#Start a filter column in the count data frame and set to 0
bulk.counts$filter    <- 0


#Set filter limits to census size
bulk.census[[populationtitle]]$base     -> lolimit.base
bulk.census[[populationtitle]]$high     -> lolimit.high
bulk.census[[populationtitle]]$low      -> lolimit.low

# 96     -> lolimit.base 
# 96     -> lolimit.high 
# 96     -> lolimit.low 


#If, marker depth is less than lolimit 
#then, change filter column to non-zero
bulk.counts$filter <- ifelse(bulk.counts$base.1 +bulk.counts$base.2 < lolimit.base,
                             bulk.counts$filter + 1,
                             bulk.counts$filter) 

bulk.counts$filter <- ifelse(bulk.counts$high.1 +bulk.counts$high.2 < lolimit.high,
                             bulk.counts$filter + 1,
                             bulk.counts$filter)

bulk.counts$filter <- ifelse(bulk.counts$low.1 +bulk.counts$low.2 < lolimit.low,
                             bulk.counts$filter + 1,
                             bulk.counts$filter)



#Subset count table to rows where filter is still zero
bulk.counts <- subset(bulk.counts, bulk.counts$filter==0)



#Calculate depth of each marker and make plots
markerdepth <- rowSums(bulk.counts)

#View four plots at a time
par(mfrow=c(2,2))
hist(markerdepth, main = "Overall")
hist(rowSums(bulk.counts[,1:2]), main = "O")
hist(rowSums(bulk.counts[,3:4]), main = "H")
hist(rowSums(bulk.counts[,5:6]), main = "L")



#-----2-Filter upper limit on marker depth within each population

#Start a filter column in the count data frame and set to 0
bulk.counts$filter    <- 0


#Set filter limits 
mean(rowSums(base.bulk)) + 2*sd(rowSums(base.bulk))     -> uplimit.base 
mean(rowSums(high.bulk)) + 2*sd(rowSums(high.bulk))     -> uplimit.high 
mean(rowSums(low.bulk))  + 2*sd(rowSums(low.bulk))      -> uplimit.low 


#If, marker depth is greater than uplimit 
#then, change filter column to non-zero
bulk.counts$filter <- ifelse(bulk.counts$base.1 +bulk.counts$base.2 > uplimit.base,
                             bulk.counts$filter + 1,
                             bulk.counts$filter) 

bulk.counts$filter <- ifelse(bulk.counts$high.1 +bulk.counts$high.2 > uplimit.high,
                             bulk.counts$filter + 1,
                             bulk.counts$filter)

bulk.counts$filter <- ifelse(bulk.counts$low.1 +bulk.counts$low.2 > uplimit.low,
                             bulk.counts$filter + 1,
                             bulk.counts$filter)



#Subset count table to rows where filter is still zero
bulk.counts <- subset(bulk.counts, bulk.counts$filter==0)


#Calculate depth of each marker and make plots
markerdepth <- rowSums(bulk.counts)


#View four plots at once
par(mfrow=c(2,2))
hist(markerdepth, main = "Overall")
hist(rowSums(bulk.counts[,1:2]), main = "O")
hist(rowSums(bulk.counts[,3:4]), main = "H")
hist(rowSums(bulk.counts[,5:6]), main = "L")


#Write plots to disk
pdf(file = paste("../results/0-coverageplot-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
par(mfrow=c(2,2))
hist(markerdepth, main = "Overall")
hist(rowSums(bulk.counts[,1:2]), main = "O")
hist(rowSums(bulk.counts[,3:4]), main = "H")
hist(rowSums(bulk.counts[,5:6]), main = "L")
dev.off()


#----------------------End Filters

#Remove filter column from count table
bulk.counts <- subset(bulk.counts, select = -c(filter))






#----------------------Calculate frequencies per marker per sample

#Do this to get the basic structure for the table of frequencies
bulk.freqs <- bulk.counts 


#Calculate frequencies for allele1 (allele 2 is symmetric)
bulk.freqs$base.1 <- bulk.counts$base.1/(bulk.counts$base.1 + bulk.counts$base.2)
bulk.freqs$high.1 <- bulk.counts$high.1/(bulk.counts$high.1 + bulk.counts$high.2)
bulk.freqs$low.1  <- bulk.counts$low.1/(bulk.counts$low.1 + bulk.counts$low.2)


#Round to 3 decimal points and reduce table to allele 1 only
bulk.freqs <- round( bulk.freqs[, seq(1, ncol(bulk.freqs), 2)] ,3)






#----------------------Mapping information 


#-----Mapping information from truncatula

#Read sam file
mt.mapping.info <- read.delim(pipe("cut  -f1,2,3,4,5 parsed/demultiplexed/SNPsCalled/SNPsMt.bwa.sam"),
                              header = FALSE)


#Remove unmapped
mt.mapping.info <- subset(mt.mapping.info, 
                          !(mt.mapping.info$V3 %in% c("","*")))


#Keep primary alignments only
mt.mapping.info <- subset(mt.mapping.info, 
                          mt.mapping.info$V2 %in% c(0,16))


#Keep high scoring alignments only (thresh = -10*log10(p(wrong alignment)))
mt.mapping.info <- subset(mt.mapping.info, 
                          mt.mapping.info$V5 >= round(-10*log10(0.001),0))


#Reduce mapping information to identifier and coordinates only
rownames(mt.mapping.info) <- mt.mapping.info$V1
mt.mapping.info           <- mt.mapping.info[, 1:4]
mt.mapping.info           <- mt.mapping.info[, -2]
colnames(mt.mapping.info) <- c("cluster", "chr", "pos")





#-----Mapping information from blasting against CADL
#Read sam file
cadl.mapping.info <- read.delim(pipe("cut  -f1,2,3,4,5 parsed/demultiplexed/SNPsCalled/SNPsCADL.bwa.sam"),
                                header = FALSE)


#Remove unmapped
cadl.mapping.info <- subset(cadl.mapping.info, 
                            !(cadl.mapping.info$V3 %in% c("","*")))


#Keep primary alignments only
cadl.mapping.info <- subset(cadl.mapping.info, 
                            cadl.mapping.info$V2 %in% c(0,16))


#Keep high scoring alignments only (thresh = -10*log10(p(wrong alignment)))
cadl.mapping.info <- subset(cadl.mapping.info, 
                            cadl.mapping.info$V5 >= round(-10*log10(0.001),0))


#Reduce mapping information to identifier and coordinates only
rownames(cadl.mapping.info) <- cadl.mapping.info$V1
cadl.mapping.info           <- cadl.mapping.info[, 1:4]
cadl.mapping.info           <- cadl.mapping.info[, -2]
colnames(cadl.mapping.info) <- c("cluster", "chr", "pos")





#-----Mapping information from blasting against mock reference genome  
mr.mapping.info <- read.table("parsed/demultiplexed/SNPsCalled/SNP_clustermap.desc.txt",
                              header = FALSE,
                              sep = "\t")


#Reduce mapping information to identifier and coordinates only
mr.mapping.info           <- subset(mr.mapping.info, select = c(V1, V4))
colnames(mr.mapping.info) <- c("id", "cluster")





#-----Combine the three sources of mapping information (by identifier column)
snp.mapping.info <- merge(mt.mapping.info, mr.mapping.info, 
                          by = "cluster",
                          all.y = TRUE)

snp.mapping.info <- merge(snp.mapping.info, cadl.mapping.info, 
                          by = "cluster",
                          all.x = TRUE)

colnames(snp.mapping.info) <- 
  c("cluster", "chr", "pos", "id", "cadl.scaf", "cadl.scaf.pos")





#----------- "Chromosome 9" (Unmapped and scaffold mapped loci)-----------#

#Add 2 new levels to chromosomes
levels(snp.mapping.info$chr) <- c(levels(snp.mapping.info$chr), "chr9", "chr0")



#If mapped to alfalfa but not truncatula then chromosome 0
snp.mapping.info$chr[is.na(snp.mapping.info$chr) & !(is.na(snp.mapping.info$cadl.scaf))] <- 
  "chr0"



#If mapped to truncatula scaffolds then chromosome 9
snp.mapping.info$chr[(snp.mapping.info$chr %in% levels(snp.mapping.info$chr)[grep("scaf*", levels(snp.mapping.info$chr))])] <- "chr9"



snp.mapping.info$pos <- ifelse(is.na(snp.mapping.info$pos), 1, snp.mapping.info$pos)

snp.mapping.info$pos <- ifelse(snp.mapping.info$chr=="chr9", 1, snp.mapping.info$pos)





#----------------------Prepare count and frequency tables

#Add identifier column to bulk matrices
bulk.counts.prep <- na.omit(cbind(id=as.character(rownames(bulk.counts)),
                                  bulk.counts))
bulk.freqs.prep  <- na.omit(cbind(id=as.character(rownames(bulk.freqs)), 
                                                  bulk.freqs))


#Merge bulk matrices with mapping information
bulk.counts.prep <- merge(snp.mapping.info, 
                          bulk.counts.prep, by = "id", all.y = TRUE)
bulk.counts.prep <- subset(bulk.counts.prep, !(is.na(bulk.counts.prep$chr)))

bulk.freqs.prep  <- merge(snp.mapping.info, 
                          bulk.freqs.prep, by = "id", all.y = TRUE)
bulk.freqs.prep <- subset(bulk.freqs.prep, !(is.na(bulk.freqs.prep$chr)))


#Mapped frequencies
mapped.counts.prep <- subset(bulk.counts.prep, !(bulk.counts.prep$chr %in% c("chr9","chr0")))
mapped.freqs.prep  <- subset(bulk.freqs.prep, !(bulk.freqs.prep$chr %in% c("chr9","chr0")))






#----------------------Write count and frequency tables

#Bulk counts
write.table(bulk.counts.prep, 
            paste("../results/bulk.counts-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

#Bulk frequencies
write.table(bulk.freqs.prep, 
            paste("../results/bulk.freqs-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

#Mapped bulk counts
write.table(mapped.counts.prep, 
            paste("../results/mapped.counts-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

#Mapped bulk frequencies
write.table(mapped.freqs.prep,
            paste("../results/mapped.freqs-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)






#----------------------Make matrices of bulk counts by bio reps within cycles


#Allele 1 counts
allele1.counts <- unfiltered.snps[, seq(1, ncol(unfiltered.snps), 2)]
ncol(allele1.counts)



#Allele 2 counts
allele2.counts <- unfiltered.snps[, seq(2, ncol(unfiltered.snps), 2)]
ncol(allele2.counts)



#Define bio reps
bulk.bio.descriptions[[populationtitle]]$base1     -> base1
bulk.bio.descriptions[[populationtitle]]$base2     -> base2
bulk.bio.descriptions[[populationtitle]]$base3     -> base3
bulk.bio.descriptions[[populationtitle]]$base4     -> base4

bulk.bio.descriptions[[populationtitle]]$high1     -> high1
bulk.bio.descriptions[[populationtitle]]$high2     -> high2
bulk.bio.descriptions[[populationtitle]]$high3     -> high3
bulk.bio.descriptions[[populationtitle]]$high4     -> high4

bulk.bio.descriptions[[populationtitle]]$low1      -> low1
bulk.bio.descriptions[[populationtitle]]$low2      -> low2
bulk.bio.descriptions[[populationtitle]]$low3      -> low3
bulk.bio.descriptions[[populationtitle]]$low4      -> low4




#Start a data frame to pull bulked bio rep counts into
bulk.bio.counts <- as.data.frame(matrix(NA, nrow = nrow(allele1.counts), ncol = 1))
colnames(bulk.bio.counts) <- "dummycolumn"




#Base bio reps
bulk.bio.counts$base.sub1.1 = rowSums(allele1.counts[, base1])
bulk.bio.counts$base.sub1.2 = rowSums(allele2.counts[, base1])

bulk.bio.counts$base.sub2.1 = rowSums(allele1.counts[, base2])
bulk.bio.counts$base.sub2.2 = rowSums(allele2.counts[, base2])

bulk.bio.counts$base.sub3.1 = rowSums(allele1.counts[, base3])
bulk.bio.counts$base.sub3.2 = rowSums(allele2.counts[, base3])

bulk.bio.counts$base.sub4.1 = rowSums(allele1.counts[, base4])
bulk.bio.counts$base.sub4.2 = rowSums(allele2.counts[, base4])




#High bio reps
bulk.bio.counts$high.sub1.1 = rowSums(allele1.counts[, high1])
bulk.bio.counts$high.sub1.2 = rowSums(allele2.counts[, high1])

bulk.bio.counts$high.sub2.1 = rowSums(allele1.counts[, high2])
bulk.bio.counts$high.sub2.2 = rowSums(allele2.counts[, high2])

bulk.bio.counts$high.sub3.1 = rowSums(allele1.counts[, high3])
bulk.bio.counts$high.sub3.2 = rowSums(allele2.counts[, high3])

bulk.bio.counts$high.sub4.1 = rowSums(allele1.counts[, high4])
bulk.bio.counts$high.sub4.2 = rowSums(allele2.counts[, high4])




#Low bio reps
bulk.bio.counts$low.sub1.1 = rowSums(allele1.counts[, low1])
bulk.bio.counts$low.sub1.2 = rowSums(allele2.counts[, low1])

bulk.bio.counts$low.sub2.1 = rowSums(allele1.counts[, low2])
bulk.bio.counts$low.sub2.2 = rowSums(allele2.counts[, low2])

bulk.bio.counts$low.sub3.1 = rowSums(allele1.counts[, low3])
bulk.bio.counts$low.sub3.2 = rowSums(allele2.counts[, low3])

bulk.bio.counts$low.sub4.1 = rowSums(allele1.counts[, low4])
bulk.bio.counts$low.sub4.2 = rowSums(allele2.counts[, low4])


#Remove dummy column
bulk.bio.counts <- subset(bulk.bio.counts, select = -c(dummycolumn))

#Remove columns with no data
present.columns  <- as.numeric(which(round(colMeans(bulk.bio.counts),0)!=0))
bulk.bio.counts <- bulk.bio.counts[,present.columns]





#----------------------Calculate frequencies per marker per sample

#Do this to get the basic structure for the table of frequencies
bulk.bio.freqs <- bulk.bio.counts 


#Calculate frequencies for allele1 (allele 2 is symmetric)

#Base bio reps
bulk.bio.freqs$base.sub1.1 <- 
  bulk.bio.counts$base.sub1.1/(bulk.bio.counts$base.sub1.1 + bulk.bio.counts$base.sub1.2)
bulk.bio.freqs$base.sub2.1 <- 
  bulk.bio.counts$base.sub2.1/(bulk.bio.counts$base.sub2.1 + bulk.bio.counts$base.sub2.2)
bulk.bio.freqs$base.sub3.1 <- 
  bulk.bio.counts$base.sub3.1/(bulk.bio.counts$base.sub3.1 + bulk.bio.counts$base.sub3.2)
bulk.bio.freqs$base.sub4.1 <- 
  bulk.bio.counts$base.sub4.1/(bulk.bio.counts$base.sub4.1 + bulk.bio.counts$base.sub4.2)


#High bio reps
bulk.bio.freqs$high.sub1.1 <- 
  bulk.bio.counts$high.sub1.1/(bulk.bio.counts$high.sub1.1 + bulk.bio.counts$high.sub1.2)
bulk.bio.freqs$high.sub2.1 <- 
  bulk.bio.counts$high.sub2.1/(bulk.bio.counts$high.sub2.1 + bulk.bio.counts$high.sub2.2)
bulk.bio.freqs$high.sub3.1 <- 
  bulk.bio.counts$high.sub3.1/(bulk.bio.counts$high.sub3.1 + bulk.bio.counts$high.sub3.2)
bulk.bio.freqs$high.sub4.1 <- 
  bulk.bio.counts$high.sub4.1/(bulk.bio.counts$high.sub4.1 + bulk.bio.counts$high.sub4.2)


#Low bio reps
bulk.bio.freqs$low.sub1.1 <- 
  bulk.bio.counts$low.sub1.1/(bulk.bio.counts$low.sub1.1 + bulk.bio.counts$low.sub1.2)
bulk.bio.freqs$low.sub2.1 <- 
  bulk.bio.counts$low.sub2.1/(bulk.bio.counts$low.sub2.1 + bulk.bio.counts$low.sub2.2)
bulk.bio.freqs$low.sub3.1 <- 
  bulk.bio.counts$low.sub3.1/(bulk.bio.counts$low.sub3.1 + bulk.bio.counts$low.sub3.2)
bulk.bio.freqs$low.sub4.1 <- 
  bulk.bio.counts$low.sub4.1/(bulk.bio.counts$low.sub4.1 + bulk.bio.counts$low.sub4.2)


#Round to 3 decimal points and reduce table to allele 1 only
bulk.bio.freqs <- round( bulk.bio.freqs[, seq(1, ncol(bulk.bio.freqs), 2)] ,3)






#----------------------Filters

#Add id column to previous tables
bulk.bio.counts.prep <- 
  cbind(id = as.character(rownames(unfiltered.snps)), bulk.bio.counts)

bulk.bio.freqs.prep  <- 
  cbind(id = as.character(rownames(unfiltered.snps)), bulk.bio.freqs)



#Retain only those markers as retained by pop level bulks 

#Merge the pop and bio rep level tables and retain only filtered SNPs
sub.mapped.counts <- 
  merge(mapped.counts.prep, bulk.bio.counts.prep, by = "id", all.x = TRUE)

sub.mapped.freqs  <- 
  merge(mapped.freqs.prep, bulk.bio.freqs.prep, by = "id", all.x = TRUE)

sub.bulk.counts   <- 
  merge(bulk.counts.prep, bulk.bio.counts.prep, by = "id", all.x = TRUE)

sub.bulk.freqs    <- 
  merge(bulk.freqs.prep, bulk.bio.freqs.prep, by = "id", all.x = TRUE)



#----------------------Write count and frequency tables

#Bulk counts
write.table(sub.bulk.counts,
            paste("../results/sub-bulk.counts-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

#Bulk frequencies
write.table(sub.bulk.freqs,
            paste("../results/sub-bulk.freqs-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

#Mapped bulk counts
write.table(sub.mapped.counts,
            paste("../results/sub-mapped.counts-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

#Mapped bulk frequencies
write.table(sub.mapped.freqs,
            paste("../results/sub-mapped.freqs-", populationtitle, ".txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)
