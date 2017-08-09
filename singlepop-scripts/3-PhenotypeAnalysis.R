#-------------------------------
# Phenotype Analysis
# Author: Gitanshu Munjal
# Affiliation: Brummer Lab, UC Davis
# Date: 7/29/2016
# Requirements: pheno-1-R.txt file
#-------------------------------




#----------------------Initialize

# #Clear workspace
# dev.off()
# rm(list=ls())                 



#Set working directory
setwd(paste("../",
            populationtitle,
            "/data/", sep = ""))



#Read phenotypes file
phenotypes <-
  read.table("../../../phenotypes/pheno-1-R.txt", header = TRUE, sep = "\t")



#Source custom tools for processing
source("../../singlepop-scripts/singlepop-rtools.R")






#----------------------Davis Evaluation (2015)

#Make factors
phenotypes$background <- factor(phenotypes$background)
phenotypes$cycle      <- factor(phenotypes$cycle)
phenotypes$row        <- factor(phenotypes$row)
phenotypes$bc         <- factor(paste(phenotypes$background, phenotypes$cycle, sep = "-"))



#Standaradize height data and remove missing rows
phenotypes$height.std <- scale(phenotypes$height, 
                               center = TRUE, scale = TRUE)
phenotypes <- na.omit(phenotypes)



#Regress out effect of field rows and pull residuals
model1              <- lm(height.std ~ row, data = phenotypes)
phenotypes$adjusted <- model1$residuals






#----------------------Data munging to get population specific data 


#Pull population specific data from overall (eg: for CUF)
phenotypes.pop <-
  subset(phenotypes, phenotypes$background == populationtitle)



#Pull population-H3 (less dormancy data) from overall into a vector
lessdormancy.pheno <- subset(phenotypes.pop, phenotypes.pop$cycle == "H3")$adjusted
length(lessdormancy.pheno)


 
#Pull population-L3 (more dormancy data) from overall into a vector
moredormancy.pheno <- subset(phenotypes.pop, phenotypes.pop$cycle == "L3")$adjusted
length(moredormancy.pheno)



#Pull population-O (starting population data) from overall into a vector
startingpop.pheno <- subset(phenotypes.pop, phenotypes.pop$cycle == "C0")$adjusted
length(startingpop.pheno)



#Re-arrange these vectors into a dataframe
phenotypes.pop <- data.frame(rbind(
  cbind(pheno = moredormancy.pheno, Population = paste(populationtitle,"-L3", sep = "")),
  cbind(pheno = lessdormancy.pheno, Population = paste(populationtitle,"-H3", sep = "")),
  cbind(pheno = startingpop.pheno, Population = paste(populationtitle,"-C0", sep = ""))
))



#Convert measurements to numeric
phenotypes.pop$pheno <- as.numeric(as.character(phenotypes.pop$pheno))



#Get means and standard deviations
data.frame(
  Population = c("H3 / Less Dorm", "Base", "L3 / More Dorm"),
  MeanHeight = c(
    mean(lessdormancy.pheno),
    mean(startingpop.pheno),
    mean(moredormancy.pheno)
  ),
  sdHeight = c(
    sd(lessdormancy.pheno),
    sd(startingpop.pheno),
    sd(moredormancy.pheno)
  )
)






#----------------------Tests for phenotypic distributions


#Are the distributions normal?
print(shapiro.test(moredormancy.pheno))
print(shapiro.test(lessdormancy.pheno))
print(shapiro.test(startingpop.pheno))



#Are variances equal?
print(var.test(startingpop.pheno, lessdormancy.pheno))
print(var.test(startingpop.pheno, moredormancy.pheno))
print(var.test(moredormancy.pheno, lessdormancy.pheno))



#Are means equal? (adjust TRUE/FALSE depending on variances)
print(t.test(startingpop.pheno, lessdormancy.pheno, var.equal = TRUE))
print(t.test(startingpop.pheno, moredormancy.pheno, var.equal = TRUE))
print(t.test(moredormancy.pheno, lessdormancy.pheno, var.equal = TRUE))






#----------------------Phenotype Figure

#Make density plot
pheno.plot <- ggplot(phenotypes.pop, aes(pheno, fill = Population)) +
  geom_density(color = "black") +
  scale_fill_manual(values = viridis(3, alpha = 0.5)[c(2,3,1)]) +
  labs(x = "Standardized autumn re-growth", y = "Density") +
  theme_classic.adjust +
  xlim(-4,4) +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  theme(axis.title.x = element_text(size = rel(1.5))) +
  theme(legend.title = element_text(size=16),
        legend.text=element_text(size=14),
        legend.position = c(0.9,0.5))



#view plot
print(pheno.plot) 



#Write plot to disk
pdf(file = paste("../results/3-phenoplot-", 
                 populationtitle, 
                 ".pdf", sep = ""), width = 11, height = 8)
print(pheno.plot)
dev.off()
