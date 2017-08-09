#-------------------------------
# ONE
# SCRIPT
# TO
# RULE 
# THEM
# ALL
#-------------------------------




# Note: find and replace all population titles to the one desired (10 changes)






#Clear work space
dev.off()
rm(list=ls())
setwd("~/singlepop/singlepop-scripts/")
populationtitle <- "CUF"

# Run the first script line by line because it is setup assuming all pops have 
# 4 biological reps (called "sub" in th script). It will throw errors on lines 
# where it is expecting a biorep that does not exist (e.g.: pop only has 2 bioreps)
parse.by.line <- parse(file = "./0-CalculateAlleleFrequencies.R")
for(i in seq_along(parse.by.line)){
  tryCatch(eval(parse.by.line[[i]]),
           error = function(e) message(
             "if this error is a replacement error in pop.sub then it is OK!", 
             as.character(e)))
}






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./1-ValidateAlleleFrequencies.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./2-AFS.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./3-PhenotypeAnalysis.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./4-Pcadapt.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./4c-Fst-nowindow.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./5-EffectivePopSize.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./6-CandidateComparison.R")






#Clear work space
dev.off()
rm(list=ls())                 
setwd("../../singlepop-scripts/")
populationtitle <- "CUF"
source("./7-CandidateMtGenes.R")





print("YAY!")




