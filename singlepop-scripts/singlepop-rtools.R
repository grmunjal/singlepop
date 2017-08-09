#------------------------------- 
# Validate Allele Frequencies
# Author: Gitanshu Munjal
# Some aux functions
#------------------------------- 






#----------------------Load packages

#Function that checks for, installs, and loads packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#List of packages used
packages <- c("ggplot2", "gridExtra", "viridis", "lattice", "lme4", "car", "purrr")

#Load packages
ipak(packages)






#----------------------Custom themes for ggplot2
theme_classic.adjust <-   
  theme_classic() + 
  theme(plot.title = element_text(size=20, 
                                  family = "Helvetica", margin = margin(t=10)),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.ticks = element_line(colour = 'black', size=0.5))


theme_classic.bold <- theme_classic.adjust +        
  theme(axis.text = element_text(colour = 'black', 
                                 size=10, face = "bold", family = "Helvetica"),
        axis.title = element_text(colour = 'black', 
                                  size=12, face = "bold", family = "Helvetica"))






#----------------------Colors
llesscolor <- viridis(8, alpha=0.4)[3]
lbasecolor <- viridis(8, alpha=0.6)[3]
lmorecolor <- viridis(8, alpha=1.0)[3]

clesscolor <- viridis(8, alpha=0.4)[5]
cbasecolor <- viridis(8, alpha=0.6)[5]
cmorecolor <- viridis(8, alpha=1.0)[5]

slesscolor <- viridis(8, alpha=0.4)[8]
sbasecolor <- viridis(8, alpha=0.6)[8]
smorecolor <- viridis(8, alpha=1.0)[8]

lesscolor <- viridis(8, alpha=0.2)[3]
basecolor <- viridis(8, alpha=0.6)[3]
morecolor <- viridis(8, alpha=1.0)[3]

plotcolors <- viridis(8)






#----------------------G-matrix function

#Function that makes a G-matrix given a matrix where:
#Every column is a marker and every row is a sample
#The data are allele frequencies

make.G.matrix <- function(M){
  
  #Remove NA markers
  M <- t(na.omit(t(M)))
  
  #Center and standardize marker matrix
  Z <- scale(M, center = TRUE, scale = TRUE)
  Z <- t(na.omit(t(Z)))
  
  #scaling factor for G-matrix (==number of markers)
  K <- ncol(Z)
  
  #G-matrix
  G.matrix <- (Z %*% t(Z) / K)
  
  return(G.matrix)
}






#----------------------Pretty heatmap

prettyHeatmap <- function(M){
  
  p.hmp <- levelplot(M[1:ncol(M), ncol(M):1], 
                     col.regions=viridis(20),
                     xlab="", ylab="")
  
  return(p.hmp)
}






#----------------------PCA for marker data

#Function that does a PCA given a matrix where:
#Every column is a marker and every row is a sample
#The data are allele frequencies

marker.pca <- function(M){
  
  #Remove NA markers
  M <- t(na.omit(t(M)))
  
  #Vector of reference allele frequencies
  P <- round(colMeans(M),3)
  
  #Matrix of reference allele frequencies
  P <- matrix(P, nrow = ncol(M), ncol = nrow(M))
  P <- t(P)
  
  #Center and standardize marker matrix
  Z <- (M - P)/(sqrt(P*(1-P)))
  
  #Peform PCA
  pca.markerdata <- prcomp(Z, center = TRUE, cor = TRUE)
  
  return(pca.markerdata)
}






#----------------------PCA for marker data (returns biplot)

#Function that does a PCA given a matrix where:
#Every column is a marker and every row is a sample
#The data are allele frequencies

prettyBiplot <- function(M){
  
  #Remove NA markers
  M <- t(na.omit(t(M)))
  
  #Vector of reference allele frequencies
  P <- round(colMeans(M),3)
  
  #Matrix of reference allele frequencies
  P <- matrix(P, nrow = ncol(M), ncol = nrow(M))
  P <- t(P)
  
  #Center and standardize marker matrix
  Z <- (M - P)/(sqrt(P*(1-P)))
  
  #Peform PCA
  pca.markerdata <- prcomp(Z, center = TRUE, cor = TRUE)
  
  #Pretty Biplot
  lambda <- pca.markerdata$sdev * sqrt(nrow(pca.markerdata$x))
  biplotdata.all  <- t(t(pca.markerdata$x)/lambda)
  
  biplot.markerdata <-
    ggplot(data.frame(biplotdata.all),
           aes(x=PC1, y=PC2, 
               label=rownames(biplotdata.all))) +
    geom_point(cex=7.7, color="black", shape=1) + 
    geom_text(nudge_x = 0.05)+
    geom_point(cex=7, aes(color=rownames(biplotdata.all))) + 
    scale_color_manual(values=c(
      rep(viridis(3, alpha = 0.5)[3],4), #high
      rep(viridis(3, alpha = 0.5)[1],4),
      rep(viridis(3, alpha = 0.5)[2],4)), #base
      guide=guide_legend(title = NULL)) +
    theme_classic.adjust +
    theme(axis.title.y = element_text(size = rel(1.5)),
          axis.title.x = element_text(size = rel(1.5)),
          legend.text  = element_text(size = rel(1.5)),
          legend.key = element_blank()) 
  
  return(biplot.markerdata)
}











#----------------------Manhattan Plot

#Function that will make a pretty manhattan plot given 
#a data frame with chr, pos, and var (column names should be exactly that)

manhattan.plot.special <- function(man.data, plot.title, thresh){
  
  #Order data by chromosomes and positions within chromosomes
  man.data <- man.data[order(man.data$chr, man.data$pos), ]
  
  #Vector of chromosome names
  chromosomes     <- levels(man.data$chr)[grep("chr*", levels(man.data$chr))]
  
  #Dummy max value for chromosome 1
  chromo.max <- 0
  
  #Useful when coloring the plot alternately for chromosomes
  chrticks <- c()
  
  for (chromosome in 1:length(chromosomes)) {
    
    #Subset chromosome
    chromo.sub <- subset(man.data, man.data$chr==chromosomes[chromosome])
    
    #Add max mapping position from last chromosome
    chromo.sub$pos <- chromo.sub$pos + chromo.max
    
    #Update this chromosome in man.data
    man.data[as.character(rownames(chromo.sub)),] <- chromo.sub
    
    #Update max mapping position value
    chromo.max <- max(chromo.sub$pos)
    
    #Where to place tick for this chromosome
    chr.tick <- round(mean(c(chromo.max, min(chromo.sub$pos))),0)
    
    #Update chromosome tick list
    chrticks <- c(chrticks, chr.tick)
  }
  
  #Alternate colors for chromosomes
  altcolor1 <- viridis(8, alpha = 0.4)[3]
  altcolor2 <- viridis(8, alpha = 0.4)[6]
  
  man.plot <- ggplot(man.data, aes(pos, var)) + 
    geom_point(cex=2, color="black", shape=1) + 
    geom_point(cex=1.6, aes(color=man.data$chr)) + 
    scale_color_manual(values=rep(c(altcolor1, altcolor2), (length(chromosomes)/2)+1)) +
    geom_point(data = na.omit(man.data), aes(pos, var)) +
    geom_abline(intercept = thresh, slope = 0, linetype=2) +
    labs(title=plot.title,x="Position", y="Statistic") + 
    scale_x_continuous(breaks=chrticks, labels=chromosomes) + 
    theme_classic.adjust + theme(legend.position="none", 
                                 plot.title = element_text(size=20, face="bold", 
                                                           hjust = 1)) 
  
  return(man.plot)
}






#----------------------Drift sim
# simulates the trajectory of a Wright-Fisher model with selection
# N: population size (number of chromosomes: N individuals for haploids and N/2 individuals for diploids)
# t: number of generations
# t0: time where the mutation appears in generations (1 by default, appears at the begining)
# j: number of A alleles appearing at time t0 in the population (1 for a de novo mutation by default)
# time sampling: 
# nb_times: number of equaly spaced time points to use between 1 ant t
# N_sample: sample size (in number of chromosomes). Can be a number or a vector of size nb_times

WF_trajectory<-function(N, t, j, t0, nb_times, N_sample, plot=TRUE, sims=100, sleep=1) { 
  
  # generations to draw sample in
  sample_times=round(seq(1,t,length=nb_times))
  
  # initialize the trajectory
  x=numeric(t) 
  x[t0]=j
  
  # set Ne
  cur_N=N    
  
  # loop over generations
  for (i in (t0+1):t)
  {
    p=x[i-1]/(cur_N)
    # calculate sampling probbilities
    prob=p/(p+(1-p))      
    
    # binomial sampling to the next generation
    x[i]=rbinom(1,cur_N,prob)
    
  }
  
  # do time sampling
  sample_freq=rbinom(nb_times,N_sample,x[sample_times]/(cur_N))
  sample_freq=round(sample_freq/N_sample,2)
  
  # plot
  dat <- matrix(nrow=length(sample_times),ncol=2)
  dat[,1] <- c(sample_times)
  dat[,2] <- c(sample_freq)
  if(plot==T)plot(dat, type="l", ylim=c(0,1), xlab="Generations", ylab="Frequency")
  
  # final frequency
  finalfrequencies <- c()
  finalfrequencies[1] <- sample_freq[length(sample_freq)]
  
  if(sims > 1){
    for(k in 2:sims){
      Sys.sleep(sleep)
      
      # loop over generations
      for (i in (t0+1):t)
      {
        p=x[i-1]/(cur_N)
        # calculate sampling probbilities
        prob=p/(p+(1-p))      
        
        # binomial sampling to the next generation
        x[i]=rbinom(1,cur_N,prob)
        
      }
      
      # do time sampling
      sample_freq=rbinom(nb_times,N_sample,x[sample_times]/(cur_N))
      sample_freq=round(sample_freq/N_sample,2)
      
      # plot
      dat <- matrix(nrow=length(sample_times),ncol=2)
      dat[,1] <- c(sample_times)
      dat[,2] <- c(sample_freq)
      if(plot==T)lines(dat, type="l", ylim=c(0,1), xlab="Generations", ylab="Frequency")
      
      # final frequency
      finalfrequencies[k] <- sample_freq[length(sample_freq)]
      
    }
    
  }
  
  # results
  return(finalfrequencies)
}






#----------------------
bulk.descriptions <- list(
  
  CUF = list(
    base = c(1:96),
    high = c(97:192),
    low = c(193:288)
  ),
  
  LAH = list(
    base = c(1:96),
    high = c(97:192),
    low = c(193:288)
  ),
  
  MES = list(
    base = c(1:12),
    high = c(13:48),
    low = c(49:96)
  ),
  
  NOR = list(
    base = c(1:24),
    high = c(25:60),
    low = c(61:96)
  ),
  
  SAR = list(
    base = c(1:96),
    high = c(97:192),
    low = c(193:288)
  ),
  
  WAQ = list(
    base = c(9:40)-8,
    high = c(41:72)-8,
    low = c(73:96)-8
  )
  
  
)






#----------------------
bulk.census <- list(
  
  CUF = list(
    base = 96,
    high = 96,
    low = 96
  ),
  
  LAH = list(
    base = 48,
    high = 96,
    low = 96
  ),
  
  MES = list(
    base = 24, 
    high = 72,
    low = 96
  ),
  
  NOR = list(
    base = 48,
    high = 72,
    low = 72
  ),
  
  SAR = list(
    base = 96,
    high = 96,
    low = 96
  ),
  
  WAQ = list(
    base = 96,
    high = 96,
    low = 72
  )
  
  
)






#----------------------
bulk.bio.descriptions <- list(
  
  CUF = list(
    base1 = 1:24,
    base2 = 1:24 + 24,
    base3 = 1:24 + 24 + 24,
    base4 = 1:24 + 24 + 24 + 24,
    
    high1 = 1:24 + 24 + 24 + 24 + 24,
    high2 = 1:24 + 24 + 24 + 24 + 24 + 24,
    high3 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24,
    high4 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    
    low1 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low2 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low3 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low4 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24
  ),
  
  LAH = list(
    base1 = 1:24,
    base2 = 1:24 + 24,
    base3 = 1:24 + 24 + 24,
    base4 = 1:24 + 24 + 24 + 24,
    
    high1 = 1:24 + 24 + 24 + 24 + 24,
    high2 = 1:24 + 24 + 24 + 24 + 24 + 24,
    high3 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24,
    high4 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    
    low1 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low2 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low3 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low4 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24
  ),
  
  MES = list(
    base1 = 1:12,
    high1 = 1:12 + 12,
    high2 = 1:12 + 12 + 12,
    high3 = 1:12 + 12 + 12 + 12,
    low1 = 1:12 + 12 + 12 + 12 + 12,
    low2 = 1:12 + 12 + 12 + 12 + 12 + 12,
    low3 = 1:12 + 12 + 12 + 12 + 12 + 12 + 12,
    low4 = 1:12 + 12 + 12 + 12 + 12 + 12 + 12 + 12
  ),
  
  NOR = list(
    base1 = 1:12,
    base2 = 1:12 + 12,
    high1 = 1:12 + 12 + 12,
    high2 = 1:12 + 12 + 12 + 12,
    high3 = 1:12 + 12 + 12 + 12 + 12,
    low1 = 1:12 + 12 + 12 + 12 + 12 + 12,
    low2 = 1:12 + 12 + 12 + 12 + 12 + 12 + 12,
    low3 = 1:12 + 12 + 12 + 12 + 12 + 12 + 12 + 12
  ),
  
  SAR = list(
    base1 = 1:24,
    base2 = 1:24 + 24,
    base3 = 1:24 + 24 + 24,
    base4 = 1:24 + 24 + 24 + 24,
    high1 = 1:24 + 24 + 24 + 24 + 24,
    high2 = 1:24 + 24 + 24 + 24 + 24 + 24,
    high3 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24,
    high4 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low1 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low2 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low3 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24,
    low4 = 1:24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24 + 24
  ),
  
  WAQ = list(
    base1 = 1:8,
    base2 = 1:8 + 8,
    base3 = 1:8 + 8 + 8,
    base4 = 1:8 + 8 + 8 + 8,
    high1 = 1:8 + 8 + 8 + 8 + 8,
    high2 = 1:8 + 8 + 8 + 8 + 8 + 8,
    high3 = 1:8 + 8 + 8 + 8 + 8 + 8 + 8,
    high4 = 1:8 + 8 + 8 + 8 + 8 + 8 + 8 + 8,
    low1 = 1:8 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 8,
    low2 = 1:8 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 8,
    low3 = 1:8 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 8
  )
  
  
)






#----------------------Fst
vectorfst.tet <- function(x1, x2, m1, m2){
  
  p1 <- x1/m1
  p2 <- x2/m2
  
  
  #Overall frequency (is ploidy independent)
  po <- (p1+p2)/2
  
  #Frequency of homozygous in overall
  o.ho1 <- po^4
  o.ho2 <- (1-po)^4
  
  #Frequency of homozygous in pop 1
  p1.ho1 <- p1^4
  p1.ho2 <- (1-p1)^4
  
  #Frequency of homozygous in pop 2
  p2.ho1 <- p2^4
  p2.ho2 <- (1-p2)^4
  
  
  ht <- 1 - o.ho1 - o.ho2
  hs <- ((1 - p1.ho1 - p1.ho2) + (1 - p2.ho1 - p2.ho2)) / 2
  Fst.tet <- (ht-hs)/ht
  return(Fst.tet)
}






