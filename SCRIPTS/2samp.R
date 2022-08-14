
library(pracma); library(ggplot2); library(ggpubr); library(reshape2); 
library(mixtools); library(fitdistrplus); library(dplyr); library(QDNAseq); 
require(gtools); library(gridExtra); library(devtools)
source_url("https://raw.githubusercontent.com/elakatos/liquidCNA/main/mixture_estimation_functions.R")
setwd("/Users/dhl/laheen Dropbox/DHL DHL/Cambridge/Internship/subclonality/")
source("SCRIPTS/3_liquidCNA.R")
load("../DATA/bam_by_patient.RData")

################################################################################

sample_num_patient <- sapply(1:80, function(x) length(bam_by_patient[[x]]))
two_samp_patients <- which(sample_num_patient == 2)

patient_ids[two_samp_patients] #patient ids of patients with two time samples

################################################################################

# segs = seg.dcn.toOrder
# cutOff = 0.025
# method = 'sd'
# base = 0

filterSegmentRatios <- function(segs, cutOff, method, base=1){
  #function editted to allow for 2 sample patients
  
  segsExtended <- cbind(segs, base) # add a further column, corresponding to the baseline
  if (!(method %in% c('minmax','sd'))){
    # if method is none of the supported ones, return the original
    print('Filtering criterion not recognised')
    return(segs)
  }
  if (method=='minmax'){
    # compute the minimal and maximal ratio in each segment, and return the ones that are far enough
    thresh.min <- 1-cutOff
    thresh.max <- 1+cutOff
    rat.max <- apply(segsExtended/base,1,max)
    rat.min <- apply(segsExtended/base,1,min)
    return(segs[(rat.max>thresh.max | rat.min<thresh.min),])
  }
  if (method=='sd'){
    # compute the standard deviation of each segment, and return the ones that are varied enough
    rat.sd <- apply(segsExtended,1,sd)
    if(nCol == 1){
      return(data.frame(Sample2 = segs[(rat.sd > cutOff),drop=F]))
    } else {
      return(segs[(rat.sd > cutOff),,drop=F])
    }
  }
}


################################################################################

p_num <- 8

seg.df <- read.delim(paste0("../DATA/2_QDNA_CNout/segment_", patient_ids[p_num], ".txt"))
cn.df <- read.delim(paste0("../DATA/2_QDNA_CNout/raw_cn_", patient_ids[p_num], ".txt"))

seg.df <- seg.df[,-(1:4)]
cn.df <- cn.df[,-(1:4)]

colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))

#re-normalise raw and segmented CN values to be centred at CN=2
reNorm <- centreSegs(seg.df)
seg.df <- as.data.frame(t(t(seg.df)/reNorm)*2)
cn.df <- as.data.frame(t(t(cn.df)/reNorm)*2)

#Plot the CN distribution of each sample to gain a quick overview
explore_plot <- ggplot(reshape2::melt(seg.df), aes(x=value,y=..scaled.., colour=variable)) +
  geom_density(adjust=1) +
  theme_bw() + scale_x_continuous(limits=c(0.5, 5)) +
  labs(x='Segment copy number',y='Density',colour='') + 
  ggtitle(paste0('Patient ', patient_ids[p_num]))
print(explore_plot)

#generate a dataframe of ensemble segments
#contiguous sections of bins that are constant in ALL samples
segchange <- sapply(1:(nrow(seg.df)-1),
                    function(i) sum(seg.df[i,]!=seg.df[(i+1),])>0)
seg.data <- data.frame(start=c(1,which(segchange)+1),
                       end=c(which(segchange),length(segchange)+1))
seg.data$length <- seg.data$end - seg.data$start+1

#Filter out small segments
seg.sub <- subset(seg.data, length>12) #12 500kb bins
cat("Patient",patient_ids[p_num], '| Total number of segments retained: ',nrow(seg.sub))

#Update segment values by fitting a normal distribution to bins in the segment
seg.fit <- getNormalFitSegments(seg.sub, cn.df)

seg.cns <- seg.fit[[1]]
names(seg.cns) <- names(seg.df)

seg.df.upd <- data.frame(matrix(NA, ncol=ncol(seg.df), nrow=nrow(seg.df)))
names(seg.df.upd) <- names(seg.df); row.names(seg.df.upd) <- row.names(seg.df)

for(i in 1:nrow(seg.sub)){
  seg.df.upd[seg.sub$start[i]:seg.sub$end[i],] <- seg.cns[i,]
}

w = c(0.8,1,1,0.15,0.05) #weights of different CN states
maxCN=8 #assumed maximum CN 
adjVec = c(0.5,0.6,0.8,0.9,1,1.2,1.3,1.5,1.8,2) #smoothing kernel adjustments
pVec = seq(0.05, 0.5, by=0.005) #range of purity values to be evaluated

#Estimated optimal purity values stored in pHat.df
#purity vals that minimise mean and median
pHat.df <- data.frame(matrix(vector(),ncol=ncol(seg.df.upd), nrow=2))
names(pHat.df) <- names(seg.df.upd)
row.names(pHat.df) <- c('mean','median')

#The plots show the error of the fit over the range of purity values 
for(i in 1:ncol(seg.df.upd)){
  #for each time sample...
  
  x <- na.omit(seg.df.upd[,i]) #get rid of NA rows (i.e., segments with NAs)
  
  pFits <- as.data.frame(sapply(adjVec,
                                function(a) sapply(pVec,
                                                   function(p) evalPurityDensity(p,w,a,x,maxCN))))
  #each row: the 1 purity value evaluated
  #each column: 1 adjVec, a smoothing kernel adjustment
  #each value: error for given purity/adjVec
  
  pFits$p <- pVec
  mins <- c(pVec[which.min(apply(pFits[,1:length(adjVec)],1,mean))],
            pVec[which.min(apply(pFits[,1:length(adjVec)],1,median))])
  #which row (purity) has minimum mean error and minimum median error?
  
  plF <- ggplot(melt(pFits,id='p'), aes(x=p, y=value)) + geom_point(size=2, alpha=0.3)+
    theme_bw() +
    geom_vline(xintercept = mins, linetype=c('dashed','dotted'), colour=c('red','blue')) +
    labs(title= paste0('Patient ', patient_ids[p_num], " | ", names(seg.df.upd)[i]))
  print(plF)
  pHat.df[,names(seg.df.upd)[i]] <- mins
}

#Correct CN values by purity
#removing CN contamination by normal
#eq being used: CN_corr = 1/p * (CN-2) + 2
pVec <- as.numeric(pHat.df[1,,drop=T])
seg.df.corr <- as.data.frame(t(t(seg.df.upd-2)*1/pVec)+2)
seg.cns.corr <- as.data.frame(t(t(seg.cns-2)*1/pVec)+2)

#Plot purity-corrected segment distribution
purity_plot <- ggplot(melt(as.data.frame(seg.df.corr)), aes(x=value, colour=variable)) +
  geom_density(adjust=1) + theme_bw() +
  scale_x_continuous(limits=c(0,8)) + geom_vline(xintercept = 1:6) +
  labs(x='Purity-corrected segment CN',y='Density',colour='') + 
  ggtitle(paste0('Patient ', patient_ids[p_num]))
print(purity_plot)

#Designate reference sample
#& calculate dCN
baseSample <- names(seg.df.corr)[1]
seg.dcn <- seg.cns.corr - seg.cns.corr[,baseSample]
seg.dcn.nonBase <- seg.dcn %>% select(-one_of(baseSample))

#Choose samples to investigate (all non-base samples by default) 
colToUse <- names(seg.dcn.nonBase)
#generate all possible permutations
nCol <- length(colToUse)
seg.dcn.toOrder <- seg.dcn.nonBase[,colToUse]
ordVec <- permutations(nCol,nCol,colToUse)

#Initial filtering to find best sample order 
#Set parameters to be used when evaluating sample order and segment monotony
filterMethod <- 'sd' #filter by saying... if sd < theta, segment = clonal
cutOffVec <- seq(0.025,0.5,by=0.005) #thetas to try for initial-filtering
epsilon <- 0.05 #error margin for monotony

#Diagnostic plot
#at different thetas shows:
#the proportion of subclonal segments (of non-clonal segments)
#total number of non-clonal (light blue)
#and subclonal (dark blue) segments
fitInfo <- list()
for (cutOff in cutOffVec){
  seg.dcn.Eval <- filterSegmentRatios(seg.dcn.toOrder, cutOff, filterMethod, 0)
  #Filter relative segment values (compared to baseline) to discard clonal 
  #(nonchanging) segments that stay relatively constant
  #Returns: a filtered set of a segments that have values exceeding the threshold
  
  best <- findBestOrder(seg.dcn.Eval, ordVec, epsilon, nCol, 0)
  #Input: segments where clonal segments are initially filtered
  #Returns: 
  #$cons: row-names of all segments considered, 
  #max: maximum fit gained, 
  #$ord: order with the maximum fit, 
  #$segs:segments monotone according to max fit
  fitInfo[[as.character(cutOff)]] <- best
}

#optimal cut-off for wanted  number of segment
minSegmentNumber <- 9
recommendCutOff <- getCutOffAuto(fitInfo.df, minSegmentNumber)
cat('Recommended cut-off is: ', recommendCutOff)

#plot DeltaCN for chosen cut off + ordering + segment classification
cutOff <- 0.025
fit <- fitInfo[[as.character(cutOff)]]

for (ind in 1:nrow(fit$ord)){
  seg.plot <- seg.dcn[,rev(c(fit$ord[ind,],baseSample))]
  seg.plot$id <- row.names(seg.plot)
  seg.plot$filtered <- seg.plot$id %in% fit$cons
  seg.plot$order <- seg.plot$id %in% fit$segs[[ind]]
  
  p <- ggplot(melt(seg.plot, id=c('id','filtered','order')), aes(x=variable, y=value, group=id, colour=paste0(filtered,' | ',order))) +
    geom_line(size=1.2, alpha=0.75) + theme_bw() +
    scale_colour_manual(values=c('grey70','#487a8b','firebrick3'), labels=c('clonal','unstable','subclonal')) +
    labs(x='',y='Change in segment CN (DeltaCN)',colour='Segment type') + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  print(p)
}

#if more than one sample order exists
ordInd <- 1
seg.dcn.toUse <- data.frame(Sample2 = seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])])

final.results = data.frame(time = "Sample2", relratio = NA)

###
#Absolute subclonal-ratio estimates
#by fitting a Gaussian mixture model of constrained means to each sample
#shared mean parameter defines the absolute ratio
# 
final.results$rat <- NA
final.results$rat_sd <- NA
for (samp in final.results$time){
  final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, 3)
}
final.results$rat_sd <- as.numeric(final.results$rat_sd)

ggplot(final.results, aes(x=as.factor(time), y=rat)) +
  geom_bar(stat='identity',colour='black',fill='firebrick3',alpha=0.8, width=0.8) +
  geom_errorbar(aes(ymin=(rat-(1.96*rat_sd))*((rat-(1.96*rat_sd))>0), ymax=rat+(1.96*rat_sd)),
                width=.1) +
  theme_bw() + theme(axis.title.x=element_blank()) +
  labs(x='Sample', y='Subclonal-ratio') +
  ggtitle(paste0('Patient ', patient_ids[p_num]))

final.results[nrow(final.results)+1,] <- c(baseSample, 0,0,0)
final.results$purity_mean <- as.character(pHat.df[1,match(final.results$time, 
                                                          names(pHat.df))])
final.results$purity_median <- as.character(pHat.df[2,match(final.results$time, 
                                                            names(pHat.df))])


################################################################################
rat_cutOff <- function(cut){
  cutOff <- cut
  fit <- fitInfo[[as.character(cutOff)]]
  
  for (ind in 1:nrow(fit$ord)){
    seg.plot <- seg.dcn[,rev(c(fit$ord[ind,],baseSample))]
    seg.plot$id <- row.names(seg.plot)
    seg.plot$filtered <- seg.plot$id %in% fit$cons
    seg.plot$order <- seg.plot$id %in% fit$segs[[ind]]
    
    p <- ggplot(melt(seg.plot, id=c('id','filtered','order')), aes(x=variable, y=value, group=id, colour=paste0(filtered,' | ',order))) +
      geom_line(size=1.2, alpha=0.75) + theme_bw() +
      scale_colour_manual(values=c('grey70','#487a8b','firebrick3'), labels=c('clonal','unstable','subclonal')) +
      labs(x='',y='Change in segment CN (DeltaCN)',colour='Segment type') + 
      ggtitle(paste0('Patient ', patient_ids[p_num]))
    print(p)
  }
  
  #if more than one sample order exists
  ordInd <- 1
  seg.dcn.toUse <- data.frame(Sample2 = seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])])
  
  final.results = data.frame(time = "Sample2", relratio = NA)
  
  ###
  #Absolute subclonal-ratio estimates
  #by fitting a Gaussian mixture model of constrained means to each sample
  #shared mean parameter defines the absolute ratio
  # 
  final.results$rat <- NA
  final.results$rat_sd <- NA
  for (samp in final.results$time){
    final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, 3)
  }
  return(final.results$rat)
}

rats <- rep(NULL, 96)

cutOffVec

for(x in 1:length(cutOffVec)){
  tryCatch({
    cat("\n Starting", cutOffVec[x], "\n")
    rats[x] <- rat_cutOff(cutOffVec[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

table(rats)

################################################################################
save(eg_2samp_final.results, rats, res1, res2, 
     res1b, res2b, file = "../DATA/weekly_june4_data")

# what does a negative absolute subclonal ratio mean?
# what would be the smartest way to choose a cutOff











