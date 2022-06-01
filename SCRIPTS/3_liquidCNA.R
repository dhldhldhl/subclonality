library(pracma); library(ggplot2); library(ggpubr); library(reshape2); 
library(mixtools); library(fitdistrplus); library(dplyr); library(QDNAseq); 
require(gtools); library(gridExtra); library(devtools)
source_url("https://raw.githubusercontent.com/elakatos/liquidCNA/main/mixture_estimation_functions.R")

#read outputs of QDNAseq
load("DATA/bam_by_patient.RData")
p_vec <- 1:length(patient_ids)
p <- p_vec[10]

seg.df <- read.delim(paste0("DATA/2_QDNA_CNout/segment_", patient_ids[p], ".txt"))
cn.df <- read.delim(paste0("DATA/2_QDNA_CNout/raw_cn_", patient_ids[p], ".txt"))

seg.df <- seg.df[,-(1:4)]
cn.df <- cn.df[,-(1:4)]

colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))

#optional
#re-normalise raw and segmented CN values to be centred at CN=2
reNorm <- centreSegs(seg.df)
seg.df <- as.data.frame(t(t(seg.df)/reNorm)*2)
cn.df <- as.data.frame(t(t(cn.df)/reNorm)*2)

#Plot the CN distribution of each sample to gain a quick overview
ggplot(reshape2::melt(seg.df), aes(x=value,y=..scaled.., colour=variable)) +
  geom_density(adjust=1) +
  theme_bw() + scale_x_continuous(limits=c(0.5, 5)) +
  labs(x='Segment copy number',y='Density',colour='')

#generate a dataframe of ensemble segments
#contiguous sections of bins that are constant in ALL samples
segchange <- sapply(1:(nrow(seg.df)-1),
                    function(i) sum(seg.df[i,]!=seg.df[(i+1),])>0)
seg.data <- data.frame(start=c(1,which(segchange)+1),
                       end=c(which(segchange),length(segchange)+1))
seg.data$length <- seg.data$end - seg.data$start+1

#Filter out small segments
seg.sub <- subset(seg.data, length>12) #12 500kb bins
cat('Total number of segments retained: ',nrow(seg.sub))

#Update segment values by fitting a normal distribution to bins in the segment
seg.fit <- getNormalFitSegments(seg.sub, cn.df)

seg.cns <- seg.fit[[1]]
names(seg.cns) <- names(seg.df)

seg.df.upd <- data.frame(matrix(NA, ncol=ncol(seg.df), nrow=nrow(seg.df)))
names(seg.df.upd) <- names(seg.df); row.names(seg.df.upd) <- row.names(seg.df)

for(i in 1:nrow(seg.sub)){
  seg.df.upd[seg.sub$start[i]:seg.sub$end[i],] <- seg.cns[i,]
}

##########################################
########### Purity Estimation ############
##########################################
w = c(0.3,1,0.8,0.15,0.05) #weights of different CN states
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
  x <- na.omit(seg.df.upd[,i])
  pFits <- as.data.frame(sapply(adjVec,
                                function(a) sapply(pVec,
                                                   function(p) evalPurityDensity(p,w,a,x,maxCN))))
  pFits$p <- pVec
  mins <- c(pVec[which.min(apply(pFits[,1:length(adjVec)],1,mean))],
            pVec[which.min(apply(pFits[,1:length(adjVec)],1,median))])
  plF <- ggplot(melt(pFits,id='p'), aes(x=p, y=value)) + geom_point(size=2, alpha=0.3)+
    theme_bw() +
    geom_vline(xintercept = mins, linetype=c('dashed','dotted'), colour=c('red','blue')) +
    labs(title=names(seg.df.upd)[i])
  print(plF)
  pHat.df[,names(seg.df.upd)[i]] <- mins
}

#Correct CN values by purity
#removing CN contamination by normal
#eq being used: CN_corr = 1/p * (CN-2) + 2
pVec <- as.numeric(pHat.df[1,,drop=T])
seg.df.corr <- as.data.frame(t(t(seg.df.upd-2)*1/pVec)+2)
seg.cns.corr <- as.data.frame(t(t(seg.cns-2)*1/pVec)+2)

#remove low purity samples (optional)
# purity_threshold <- 0.05
# abovePurTh <- pVec >= purity_threshold
# seg.df.corr <- seg.df.corr[,abovePurTh]
# seg.cns.corr <- seg.cns.corr[,abovePurTh]

#Plot purity-corrected segment distribution
ggplot(melt(as.data.frame(seg.df.corr)), aes(x=value, colour=variable)) +
  geom_density(adjust=1) + theme_bw() +
  scale_x_continuous(limits=c(0,8)) + geom_vline(xintercept = 1:6) +
  labs(x='Purity-corrected segment CN',y='Density',colour='')

#Designate reference sample
#& calculate dCN
baseSample <- 'Sample1'
seg.dcn <- seg.cns.corr - seg.cns.corr[,baseSample]
seg.dcn.nonBase <- seg.dcn %>% select(-one_of(baseSample))

#Choose samples to investigate (all non-base samples by default) 
colToUse <- names(seg.dcn.nonBase)
#generate all possible permutations
nCol <- length(colToUse)
seg.dcn.toOrder <- seg.dcn.nonBase[,colToUse]
ordVec <- permutations(nCol,nCol,colToUse)

#Set parameters to be used when evaluating sample order and segment monotony
filterMethod <- 'sd'#exclude clonal, sd < theta
cutOffVec <- seq(0.025,0.5,by=0.005) #cut-off values for pre-filtering
epsilon <- 0.05 #error margin for monotony

#Diagnostic plot
fitInfo <- list()
for (cutOff in cutOffVec){
  seg.dcn.Eval <- filterSegmentRatios(seg.dcn.toOrder, cutOff, filterMethod, 0)
  best <- findBestOrder(seg.dcn.Eval, ordVec, epsilon, nCol, 0)
  fitInfo[[as.character(cutOff)]] <- best
}

fitInfo.df <- data.frame(cutOff = cutOffVec,
                         maxFit = sapply(fitInfo, function(x) x$max),
                         segsAboveCut = sapply(fitInfo, function(x) length(x$cons)),
                         segsInOrder = sapply(fitInfo, function(x) max(sapply(x$segs, length))))
p1 <- ggplot(fitInfo.df, aes(x=cutOff, y=maxFit)) + geom_line(size=2,colour='darkgreen') + theme_bw() +
  labs(x='Clonal cut-off',y='Subclonal proportion (of non-clonal segments)')
p2 <- ggplot(fitInfo.df, aes(x=cutOff, y=segsAboveCut)) + geom_line(size=2,colour='deepskyblue3') +
  geom_line(aes(y=segsInOrder), size=2,colour='dodgerblue4') + theme_bw() +
  labs(x='Clonal cut-off',y='Total number of non-clonal and subclonal segments')
grid.arrange(p1,p2,nrow=1)
#p1: proportion of subclonal segments
#p2: total number of non-clonal (light blue) and subclonal (dark blue) segments

#A typical cut-off ~0.1 is adequate
#too low or too high cut-off will select for noise
#ideal cut-off should have a reasonable number of subclonal segments (5-20) 
#while maximises the subclonal proportion

#optimal cut-off for wanted  number of segment
minSegmentNumber <- 12
recommendCutOff <- getCutOffAuto(fitInfo.df, minSegmentNumber)
cat('Recommended cut-off is: ', recommendCutOff)

#plot DeltaCN for chosen cut off + ordering + segment classification
cutOff <- recommendCutOff
fit <- fitInfo[[as.character(cutOff)]]

for (ind in 1:nrow(fit$ord)){
  seg.plot <- seg.dcn[,rev(c(fit$ord[ind,],baseSample))]
  seg.plot$id <- row.names(seg.plot)
  seg.plot$filtered <- seg.plot$id %in% fit$cons
  seg.plot$order <- seg.plot$id %in% fit$segs[[ind]]
  
  p <- ggplot(melt(seg.plot, id=c('id','filtered','order')), aes(x=variable, y=value, group=id, colour=paste0(filtered,' | ',order))) +
    geom_line(size=1.2, alpha=0.75) + theme_bw() +
    scale_colour_manual(values=c('grey70','#487a8b','firebrick3'), labels=c('clonal','unstable','subclonal')) +
    labs(x='',y='Change in segment CN (DeltaCN)',colour='Segment type')
  print(p)
}

#if more than one sample order exists
ordInd <- 1
seg.dcn.toUse <- seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])]


##########################################
######## Derive  subclonal-ratio #########
##########################################

###
#RELATIVE subclonal-ratio estimates
#via segment-by-segment comparison to the highest subclonal-ratio sample
topSample <- names(seg.dcn.toUse)[ncol(seg.dcn.toUse)]
toEstimate <- setdiff(names(seg.dcn.toUse), c(topSample,baseSample))

final.ratios <- estimateRSegmentRatio(seg.dcn.toUse, toEstimate, topSample, 1)

final.medians <- aggregate(final.ratios$value, 
                           by=list(final.ratios$variable), median)

#plot estimated relative ratios
final.medians$xpos <- 1:nrow(final.medians)

ggplot(final.ratios, aes(y=value, x=variable, weight=1)) + geom_violin(fill='firebrick3',alpha=0.2) +
  theme_bw() + scale_y_continuous(limits=c(0,1.2)) +
  geom_jitter(width=0.1, height=0, colour='firebrick3', size=2, alpha=0.6) +
  #geom_point(data=true.df, aes(x=V3, y=V2), colour='firebrick3', size=3) +
  geom_segment(data=final.medians,aes(x=xpos-0.15,xend=xpos+0.15,y=x,yend=x), size=1.5) +
  labs(x='',y=paste0('Subclonal-ratio compared to ',topSample))

#final results table of median relative ratio values
final.results <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
names(final.results) <- c('time','relratio')
final.results$time <- as.character(final.results$time)
final.results[nrow(final.results)+1,] <- c(topSample, 1)

###
#Absolute subclonal-ratio estimates
#by fitting a Gaussian mixture model of constrained means to each sample
#shared mean parameter defines the absolute ratio

# final.results$rat <- NA
# final.results$rat_sd <- NA
# for (samp in final.results$time){
#   final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, 3)
# }
# final.results$rat_sd <- as.numeric(final.results$rat_sd)
# 
# ggplot(final.results, aes(x=as.factor(time), y=rat)) +
#   geom_bar(stat='identity',colour='black',fill='firebrick3',alpha=0.8, width=0.8) +
#   geom_errorbar(aes(ymin=(rat-(1.96*rat_sd))*((rat-(1.96*rat_sd))>0), ymax=rat+(1.96*rat_sd)),
#                 width=.1) +
#   theme_bw() + theme(axis.title.x=element_blank()) +
#   labs(x='Sample', y='Subclonal-ratio')

##########################################
############# Final results ##############
##########################################
final.results[nrow(final.results)+1,] <- c(baseSample, 0,0,0)
final.results$purity_mean <- as.numeric(pHat.df[1,match(final.results$time, 
                                                        names(pHat.df))])
final.results$purity_median <- as.numeric(pHat.df[2,match(final.results$time, 
                                                          names(pHat.df))])

cat('Final result table:\n'); print(final.results)



################################################################################
#STEP A
#get estimated purity for all samples
#read outputs of QDNAseq
load("DATA/bam_by_patient.RData")
p_vec <- 1:length(patient_ids)

get_purity <- function(p){
  print(paste0("Processing ", p, "/80"))
  seg.df <- read.delim(paste0("DATA/2_QDNA_CNout/segment_", patient_ids[p], ".txt"))
  cn.df <- read.delim(paste0("DATA/2_QDNA_CNout/raw_cn_", patient_ids[p], ".txt"))
  
  seg.df <- seg.df[,-(1:4)]
  cn.df <- cn.df[,-(1:4)]
  
  colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
  colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))
  
  #re-normalise raw and segmented CN values to be centred at CN=2
  reNorm <- centreSegs(seg.df)
  seg.df <- as.data.frame(t(t(seg.df)/reNorm)*2)
  cn.df <- as.data.frame(t(t(cn.df)/reNorm)*2)
  
  #generate a dataframe of ensemble segments
  #contiguous sections of bins that are constant in ALL samples
  segchange <- sapply(1:(nrow(seg.df)-1),
                      function(i) sum(seg.df[i,]!=seg.df[(i+1),])>0)
  seg.data <- data.frame(start=c(1,which(segchange)+1),
                         end=c(which(segchange),length(segchange)+1))
  seg.data$length <- seg.data$end - seg.data$start+1
  
  #Filter out small segments
  seg.sub <- subset(seg.data, length>12) #12 500kb bins
  
  #Update segment values by fitting a normal distribution to bins in the segment
  seg.fit <- getNormalFitSegments(seg.sub, cn.df)
  
  seg.cns <- seg.fit[[1]]
  names(seg.cns) <- names(seg.df)
  
  seg.df.upd <- data.frame(matrix(NA, ncol=ncol(seg.df), nrow=nrow(seg.df)))
  names(seg.df.upd) <- names(seg.df); row.names(seg.df.upd) <- row.names(seg.df)
  
  for(i in 1:nrow(seg.sub)){
    seg.df.upd[seg.sub$start[i]:seg.sub$end[i],] <- seg.cns[i,]
  }
  
  # Purity Estimation 
  w = c(0.3,1,0.8,0.15,0.05) #weights of different CN states
  maxCN=8 #assumed maximum CN 
  adjVec = c(0.5,0.6,0.8,0.9,1,1.2,1.3,1.5,1.8,2) #smoothing kernel adjustments
  pVec = seq(0.05, 0.5, by=0.005) #range of purity values to be evaluated
  
  #Estimated optimal purity values stored in pHat.df
  #purity vals that minimise mean and median
  pHat.df <- data.frame(matrix(vector(),ncol=ncol(seg.df.upd), nrow=2))
  names(pHat.df) <- names(seg.df.upd)
  row.names(pHat.df) <- c('mean','median')
  
  for(i in 1:ncol(seg.df.upd)){
    x <- na.omit(seg.df.upd[,i])
    pFits <- as.data.frame(sapply(adjVec,
                                  function(a) sapply(pVec,
                                                     function(p) evalPurityDensity(p,w,a,x,maxCN))))
    pFits$p <- pVec
    mins <- c(pVec[which.min(apply(pFits[,1:length(adjVec)],1,mean))],
              pVec[which.min(apply(pFits[,1:length(adjVec)],1,median))])
    plF <- ggplot(melt(pFits,id='p'), aes(x=p, y=value)) + geom_point(size=2, alpha=0.3)+
      theme_bw() +
      geom_vline(xintercept = mins, linetype=c('dashed','dotted'), colour=c('red','blue')) +
      labs(title=names(seg.df.upd)[i])
    print(plF)
    pHat.df[,names(seg.df.upd)[i]] <- mins
  }
  
  return(pHat.df)
}

all_patient_purities <- vector(mode = "list", length = 80)
for(x in 1:length(p_vec)){
  tryCatch({
    all_patient_purities[[x]] <- get_purity(x)
  }, error=function(e){cat("ERROR at:", x, "\n")})
}
names(all_patient_purities) <- paste0("patient_", patient_ids)

save(all_patient_purities, 
     file = "DATA/all_patient_purities.RData")

#STEP B
load("DATA/all_patient_purities.RData")
binded_all <- bind_cols(all_patient_purities)

plot(density(unlist(binded_all[1,])),
     main = "Density plot of mean purities across all samples \n 
     in all patients")



