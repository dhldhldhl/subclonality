# library(pracma); library(ggplot2); library(ggpubr); library(reshape2); 
# library(mixtools); library(fitdistrplus); library(dplyr); library(QDNAseq); 
# require(gtools); library(gridExtra); library(devtools)
# source_url("https://raw.githubusercontent.com/elakatos/liquidCNA/main/mixture_estimation_functions.R")

################################################################################

# #setwd
# wd <- "/Users/dhl/laheen Dropbox/DHL DHL/Cambridge/Internship/subclonality"
# setwd(wd)

################################################################################
#load("../DATA/time.RData")

em_steps <- 1500

run_liquidCNA <- function(p_num, timely = FALSE, batch1 = TRUE, nbatch = NULL){
  #timely: == TRUE is patient has time samples > 6
  #batch1: == TURE is first batch, == FALSE is second batch
  #nbatch: is a number indicating number of time samples in the batch
  
  seg.df <- read.delim(paste0("../DATA/2_QDNA_CNout/segment_", patient_ids[p_num], ".txt"))
  cn.df <- read.delim(paste0("../DATA/2_QDNA_CNout/raw_cn_", patient_ids[p_num], ".txt"))
  
  seg.df <- seg.df[,-(1:4)]
  cn.df <- cn.df[,-(1:4)]
  
  colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
  colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))
  
  #####
  # # option1
  # if(timely == TRUE){
  #   if(batch1 == TRUE){
  #     seg.df <- seg.df[,1:nbatch]
  #     cn.df <- cn.df[,1:nbatch]
  #   } else {
  #     tail.id <- tail(1:ncol(seg.df), nbatch)
  #     col.id <- c(base.id,tail.id)
  #     seg.df <- seg.df[,col.id]
  #     cn.df <- cn.df[,col.id]
  #   }
  # }
  
  # # option2
  # if(timely == TRUE){
  #   if(batch1 == TRUE){
  #     seg.df <- seg.df[,1:nbatch]
  #     cn.df <- cn.df[,1:nbatch]
  #   } else {
  #     tail.id <- tail(1:ncol(seg.df), nbatch)
  #     col.id <- c(base.id,tail.id)
  #     seg.df <- seg.df[,col.id]
  #     cn.df <- cn.df[,col.id]
  #   }
  # }
  
  # # option 3
  if(timely == TRUE){
    if(batch1 == TRUE){
      seg.df <- seg.df[,1:nbatch]
      cn.df <- cn.df[,1:nbatch]
    } else {
      tail.id <- tail(1:ncol(seg.df), nbatch)
      col.id <- c(base.id,tail.id)
      seg.df <- seg.df[,col.id]
      cn.df <- cn.df[,col.id]
    }
  }
  #####
  
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
  
  ##########################################
  ########### Purity Estimation ############
  ##########################################
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
  
  #Designate reference sample & calculate dCN
  
  #baseSample has been set such that sample with minimum CNS is baseSample
  means.cns <- sapply(1:ncol(seg.cns.corr), function(x) mean(seg.cns.corr[,x]))
  baseIndex <- match(min(means.cns), means.cns)
  baseSample <- names(seg.df.corr)[baseIndex]
  if(timely == TRUE & batch1 == FALSE){
    #unless its the second batch for timely patients
    baseSample <- names(seg.df.corr)[1] 
  }
  
  seg.dcn <- seg.cns.corr - seg.cns.corr[,baseSample]
  seg.dcn.nonBase <- seg.dcn %>% select(-one_of(baseSample))
  
  #Choose samples to investigate (all non-base samples) 
  colToUse <- names(seg.dcn.nonBase)
  #generate all possible permutations
  nCol <- length(colToUse)
  seg.dcn.toOrder <- seg.dcn.nonBase[,colToUse]
  ordVec <- permutations(nCol,nCol,colToUse)
  
  #Initial filtering to find best sample order 
  #Set parameters to be used when evaluating sample order and segment monotony
  filterMethod <- 'sd' #filter by saying... if sd < theta, segment = clonal
  cutOffVec <- seq(0.025,0.3,by=0.005) #thetas to try for initial-filtering (of clonal vs non-clonal)
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
  
  fitInfo.df <- data.frame(cutOff = cutOffVec,
                           maxFit = sapply(fitInfo, function(x) x$max),
                           segsAboveCut = sapply(fitInfo, function(x) length(x$cons)),
                           segsInOrder = sapply(fitInfo, function(x) max(sapply(x$segs, length))))
  p1 <- ggplot(fitInfo.df, aes(x=cutOff, y=maxFit)) + geom_line(size=2,colour='darkgreen') + theme_bw() +
    labs(x='Clonal cut-off',y='Subclonal proportion (of non-clonal segments)') + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  p2 <- ggplot(fitInfo.df, aes(x=cutOff, y=segsAboveCut)) + geom_line(size=2,colour='deepskyblue3') +
    geom_line(aes(y=segsInOrder), size=2,colour='dodgerblue4') + theme_bw() +
    labs(x='Clonal cut-off',y='Total number of non-clonal and subclonal segments') + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  grid.arrange(p1,p2,nrow=1)
  #p1: proportion of subclonal segments
  #p2: total number of non-clonal (light blue) and subclonal (dark blue) segments
  
  #A typical cut-off ~0.1 is adequate
  #too low or too high cut-off will select for noise
  #ideal cut-off should have a reasonable number of subclonal segments (5-20) 
  #while maximises the subclonal proportion
  
  #optimal cut-off for wanted  number of segment
  minSegmentNumber <- 9
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
    
    p <- ggplot(melt(seg.plot, id=c('id','filtered','order')), 
                aes(x=variable, y=value, group=id, colour=paste0(filtered,' | ',order))) +
      geom_line(size=1.2, alpha=0.75) + theme_bw() +
      scale_colour_manual(values=c('grey70','#487a8b','firebrick3'), 
                          labels=c('clonal','unstable','subclonal')) +
      labs(x='',y='Change in segment CN (DeltaCN)',colour='Segment type') + 
      ggtitle(paste0('Patient ', patient_ids[p_num]))
    print(p)
  }
  
  #if more than one sample order exists
  ordInd <- 1
  seg.dcn.toUse <- seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])]
  
  big.dCN <- sapply(1:ncol(seg.dcn.toUse), function(x) which(seg.dcn.toUse[,x]>4))
  big.dCN <- unlist(big.dCN)
  if(length(big.dCN) > 0){
    seg.dcn.toUse <- seg.dcn.toUse[-unique(big.dCN),]
  }
  ##########################################
  ######## Derive  subclonal-ratio #########
  ##########################################
  
  estimateRGaussianFit <- function(seg.rel.toUse, topSample, final.medians, gaussSigma=0.1, nStates=3, lambdaOverwrite=0.2, verbose=F){
    seg.top <- seg.rel.toUse[,topSample]
    start <- mean(seg.top[seg.top<0])
    if (is.nan(start)){
      # estimate the initial value for fitting from segments with DCN between 0 and 1
      start <- -1* mean(seg.top[seg.top>0 & seg.top<1])
    }
    if (is.nan(start)){ # if no starting point can be found, start from a default of 25%
      start <- -0.25
    }
    # By default the following DCN values are used: -2, -1, 1, 2, 3
    mixfit <- normalmixEM(as.numeric(seg.top),
                          lambda = lambdaOverwrite,
                          mean.constr = c("2a","a","-a","-2a","-3a"),
                          mu = c(2,1,-1,-2,-3)*start,
                          sigma=gaussSigma)
    # only return the obtained value if converged in less than 800 iterations
    if (length(mixfit$all.loglik)<em_steps){
      final.medians$rat[final.medians$time==topSample] <- -1*(mixfit$mu[2])
      #sigma depends on the number of segments used and that depends on how many states they were distributed across
      final.medians$rat_sd[final.medians$time==topSample] <- (mixfit$sigma[1])/sqrt(nrow(seg.rel.toUse)/nStates)
    }
    if (verbose){
      print(mixfit$lambda)
    }
    return(final.medians)
  }
  
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
    labs(x='',y=paste0('Subclonal-ratio compared to ',topSample)) + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  
  #final results table of median relative ratio values
  final.results <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
  names(final.results) <- c('time','relratio')
  final.results$time <- as.character(final.results$time)
  final.results[nrow(final.results)+1,] <- c(topSample, 1)
  
  ###
  #Absolute subclonal-ratio estimates
  #by fitting a Gaussian mixture model of constrained means to each sample
  #shared mean parameter defines the absolute ratio
  # 
  final.results$rat <- NA
  final.results$rat_sd <- NA
  for(samp in final.results$time){
    tryCatch({
      ##################
      final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, nStates = 3)
      ##################
    }, error=function(e){cat("\n ERROR at absolute estimation \n")})
  }
  final.results$rat_sd <- as.numeric(final.results$rat_sd)
  
  ggplot(final.results, aes(x=as.factor(time), y=rat)) +
    geom_bar(stat='identity',colour='black',fill='firebrick3',alpha=0.8, width=0.8) +
    geom_errorbar(aes(ymin=(rat-(1.96*rat_sd))*((rat-(1.96*rat_sd))>0), ymax=rat+(1.96*rat_sd)),
                  width=.1) +
    theme_bw() + theme(axis.title.x=element_blank()) +
    labs(x='Sample', y='Subclonal-ratio') +
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  
  ##########################################
  ############# Final results ##############
  ##########################################
  # final.results <- rbind(c(baseSample, 0,0,0), final.results)
  final.results[nrow(final.results)+1,] <- c(baseSample, 0,0,0)
  final.results$purity_mean <- as.character(pHat.df[1,match(final.results$time, 
                                                            names(pHat.df))])
  final.results$purity_median <- as.character(pHat.df[2,match(final.results$time, 
                                                              names(pHat.df))])
  final.results$seg_used <- nrow(seg.dcn.toUse)
  final.results$cutOff <- recommendCutOff
  return(final.results)
}

old_run_liquidCNA <- function(p_num, timely = FALSE, batch1 = TRUE, nbatch = NULL){
  #timely: == TRUE is patient has time samples > 6
  #batch1: == TURE is first batch, == FALSE is second batch
  #nbatch: is a number indicating number of time samples in the batch
  
  seg.df <- read.delim(paste0("../DATA/2_QDNA_CNout/segment_", patient_ids[p_num], ".txt"))
  cn.df <- read.delim(paste0("../DATA/2_QDNA_CNout/raw_cn_", patient_ids[p_num], ".txt"))
  
  seg.df <- seg.df[,-(1:4)]
  cn.df <- cn.df[,-(1:4)]
  
  colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
  colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))
  
  #####
  # # option 0 
  # if(timely == TRUE){
  #   if(batch1 == TRUE){
  #     seg.df <- seg.df[,1:nbatch]
  #     cn.df <- cn.df[,1:nbatch]
  #   } else {
  #     seg.df <- seg.df[,tail(1:ncol(seg.df), nbatch)]
  #     cn.df <- cn.df[,tail(1:ncol(cn.df), nbatch)]
  #   }
  # }
  
  # # option1
  # if(timely == TRUE){
  #   if(batch1 == TRUE){
  #     seg.df <- seg.df[,1:nbatch]
  #     cn.df <- cn.df[,1:nbatch]
  #   } else {
  #     tail.id <- tail(1:ncol(seg.df), nbatch)
  #     col.id <- c(base.id,tail.id)
  #     seg.df <- seg.df[,col.id]
  #     cn.df <- cn.df[,col.id]
  #   }
  # }
  
  # option2
  # if(timely == TRUE){
  #   if(batch1 == TRUE){
  #     seg.df <- seg.df[,1:nbatch]
  #     cn.df <- cn.df[,1:nbatch]
  #   } else {
  #     tail.id <- tail(1:ncol(seg.df), nbatch)
  #     col.id <- c(base.id,tail.id)
  #     seg.df <- seg.df[,col.id]
  #     cn.df <- cn.df[,col.id]
  #   }
  # }
  
  # # option 3
  if(timely == TRUE){
    base.id <- 1
    if(batch1 == TRUE){
      seg.df <- seg.df[,1:nbatch]
      cn.df <- cn.df[,1:nbatch]
    } else {
      tail.id <- tail(1:ncol(seg.df), nbatch)
      col.id <- c(base.id,tail.id)
      seg.df <- seg.df[,col.id]
      cn.df <- cn.df[,col.id]
    }
  }
  #####
  
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
  
  ##########################################
  ########### Purity Estimation ############
  ##########################################
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
  
  #baseSample has been reset such that sample with minimum CNS is baseSample
  #baseSample <- names(seg.df.corr)[1] #this was previous
  means.cns <- sapply(1:ncol(seg.cns.corr), function(x) mean(seg.cns.corr[,x]))
  baseIndex <- match(min(means.cns), means.cns)
  baseSample <- names(seg.df.corr)[baseIndex]
  
  seg.dcn <- seg.cns.corr - seg.cns.corr[,baseSample]
  seg.dcn.nonBase <- seg.dcn %>% select(-one_of(baseSample))
  
  #Choose samples to investigate (all non-base samples by default) 
  colToUse <- names(seg.dcn.nonBase)
  #generate all possible permutations
  nCol <- length(colToUse)
  seg.dcn.toOrder <- seg.dcn.nonBase[,colToUse]
  ordVec <- permutations(nCol,nCol,colToUse)
  
  #Permutation is too long
  # if( nrow(ordVec) > 1000 ){
  # 
  #   chrono <- which(apply(ordVec, 1, function(x) identical(x, colToUse)))
  #   nonchrono <- (1:nrow(ordVec))[-chrono]
  #   rows <- sample(nonchrono, 3000, replace = F)
  #   subOrdVec <- ordVec[rows,]
  #   subOrdVec <- rbind(subOrdVec, ordVec[chrono,])
  #   ordVec <- subOrdVec
  # 
  # }
  
  #Initial filtering to find best sample order 
  #Set parameters to be used when evaluating sample order and segment monotony
  filterMethod <- 'sd' #filter by saying... if sd < theta, segment = clonal
  cutOffVec <- seq(0.025,0.3,by=0.005) #thetas to try for initial-filtering
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
  
  fitInfo.df <- data.frame(cutOff = cutOffVec,
                           maxFit = sapply(fitInfo, function(x) x$max),
                           segsAboveCut = sapply(fitInfo, function(x) length(x$cons)),
                           segsInOrder = sapply(fitInfo, function(x) max(sapply(x$segs, length))))
  p1 <- ggplot(fitInfo.df, aes(x=cutOff, y=maxFit)) + geom_line(size=2,colour='darkgreen') + theme_bw() +
    labs(x='Clonal cut-off',y='Subclonal proportion (of non-clonal segments)') + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  p2 <- ggplot(fitInfo.df, aes(x=cutOff, y=segsAboveCut)) + geom_line(size=2,colour='deepskyblue3') +
    geom_line(aes(y=segsInOrder), size=2,colour='dodgerblue4') + theme_bw() +
    labs(x='Clonal cut-off',y='Total number of non-clonal and subclonal segments') + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  grid.arrange(p1,p2,nrow=1)
  #p1: proportion of subclonal segments
  #p2: total number of non-clonal (light blue) and subclonal (dark blue) segments
  
  #A typical cut-off ~0.1 is adequate
  #too low or too high cut-off will select for noise
  #ideal cut-off should have a reasonable number of subclonal segments (5-20) 
  #while maximises the subclonal proportion
  
  #optimal cut-off for wanted  number of segment
  minSegmentNumber <- 9
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
    
    p <- ggplot(melt(seg.plot, id=c('id','filtered','order')), 
                aes(x=variable, y=value, group=id, colour=paste0(filtered,' | ',order))) +
      geom_line(size=1.2, alpha=0.75) + theme_bw() +
      scale_colour_manual(values=c('grey70','#487a8b','firebrick3'), 
                          labels=c('clonal','unstable','subclonal')) +
      labs(x='',y='Change in segment CN (DeltaCN)',colour='Segment type') + 
      ggtitle(paste0('Patient ', patient_ids[p_num]))
    print(p)
  }
  
  #if more than one sample order exists
  ordInd <- 1
  seg.dcn.toUse <- seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])]
  
  
  ##########################################
  ######## Derive  subclonal-ratio #########
  ##########################################
  
  estimateRGaussianFit <- function(seg.rel.toUse, topSample, final.medians, gaussSigma=0.1, nStates=3, lambdaOverwrite=0.2, verbose=F){
    seg.top <- seg.rel.toUse[,topSample]
    start <- mean(seg.top[seg.top<0])
    if (is.nan(start)){
      # estimate the initial value for fitting from segments with DCN between 0 and 1
      start <- -1* mean(seg.top[seg.top>0 & seg.top<1])
    }
    if (is.nan(start)){ # if no starting point can be found, start from a default of 25%
      start <- -0.25
    }
    # By default the following DCN values are used: -2, -1, 1, 2, 3
    mixfit <- normalmixEM(as.numeric(seg.top),
                          lambda = lambdaOverwrite,
                          mean.constr = c("2a","a","-a","-2a","-3a"),
                          mu = c(2,1,-1,-2,-3)*start,
                          sigma=gaussSigma)
    # only return the obtained value if converged in less than 800 iterations
    if (length(mixfit$all.loglik)<em_steps){
      final.medians$rat[final.medians$time==topSample] <- -1*(mixfit$mu[2])
      #sigma depends on the number of segments used and that depends on how many states they were distributed across
      final.medians$rat_sd[final.medians$time==topSample] <- (mixfit$sigma[1])/sqrt(nrow(seg.rel.toUse)/nStates)
    }
    if (verbose){
      print(mixfit$lambda)
    }
    return(final.medians)
  }
  
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
    labs(x='',y=paste0('Subclonal-ratio compared to ',topSample)) + 
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  
  #final results table of median relative ratio values
  final.results <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
  names(final.results) <- c('time','relratio')
  final.results$time <- as.character(final.results$time)
  final.results[nrow(final.results)+1,] <- c(topSample, 1)
  
  ###
  #Absolute subclonal-ratio estimates
  #by fitting a Gaussian mixture model of constrained means to each sample
  #shared mean parameter defines the absolute ratio
  # 
  final.results$rat <- NA
  final.results$rat_sd <- NA
  for(samp in final.results$time){
    tryCatch({
      ##################
      final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, nStates = 3)
      ##################
    }, error=function(e){cat("\n ERROR at absolute estimation \n")})
  }
  final.results$rat_sd <- as.numeric(final.results$rat_sd)
  
  ggplot(final.results, aes(x=as.factor(time), y=rat)) +
    geom_bar(stat='identity',colour='black',fill='firebrick3',alpha=0.8, width=0.8) +
    geom_errorbar(aes(ymin=(rat-(1.96*rat_sd))*((rat-(1.96*rat_sd))>0), ymax=rat+(1.96*rat_sd)),
                  width=.1) +
    theme_bw() + theme(axis.title.x=element_blank()) +
    labs(x='Sample', y='Subclonal-ratio') +
    ggtitle(paste0('Patient ', patient_ids[p_num]))
  
  ##########################################
  ############# Final results ##############
  ##########################################
  # final.results <- rbind(c(baseSample, 0,0,0), final.results)
  final.results[nrow(final.results)+1,] <- c(baseSample, 0,0,0)
  final.results$purity_mean <- as.character(pHat.df[1,match(final.results$time, 
                                                            names(pHat.df))])
  final.results$purity_median <- as.character(pHat.df[2,match(final.results$time, 
                                                              names(pHat.df))])
  final.results$seg_used <- nrow(seg.dcn.toUse)
  final.results$cutOff <- recommendCutOff
  return(final.results)
}

run_2samp_liquidCNA <- function(p_num){
  cut <- 0.16
  #update filterSegmentRatios() for 2 time patients
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
  means.cns <- sapply(1:ncol(seg.cns.corr), function(x) mean(seg.cns.corr[,x]))
  baseIndex <- match(min(means.cns), means.cns)
  baseSample <- names(seg.df.corr)[baseIndex]
  seg.dcn <- seg.cns.corr - seg.cns.corr[,baseSample]
  seg.dcn.nonBase <- seg.dcn %>% select(-one_of(baseSample))
  nonBaseSample <- colnames(seg.dcn.nonBase)
  
  #Choose samples to investigate (all non-base samples by default) 
  colToUse <- names(seg.dcn.nonBase)
  #generate all possible permutations
  nCol <- length(colToUse)
  seg.dcn.toOrder <- seg.dcn.nonBase[,colToUse]
  ordVec <- permutations(nCol,nCol,colToUse)
  
  #Initial filtering to find best sample order 
  #Set parameters to be used when evaluating sample order and segment monotony
  filterMethod <- 'sd' #filter by saying... if sd < theta, segment = clonal
  cutOffVec <- seq(0,0.35,by=0.005) #thetas to try for initial-filtering
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
  # recommendCutOff <- getCutOffAuto(fitInfo.df, minSegmentNumber)
  # cat('Recommended cut-off is: ', recommendCutOff)
  
  #plot DeltaCN for chosen cut off + ordering + segment classification
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
  seg.dcn.toUse <- data.frame(nonBaseSample = seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])])
  colnames(seg.dcn.toUse) <- nonBaseSample
  
  big.dCN <- sapply(1:ncol(seg.dcn.toUse), function(x) which(seg.dcn.toUse[,x]>4))
  big.dCN <- unlist(big.dCN)
  if(length(big.dCN) > 0){
    seg.dcn.toUse <- seg.dcn.toUse[-unique(big.dCN),]
  }
  seg.dcn.toUse <- data.frame(seg.dcn.toUse)
  colnames(seg.dcn.toUse) <- nonBaseSample
  final.results = data.frame(time = nonBaseSample, relratio = NA)
  
  ###
  #Absolute subclonal-ratio estimates
  #by fitting a Gaussian mixture model of constrained means to each sample
  #shared mean parameter defines the absolute ratio
  # 
  
  estimateRGaussianFit <- function(seg.rel.toUse, topSample, final.medians, gaussSigma=0.1, nStates=3, lambdaOverwrite=0.2, verbose=F){
    seg.top <- seg.rel.toUse[,topSample]
    start <- mean(seg.top[seg.top<0])
    if (is.nan(start)){
      # estimate the initial value for fitting from segments with DCN between 0 and 1
      start <- -1* mean(seg.top[seg.top>0 & seg.top<1])
    }
    if (is.nan(start)){ # if no starting point can be found, start from a default of 25%
      start <- -0.25
    }
    # By default the following DCN values are used: -2, -1, 1, 2, 3
    mixfit <- normalmixEM(as.numeric(seg.top),
                          lambda=lambdaOverwrite,
                          mean.constr = c("2a","a","-a","-2a","-3a"),
                          mu=c(2,1,-1,-2,-3)*start,sigma=gaussSigma)
    # only return the obtained value if converged in less than 800 iterations
    if (length(mixfit$all.loglik)<em_steps){
      final.medians$rat[final.medians$time==topSample] <- -1*(mixfit$mu[2])
      #sigma depends on the number of segments used and that depends on how many states they were distributed across
      final.medians$rat_sd[final.medians$time==topSample] <- (mixfit$sigma[1])/sqrt(nrow(seg.rel.toUse)/nStates)
    }
    if (verbose){
      print(mixfit$lambda)
    }
    return(final.medians)
  }
  
  final.results$rat <- NA
  final.results$rat_sd <- NA
  
  for (samp in final.results$time){
    final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, nStates = 2)
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
  final.results$seg_used <- nrow(seg.dcn.toUse)
  final.results$cutOff <- cut
  return(final.results)
}

plot.3C <- function(timeSample, cn.df, seg.df, seg.sub, seg.plot, patient_id = NULL){
  p.id <- patient_id
  
  ts <- timeSample
  
  ################################################
  #Which of all ensemble segments clonal, unstable or sub-clonal?
  
  #subclonal
  subclonal.seg <- seg.plot[seg.plot$filtered == TRUE & seg.plot$order == TRUE,]
  subclonal.seg.id <- as.numeric(rownames(subclonal.seg))
  
  #unstable
  unstable.seg <- seg.plot[seg.plot$filtered == TRUE & seg.plot$order == FALSE,]
  unstable.seg.id <- as.numeric(rownames(unstable.seg))
  
  #clonal
  clonal.seg <- seg.plot[seg.plot$filtered == FALSE & seg.plot$order == FALSE,]
  clonal.seg.id <- as.numeric(rownames(clonal.seg))
  
  ################################################
  #cns.of.bins() to get CNS of bins involved in the three segment types
  
  cns.of.bins <- function(x.seg.id){
    #requires: cn.df; seg.df; seg.sub; seg.plot
    #Input: x.seg.id <- id (row number) out of all ensemble segments (seg.sub)  for clonal, unstable or sub-clonal segments
    #Output: cns of the bins of segments that are either clonal, unstable or sub-clonal.
    
    #make an empty vector length of the genome (total number of bins)
    bins <- rep(NA, nrow(cn.df))
    
    if(length(x.seg.id) == 0){
      return(bins)
    } else{
      #fill in the bin
      for(x in 1:length(x.seg.id)){
        seg.loci.df <- seg.sub[x.seg.id,]
        start <- as.numeric(seg.loci.df[x,][1])
        end <- as.numeric(seg.loci.df[x,][2])
        bins[start:end] <- seg.df[start,ts]
      }
      return(bins)
    }
  }
  
  subclonal.bins.cns <- cns.of.bins(subclonal.seg.id)
  unstable.bins.cns <- cns.of.bins(unstable.seg.id)
  clonal.bins.cns <- cns.of.bins(clonal.seg.id)
  
  ################################################
  #Plot
  
  plot(cn.df[,ts], col = "gray", xlab = "Bin", ylab = "CN states",
       main= paste0("Segments used by LiquidCNA to calculate subclonality. \n Sample", ts," from Patient ", p.id))
  points(seg.df[,ts], col = "lightblue")
  points(unstable.bins.cns, col = "blue")
  points(clonal.bins.cns, col = "black")
  points(subclonal.bins.cns, col = "firebrick3")
  legend("topright", legend=c("sub-clonal", "unstable", "clonal", "segment-wise", "bin-wise"), 
         col=c("firebrick3", "blue", "black", "lightblue", "gray"), pch=21, cex=0.8)
}

append.RECIST <- function(){
  #RECIST$ichorCNA <- NA
  RECIST$rat <- NA
  RECIST$purity_mean <- NA
  
  #Fixing error in dates
  RECIST[RECIST$Patient_ID==3357,]$Date[1] <- "2011-11-29"
  RECIST[RECIST$Patient_ID==1483,]$Date[3] <- "2011-12-01"
  
  for(p in 1:length(patient_ids)){
    tryCatch({
      #get one patient
      p.id <- patient_ids[p]
      #subset of ichorCNA for this patient
      patient.ichor <- ichorCNA[ichorCNA$Patient_ID == p.id,]
      patient.ichor <- patient.ichor[order(patient.ichor$Sample_ID, decreasing = F),] #order by date
      patient.ichor.dates <- patient.ichor$Date
      
      #subset of RECIST for this patient
      patient.RECIST <- RECIST[RECIST$Patient_ID == p.id,]
      patient.RECIST.date <- patient.RECIST$Blood.draw
      
      patient.RECIST.date.new <- structure(numeric(length(patient.RECIST.date)), class="Date")
      
      for(x in 1:length(patient.RECIST.date)){
        new.date <- as.Date(patient.RECIST.date[x], format = "%m/%d/%y")
        patient.RECIST.date.new[x] <- new.date
      }
      
      #For which timeSamples do we have RECIST data?
      samp.w.RECIST <- match(patient.RECIST.date.new, patient.ichor.dates)
      
      #Get patients' liquidCNA results.
      p.res <- liquidCNA_results[[p]]
      #order the samples so it matches RECIST ordering
      p.res.ord <- p.res[order(p.res$time, decreasing = F),]
      p.res.ord.rat <- as.numeric(p.res.ord$rat)
      p.res.ord.purity <- as.numeric(p.res.ord$purity_mean)
      
      #Get patients' ichor results.
      #p.ichor.res <- ichor.purities[[p]]
      
      #rat and purity for samples with RECIST
      rats <- p.res.ord.rat[samp.w.RECIST]
      purities <- p.res.ord.purity[samp.w.RECIST]
      #purities.ichor <- p.ichor.res[samp.w.RECIST]
      
      #append
      RECIST.rownum <- (rownames(patient.RECIST))
      if(length(RECIST.rownum) != 0){
        RECIST[RECIST.rownum,]$purity_mean <- purities
        RECIST[RECIST.rownum,]$rat <- rats
        #RECIST[RECIST.rownum,]$ichorCNA <- purities.ichor
      }
    }, error=function(e){cat("ERROR at:", p, "\n")})
  }
  ###############
  #Calculating time
  #First meta
  First.Meta <- RECIST$Date.Meta
  First.Meta.new <- structure(numeric(length(First.Meta)), class="Date")
  
  for(x in 1:length(First.Meta.new)){
    new.date <- as.Date(First.Meta[x], format = "%m/%d/%y")
    First.Meta.new[x] <- new.date
  }
  #Scan date
  Scan.date <- RECIST$Date
  time <- Scan.date - First.Meta.new
  RECIST <- cbind(RECIST,time)
  ##############################
  #new metastasis
  RECIST$New.Met <- NA
  
  ts.met <- function(row_num){
    p.RECIST <- RECIST[row_num,]
    if(p.RECIST$time[1] <= 10){
      m = FALSE
    } else {
      m <- "y" %in% unname(unlist(p.RECIST[,5:14]))
    }
    return(m)
  }
  
  ts.metastasis <- sapply(1:nrow(RECIST), function(x) ts.met(x))
  
  RECIST$New.Met <- ts.metastasis
  ##############################
  return(RECIST)
}

plot.estimation.w.RECIST <- function(p){ 
  p.id <- patient_ids[p]
  p.result <- liquidCNA_results[[p]]
  p.RECIST <- RECIST[which(RECIST$Patient_ID == p.id),]
  
  time.sample <- p.result$time
  sample.rat <- as.numeric(p.result$rat)
  sample.purity <- as.numeric(p.result$purity_mean)
  p.df.unord <- data.frame(time = time.sample, subclone = sample.rat, purity = sample.purity)
  p.df <- p.df.unord[order(p.df.unord$time),]
  
  ######
  #For which timeSamples do we have RECIST data?
  #subset of ichorCNA for this patient
  patient.ichor <- ichorCNA[ichorCNA$Patient_ID == p.id,]
  patient.ichor <- patient.ichor[order(patient.ichor$Sample_ID, decreasing = F),] #order by date
  patient.ichor.dates <- patient.ichor$Date
  
  #subset of RECIST for this patient
  p.RECIST.blood <- p.RECIST$Blood.draw
  p.RECIST.blood.new <- structure(numeric(length(p.RECIST.blood)), class="Date")
  
  p.RECIST.scan.date <- p.RECIST$Scan.date
  p.RECIST.scan.date.new <- structure(numeric(length(p.RECIST.scan.date)), class="Date")
  
  for(x in 1:length(p.RECIST.blood)){
    new.date <- as.Date(p.RECIST.blood[x], format = "%m/%d/%y")
    p.RECIST.blood.new[x] <- new.date
    
    new.scan <- as.Date(p.RECIST.scan.date[x], format = "%m/%d/%y")
    p.RECIST.scan.date.new[x] <- new.scan
  }
  
  #For which timeSamples do we have RECIST data?
  samp.w.RECIST <- match(p.RECIST.blood.new, patient.ichor.dates)
  
  #######
  #Append date p.df
  blood.date <-  patient.ichor.dates
  p.df <- cbind(p.df, blood.date)
  
  #Scan.date
  p.df$scan.date <- NA
  p.df$scan.date = as.Date(as.character(p.df$scan.date),format = "%m/%d/%y")
  p.df[samp.w.RECIST,]$scan.date <- p.RECIST.scan.date.new
  
  #Append boolean of RECIST to p.df
  RECIST.data <-  1:nrow(p.df) %in% samp.w.RECIST
  p.df <- cbind(p.df, RECIST.data)
  
  #Append progression to p.df
  progression <- rep(NA, nrow(p.df))
  progression[samp.w.RECIST] <- p.RECIST$Progression
  p.df <- cbind(p.df, progression)
  
  #Append Metastasis to p.df
  new.met <- rep(NA, nrow(p.df))
  new.met[samp.w.RECIST] <- p.RECIST$New.Met
  p.df <- cbind(p.df, new.met)
  
  #Get range of dates for x-axis to plot
  all.dates <- c(p.df$scan.date, p.df$blood.date)
  if(any(is.na(all.dates))){
    all.dates <- all.dates[-which(is.na(all.dates))]
  }
  xlim.min <- min(pmin(all.dates))
  xlim.max <- max(pmax(all.dates))
  
  
  #PLOT!
  par(mar = c(5, 5, 3.3, 4))
  #plot subclone estimation against date
  plot(subclone ~ blood.date, p.df, xaxt = "n", type = "b", pch = 19,
       ylab = "Subclonality & Purity estimates", xlab = "Dates",
       xlim = c(xlim.min, xlim.max),
       main = paste0("Patient " , p.id, " | LiquidCNA estimations with RECIST \n"),
       col = "indianred")
  axis(1, p.df$blood.date, format(p.df$blood.date, "%d %b %Y"), cex.axis = .7)
  #add purity estimation
  lines(purity ~ blood.date, p.df, xaxt = "n", type = "b", pch = 18,col = "steelblue")
  
  #add RECIST boolean
  p.df.w.RECIST <- p.df[p.df$RECIST.data,]
  # abline(v = p.df.w.RECIST$scan.date, lty= 2, col = ifelse(p.df.w.RECIST$progression == "YES", "yellowgreen", "lightpink3"))
  abline(v = p.df.w.RECIST$scan.date, lty= 2, col = ifelse(p.df.w.RECIST$new.met == TRUE, "yellowgreen", "lightpink3"))
  
  #legend
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 1.7, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  add_legend("topright", 
             legend=c("Subclonality", "Purity", "YES", "NO"), 
             pch=c(19,18,NA,NA), 
             lty=c(1,1,2,2),
             col=c("indianred", "steelblue", "yellowgreen", "lightpink3"),
             horiz=TRUE, 
             bty='n', 
             #x.intersp = 2,
             cex=0.8)
}

################################################################################

# #read outputs of QDNAseq
# load("../DATA/bam_by_patient.RData")
# p_num <- 80
# 
# seg.df <- read.delim(paste0("../DATA/2_QDNA_CNout/segment_", patient_ids[p_num], ".txt"))
# cn.df <- read.delim(paste0("../DATA/2_QDNA_CNout/raw_cn_", patient_ids[p_num], ".txt"))
# 
# seg.df <- seg.df[,-(1:4)]
# cn.df <- cn.df[,-(1:4)]
# 
# colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
# colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))
# 
# #optional
# #re-normalise raw and segmented CN values to be centred at CN=2
# reNorm <- centreSegs(seg.df)
# seg.df <- as.data.frame(t(t(seg.df)/reNorm)*2)
# cn.df <- as.data.frame(t(t(cn.df)/reNorm)*2)
# 
# #Plot the CN distribution of each sample to gain a quick overview
# explore_plot <- ggplot(reshape2::melt(seg.df), aes(x=value,y=..scaled.., colour=variable)) +
#   geom_density(adjust=1) +
#   theme_bw() + scale_x_continuous(limits=c(0.5, 5)) +
#   labs(x='Segment copy number',y='Density',colour='') + 
#   ggtitle(paste0('Patient ', patient_ids[p_num]))
# print(explore_plot)
# 
# #generate a dataframe of ensemble segments
# #contiguous sections of bins that are constant in ALL samples
# segchange <- sapply(1:(nrow(seg.df)-1),
#                     function(i) sum(seg.df[i,]!=seg.df[(i+1),])>0)
# seg.data <- data.frame(start=c(1,which(segchange)+1),
#                        end=c(which(segchange),length(segchange)+1))
# seg.data$length <- seg.data$end - seg.data$start+1
# 
# #Filter out small segments
# seg.sub <- subset(seg.data, length>12) #12 500kb bins
# cat('Total number of segments retained: ',nrow(seg.sub))
# 
# #Update segment values by fitting a normal distribution to bins in the segment
# seg.fit <- getNormalFitSegments(seg.sub, cn.df)
# 
# seg.cns <- seg.fit[[1]]
# names(seg.cns) <- names(seg.df)
# 
# seg.df.upd <- data.frame(matrix(NA, ncol=ncol(seg.df), nrow=nrow(seg.df)))
# names(seg.df.upd) <- names(seg.df); row.names(seg.df.upd) <- row.names(seg.df)
# 
# for(i in 1:nrow(seg.sub)){
#   seg.df.upd[seg.sub$start[i]:seg.sub$end[i],] <- seg.cns[i,]
# }
# 
# ##########################################
# ########### Purity Estimation ############
# ##########################################
# w = c(0.8,1,1,0.15,0.05) #weights of different CN states
# maxCN=8 #assumed maximum CN 
# adjVec = c(0.5,0.6,0.8,0.9,1,1.2,1.3,1.5,1.8,2) #smoothing kernel adjustments
# pVec = seq(0.05, 0.5, by=0.005) #range of purity values to be evaluated
# 
# #Estimated optimal purity values stored in pHat.df
# #purity vals that minimise mean and median
# pHat.df <- data.frame(matrix(vector(),ncol=ncol(seg.df.upd), nrow=2))
# names(pHat.df) <- names(seg.df.upd)
# row.names(pHat.df) <- c('mean','median')
# 
# #The plots show the error of the fit over the range of purity values 
# for(i in 1:ncol(seg.df.upd)){
#   #for each time sample...
#   
#   x <- na.omit(seg.df.upd[,i]) #get rid of NA rows (i.e., segments with NAs)
#   
#   pFits <- as.data.frame(sapply(adjVec,
#                                 function(a) sapply(pVec,
#                                                    function(p) evalPurityDensity(p,w,a,x,maxCN))))
#   #each row: the 1 purity value evaluated
#   #each column: 1 adjVec, a smoothing kernel adjustment
#   #each value: error for given purity/adjVec
#   
#   pFits$p <- pVec
#   mins <- c(pVec[which.min(apply(pFits[,1:length(adjVec)],1,mean))],
#             pVec[which.min(apply(pFits[,1:length(adjVec)],1,median))])
#   #which row (purity) has minimum mean error and minimum median error?
#   
#   plF <- ggplot(melt(pFits,id='p'), aes(x=p, y=value)) + geom_point(size=2, alpha=0.3)+
#     theme_bw() +
#     geom_vline(xintercept = mins, linetype=c('dashed','dotted'), colour=c('red','blue')) +
#     labs(title= paste0('Patient ', patient_ids[p_num], " | ", names(seg.df.upd)[i]))
#   print(plF)
#   pHat.df[,names(seg.df.upd)[i]] <- mins
# }
# 
# #Correct CN values by purity
# #removing CN contamination by normal
# #eq being used: CN_corr = 1/p * (CN-2) + 2
# pVec <- as.numeric(pHat.df[1,,drop=T])
# seg.df.corr <- as.data.frame(t(t(seg.df.upd-2)*1/pVec)+2)
# seg.cns.corr <- as.data.frame(t(t(seg.cns-2)*1/pVec)+2)
# 
# #remove low purity samples if number of number of samples is many
# if( length(pVec) > 6 ){
#   purity_threshold <- 0.1
#   abovePurTh <- pVec >= purity_threshold
#   seg.df.corr <- seg.df.corr[,abovePurTh]
#   seg.cns.corr <- seg.cns.corr[,abovePurTh]
# }
# 
# #Plot purity-corrected segment distribution
# purity_plot <- ggplot(melt(as.data.frame(seg.df.corr)), aes(x=value, colour=variable)) +
#   geom_density(adjust=1) + theme_bw() +
#   scale_x_continuous(limits=c(0,8)) + geom_vline(xintercept = 1:6) +
#   labs(x='Purity-corrected segment CN',y='Density',colour='') + 
#   ggtitle(paste0('Patient ', patient_ids[p_num]))
# print(purity_plot)
# 
# #Designate reference sample
# #& calculate dCN
# baseSample <- names(seg.df.corr)[1]
# seg.dcn <- seg.cns.corr - seg.cns.corr[,baseSample]
# seg.dcn.nonBase <- seg.dcn %>% select(-one_of(baseSample))
# 
# #Choose samples to investigate (all non-base samples by default) 
# colToUse <- names(seg.dcn.nonBase)
# #generate all possible permutations
# nCol <- length(colToUse)
# seg.dcn.toOrder <- seg.dcn.nonBase[,colToUse]
# ordVec <- permutations(nCol,nCol,colToUse)
# 
# #Permutation is too long
# if( nrow(ordVec) > 1000 ){
#   
#   chrono <- which(apply(ordVec, 1, function(x) identical(x, colToUse)))
#   nonchrono <- (1:nrow(ordVec))[-chrono]
#   rows <- sample(nonchrono, 1000, replace = F)
#   subOrdVec <- ordVec[rows,]
#   subOrdVec <- rbind(subOrdVec, ordVec[chrono,])
#   ordVec <- subOrdVec
#   
# }
# 
# #Initial filtering to find best sample order 
# #Set parameters to be used when evaluating sample order and segment monotony
# filterMethod <- 'sd' #filter by saying... if sd < theta, segment = clonal
# cutOffVec <- seq(0.025,0.5,by=0.005) #thetas to try for initial-filtering
# epsilon <- 0.05 #error margin for monotony
# 
# #Diagnostic plot
#   #at different thetas shows:
#     #the proportion of subclonal segments (of non-clonal segments)
#     #total number of non-clonal (light blue)
#     #and subclonal (dark blue) segments
# fitInfo <- list()
# for (cutOff in cutOffVec){
#   seg.dcn.Eval <- filterSegmentRatios(seg.dcn.toOrder, cutOff, filterMethod, 0)
#   #Filter relative segment values (compared to baseline) to discard clonal 
#     #(nonchanging) segments that stay relatively constant
#     #Returns: a filtered set of a segments that have values exceeding the threshold
#   
#   best <- findBestOrder(seg.dcn.Eval, ordVec, epsilon, nCol, 0)
#   #Input: segments where clonal segments are initially filtered
#   #Returns: 
#     #$cons: row-names of all segments considered, 
#     #max: maximum fit gained, 
#     #$ord: order with the maximum fit, 
#     #$segs:segments monotone according to max fit
#   fitInfo[[as.character(cutOff)]] <- best
# }
# 
# fitInfo.df <- data.frame(cutOff = cutOffVec,
#                          maxFit = sapply(fitInfo, function(x) x$max),
#                          segsAboveCut = sapply(fitInfo, function(x) length(x$cons)),
#                          segsInOrder = sapply(fitInfo, function(x) max(sapply(x$segs, length))))
# p1 <- ggplot(fitInfo.df, aes(x=cutOff, y=maxFit)) + geom_line(size=2,colour='darkgreen') + theme_bw() +
#   labs(x='Clonal cut-off',y='Subclonal proportion (of non-clonal segments)') + 
#   ggtitle(paste0('Patient ', patient_ids[p_num]))
# p2 <- ggplot(fitInfo.df, aes(x=cutOff, y=segsAboveCut)) + geom_line(size=2,colour='deepskyblue3') +
#   geom_line(aes(y=segsInOrder), size=2,colour='dodgerblue4') + theme_bw() +
#   labs(x='Clonal cut-off',y='Total number of non-clonal and subclonal segments') + 
#   ggtitle(paste0('Patient ', patient_ids[p_num]))
# grid.arrange(p1,p2,nrow=1)
# #p1: proportion of subclonal segments
# #p2: total number of non-clonal (light blue) and subclonal (dark blue) segments
# 
# #A typical cut-off ~0.1 is adequate
# #too low or too high cut-off will select for noise
# #ideal cut-off should have a reasonable number of subclonal segments (5-20) 
# #while maximises the subclonal proportion
# 
# #optimal cut-off for wanted  number of segment
# minSegmentNumber <- 9
# recommendCutOff <- getCutOffAuto(fitInfo.df, minSegmentNumber)
# cat('Recommended cut-off is: ', recommendCutOff)
# 
# #plot DeltaCN for chosen cut off + ordering + segment classification
# cutOff <- recommendCutOff
# fit <- fitInfo[[as.character(cutOff)]]
# 
# for (ind in 1:nrow(fit$ord)){
#   print(ind)
#   seg.plot <- seg.dcn[,rev(c(fit$ord[ind,],baseSample))]
#   seg.plot$id <- row.names(seg.plot)
#   seg.plot$filtered <- seg.plot$id %in% fit$cons
#   seg.plot$order <- seg.plot$id %in% fit$segs[[ind]]
#   
#   p <- ggplot(melt(seg.plot, id=c('id','filtered','order')), aes(x=variable, y=value, group=id, colour=paste0(filtered,' | ',order))) +
#     geom_line(size=1.2, alpha=0.75) + theme_bw() +
#     scale_colour_manual(values=c('grey70','#487a8b','firebrick3'), labels=c('clonal','unstable','subclonal')) +
#     labs(x='',y='Change in segment CN (DeltaCN)',colour='Segment type') + 
#     ggtitle(paste0('Patient ', patient_ids[p_num]))
#   print(p)
# }
# 
# #if more than one sample order exists
# ordInd <- 1
# seg.dcn.toUse <- seg.dcn[fit$segs[[ordInd]],rev(fit$ord[ordInd,])]
# 
# 
# ##########################################
# ######## Derive  subclonal-ratio #########
# ##########################################
# 
# ###
# #RELATIVE subclonal-ratio estimates
# #via segment-by-segment comparison to the highest subclonal-ratio sample
# topSample <- names(seg.dcn.toUse)[ncol(seg.dcn.toUse)]
# toEstimate <- setdiff(names(seg.dcn.toUse), c(topSample,baseSample))
# 
# final.ratios <- estimateRSegmentRatio(seg.dcn.toUse, toEstimate, topSample, 1)
# 
# final.medians <- aggregate(final.ratios$value, 
#                            by=list(final.ratios$variable), median)
# 
# #plot estimated relative ratios
# final.medians$xpos <- 1:nrow(final.medians)
# 
# ggplot(final.ratios, aes(y=value, x=variable, weight=1)) + geom_violin(fill='firebrick3',alpha=0.2) +
#   theme_bw() + scale_y_continuous(limits=c(0,1.2)) +
#   geom_jitter(width=0.1, height=0, colour='firebrick3', size=2, alpha=0.6) +
#   #geom_point(data=true.df, aes(x=V3, y=V2), colour='firebrick3', size=3) +
#   geom_segment(data=final.medians,aes(x=xpos-0.15,xend=xpos+0.15,y=x,yend=x), size=1.5) +
#   labs(x='',y=paste0('Subclonal-ratio compared to ',topSample)) + 
#   ggtitle(paste0('Patient ', patient_ids[p_num]))
# 
# #final results table of median relative ratio values
# final.results <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
# names(final.results) <- c('time','relratio')
# final.results$time <- as.character(final.results$time)
# final.results[nrow(final.results)+1,] <- c(topSample, 1)
# 
# ###
# #Absolute subclonal-ratio estimates
# #by fitting a Gaussian mixture model of constrained means to each sample
# #shared mean parameter defines the absolute ratio
# 
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
#   labs(x='Sample', y='Subclonal-ratio') +
#   ggtitle(paste0('Patient ', patient_ids[p_num]))
# 
# ##########################################
# ############# Final results ##############
# ##########################################
# final.results[nrow(final.results)+1,] <- c(baseSample, 0,0,0)
# final.results$purity_mean <- as.numeric(pHat.df[1,match(final.results$time, 
#                                                         names(pHat.df))])
# final.results$purity_median <- as.numeric(pHat.df[2,match(final.results$time, 
#                                                           names(pHat.df))])
# 
# cat('Final result table:\n'); print(final.results)


