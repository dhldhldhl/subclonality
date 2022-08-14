estimateRSegmentRatio(seg.dcn.toUse, toEstimate, topSample, 1)
seg.ratios = seg.dcn.toUse
topSamples = topSample
w = 1
samp = "Sample2"
top = topSamples
final.ratios <- data.frame(matrix(vector()))
  
  for (samp in toEstimate){
    tmp <- c()
    for (top in topSamples){
      
      ##tmp <- c(tmp,(seg.ratios[,samp])/(seg.ratios[,top]))
      
      
      tmp <- c(tmp,(seg.ratios[,samp])/nrow(seg.ratios))
      
      #simple ratio of purity corrected CNs of segment
      # = (non-top-sample / top-sample)
      
    }
    final.ratios <- rbind(final.ratios, data.frame(value=tmp, variable=samp, weight=w))
    
  }


#mean
mean(seg.dcn.toUse[,1])
mean(seg.dcn.toUse[,2])
mean(seg.dcn.toUse[,1]) / mean(seg.dcn.toUse[,2])

#median
median(seg.dcn.toUse[,1]/nrow(seg.dcn.toUse))
median(seg.dcn.toUse[,2]/nrow(seg.dcn.toUse))
median(seg.dcn.toUse[,1]/nrow(seg.dcn.toUse)) / median(seg.dcn.toUse[,2]/nrow(seg.dcn.toUse))

#I think in the package CNs of segments are ratioed against topSample as it allows
#CN to be evaluated segment by segment,
#thus dCN is calculated first then averaged
#rather than averaged then dCN,
#I am guessing that this will result less bias for extreme (large/small) CN of some segments.
#which must be why they are doing medians > means as well.
#can look at this further by looking at synthetic data?


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
  labs(x='',y=paste0('Subclonal-ratio compared to ',topSample)) + 
  ggtitle(paste0('Patient ', patient_ids[p_num]))

#final results table of median relative ratio values
final.results <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
names(final.results) <- c('time','relratio')
final.results$time <- as.character(final.results$time)
final.results[nrow(final.results)+1,] <- c(topSample, 1)

###
# Absolute subclonal-ratio estimates
#by fitting a Gaussian mixture model of constrained means to each sample
#shared mean parameter defines the absolute ratio

final.results$rat <- NA
final.results$rat_sd <- NA
for (samp in final.results$time){
  final.results <- estimateRGaussianFit(seg.dcn.toUse, samp, final.results, 3) #what is this default of three? gaussSigma?
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
final.results[nrow(final.results)+1,] <- c(baseSample, 0,0,0)
final.results$purity_mean <- as.character(pHat.df[1,match(final.results$time, 
                                                          names(pHat.df))])
final.results$purity_median <- as.character(pHat.df[2,match(final.results$time, 
                                                            names(pHat.df))])
return(final.results)



################################################################################
seg.rel.toUse = seg.dcn.toUse
topSample = "Sample3"
final.medians = final.results
gaussSigma=3
nStates=3
lambdaOverwrite=0.2
verbose=F

estimateRGaussianFit <- function(seg.rel.toUse, topSample, 
                                 final.medians, gaussSigma=0.1, nStates=2, 
                                 lambdaOverwrite=0.2, verbose=F){
  seg.top <- seg.rel.toUse[,topSample]
  start <- mean(seg.top[seg.top<0])
  #average of deltaCNs of top sample that are negative
  
  if (is.nan(start)){
    # estimate the initial value for fitting from segments with DCN between 0 and 1
    start <- -1* mean(seg.top[seg.top>0 & seg.top<1])
  }
  if (is.nan(start)){ 
    # if no starting point can be found, start from a default of 25%
    start <- -0.25
  }
  
  # By default the following DCN values are used: -2, -1, 1, 2, 3
  #function for fitting of mixture distributions
  mixfit <- normalmixEM(as.numeric(seg.top),
                        lambda=lambdaOverwrite,
                        mean.constr = c("2a","a","-a","-2a","-3a"),
                        mu=c(2,1,-1,-2,-3)*start,
                        sigma=gaussSigma)
  plot(mixfit,which = 2)
  
  mixfit2 <- normalmixEM(as.numeric(seg.top), k = 5)
  plot(mixfit2, which = 2)
  
  mixfit_res <- t(data.frame(lambda = mixfit$lambda, mu = mixfit$mu, mixfit$sigma))
  
  mixfit_res <- t(data.frame(lambda = mixfit$lambda, mu = mixfit$mu, mixfit$sigma))
  
  # only return the obtained value if converged in less than 800 iterations
  if (length(mixfit$all.loglik)<800){
    final.medians$rat[final.medians$time==topSample] <- -1*(mixfit$mu[2])
    #sigma depends on the number of segments used and that depends on how many states they were distributed across
    final.medians$rat_sd[final.medians$time==topSample] <- (mixfit$sigma[1])/sqrt(nrow(seg.rel.toUse)/nStates)
  }
  if (verbose){
    print(mixfit$lambda)
  }
  return(final.medians)
}

#Questions
#What is start? why mean of just negative dCNs?

mixfit2 <- normalmixEM(seg.top)





