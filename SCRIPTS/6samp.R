timely_patients <- which(sample_num_patient > 6)
timely_bam_num <- sample_num_patient[timely_patients]

nbatch1 <- ceiling(timely_bam_num/2)
nbatch2 <- timely_bam_num - nbatch1


res1 <- vector(mode = "list", length = length(timely_patients))
res2 <- vector(mode = "list", length = length(timely_patients))

for(x in 1:length(timely_patients)){
  tryCatch({
    cat("\n Starting", x, "\n")
    res1[[x]] <- run_liquidCNA(timely_patients[x], 
                               timely = TRUE,
                               batch1 = TRUE,
                               nbatch = nbatch1[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

for(x in 1:length(timely_patients)){
  tryCatch({
    cat("\n Starting", x, "\n")
    res2[[x]] <- run_liquidCNA(timely_patients[x], 
                               timely = TRUE,
                               batch1 = FALSE,
                               nbatch = nbatch2[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}



##############################################################################
estimateRGaussianFit(seg.dcn.toUse, samp, final.results, 3)
seg.rel.toUse = seg.dcn.toUse
topSample = "Sample4"
final.medians = final.results
nStates = 2 #what is nStates
gaussSigma=0.1 #why is gaussSigma 3
lambdaOverwrite=0.2
verbose=F

estimateRGaussianFit <- function(seg.rel.toUse, topSample, 
                                 final.medians, gaussSigma=0.1, nStates=3, 
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
  mixfit <- normalmixEM(as.numeric(seg.top),
                        lambda=lambdaOverwrite,
                        mean.constr = c("2a","a","-a","-2a","-3a"),
                        mu=c(2,1,-1,-2,-3)*start,sigma=gaussSigma)
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
