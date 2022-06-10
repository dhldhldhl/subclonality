x <- na.omit(seg.df.upd[,i]) #get rid of NA rows (i.e., segments with NAs) # x is a a vector
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

####################
#eval_purity_density
function(p,w,a,x,maxCN=8,nu=7,nd=4){
  #p is purity value being evaluated
  #w is weight of different CN states
  #a is each element of adjVec
  #x are CN for all non-NA segments for one time sample
  
  seg.dist <- density(x[x>1.25 & x < 4.25], adjust=a)
  #get rid of extreme CNs
  #then kernel density est. with adjust smoothing bandwidth by value of adjVec
  
  z <- seg.dist$x[findpeaks(seg.dist$y,nups=nu, ndowns=nd)[,2]]
  #nups: minimum number of increasing steps before a peak is reached #=7
  #ndowns: minimum number of decreasing steps after the peak #=4
  #second column gives the position/index where the maximum is reached
  #thus, z is the peaks of CN values found by peak finding algorithm
  
  zhat2 <- sapply(1:length(w), function(i) min( ( ( 2 + p*(i-2) ) - z ) ^2 ) )
  #(2 + p*(i-2)) is the expected peak for CN i #page 18 Eq 5.
    #expectation under *assumption* that:
      #C(A) take only integer values
        #thereby, peaks being at regular intervals of p
  #(( 2 + p*(i-2) ) -z ) ^2) = sqrd dist. of expected and observed peaks
    #for each i, z is a vector ==> output to be vector
    #each output element, is squared distance to each of the observed peaks
    #minimum means, we only extract the error to the closest observed peak
  #thus, zhat2 is a vector of squared distances 
    #between each observed peaks with the nearest expected peaks
  
  if (min(((2+p*(maxCN-2))-z)^2) < p){
    zhat2 = zhat2+1
  }
  
  return(sqrt(sum(zhat2*w)))
  #the summed squared distance of observed and expected peaks
    #zhat2 is the error for each peak (i.e. CN values)
    #high ploidy is measured less accurately thus is given lower weight
}

###############################################################################
#evaluating segment subclonality

#find_best_order()
function(seg.ratios.Eval, ordVec, threshold, nCol, base=1){
  #findBestOrder(seg.dcn.Eval, ordVec, epsilon, nCol, 0)
  #so threhold here = epsilon
  if (nrow(seg.ratios.Eval)==0){
    return(list(cons=row.names(seg.ratios.Eval),
                max=NA,
                ord=NA,
                segs=list(NULL)))
  }
  if (nCol<2){
    return(list(cons=row.names(seg.ratios.Eval),
                max=1,
                ord=ordVec[1,,drop=F],
                segs=list(row.names(seg.ratios.Eval))))
  }
  fitTotal <- c()
  fitSegments <- list()
  for (ind in 1:nrow(ordVec)){
    ord <- ordVec[ind,]
    #for each permutation of time point samples...
    
    fitCounts <- apply(seg.ratios.Eval[,ord],1,function(z) evalOrder(z,threshold,nCol, base))
      #for each segment...
      #evaluate if change in CN across order is monotone or not
    
    fitTotal[[ind]] <- sum(fitCounts)/(nrow(seg.ratios.Eval))
    fitSegments[[ind]] <- row.names(seg.ratios.Eval)[fitCounts==1]
  }
  
  # return the order with the best fit and segments that match it
  inds <- which(fitTotal==max(unlist(fitTotal)))
  return(list(cons=row.names(seg.ratios.Eval),
              max=max(unlist(fitTotal)),
              ord=ordVec[inds,,drop=F],
              segs=fitSegments[inds]))
}

# Evaluate a monotony of a sample order within a segment
# Input:
# - x = vector of segment CN values, ordered
# - th = error-threshold to allow for uncertain measurements of equal values <-- EPSILON
# - nCol = total number of samples observed and evaluated
# Returns: 1 or 0 depending if the segment is monotone ordered with threshold th
evalOrder <- function(x, th, nCol, base=1){
  diff <- x - c(x[2:length(x)],base)
    #diff in segment CN between one time point and the next (in the order)
  if ( sum(x > base)>nCol/2 ){
    #Eq 11a
    orderFit <- sum(diff> -1*th)
    
  } else {
    #Eq 11b
    orderFit <- sum(diff< th)
    
    }
  return(sum(orderFit==nCol))
}


################################################################################
# Find the best (highest proportion of subclonal vs random segments) cutoff for clonal segments
# given a required minimal number of subclonal segments to recover (if known in advance)
# Input:
# - fitInfo.df = output of fitting at a range of clonal cutoff values
# - minNum = minimum number of subclonal (monotone) segments required to be returned
# Returns: cutoff value meeting requirement and maximising the fit
getCutOffAuto <- function(fitInfo.df, minNum){
  x <- subset(fitInfo.df, segsInOrder > minNum)
  if(nrow(x)>0){
    c <- x$cutOff[which.max(x$maxFit)]
  }else{
    c <- fitInfo.df$cutOff[which.max(fitInfo.df$maxFit)]
  }
  return(c)
}



