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

#eval_purity_density

function(p,w,a,x,maxCN=8,nu=7,nd=4){
  #p is purity value being evaluated
  #w is weight of different CN states
  #a is each element of adjVec
  
  seg.dist <- density(x[x>1.25 & x < 4.25], adjust=a)
  #get rid of extreme CNs
  #then kernel density est. with adjust smoothing bandwidth by value of adjVec
  
  z <- seg.dist$x[findpeaks(seg.dist$y,nups=nu, ndowns=nd)[,2]]
  #nups: minimum number of increasing steps before a peak is reached
  #ndowns: minimum number of decreasing steps after the peak
  #second column gives the position/index where the maximum is reached
  #thus, z is the peaks of CN values found by peak finding algorithm
  
  zhat2 <- sapply(1:length(w), function(i) min( ( ( 2 + p*(i-2) ) - z ) ^2 ) )
  #(2 + p*(i-2)) is the expected peak for CN i
    #expectation under assumption that:
      #C(A) take only integer values
        #thereby, peaks being at regular intervals of p
  #(( 2 + p*(i-2) ) -z ) ^2) = sqrd dist. of expected and observed peaks
    #for each i, z is a vector ==> output to be vector
    #each output element, is squared distance to each of the observed peaks
    #minimum means, we only extract the error to the closest observed peak
  #thus, zhat2 is a vector of squared distanceas 
    #between each observed peaks with the nearest expected peaks
  
  if (min(((2+p*(maxCN-2))-z)^2) < p){
    zhat2 = zhat2+1
  }
  
  
  return(sqrt(sum(zhat2*w)))
  #the summed squared distance of observed and expected peaks
    #zhat2 is the error for each peak (i.e. CN values)
    #high ploidy is measured less accurately thus is given lower weight
}