###############################################################################
#load packages
library(future)
library(QDNAseq)
library(QDNAseq.hg38)
###############################################################################
all_readCounts <- list.files("./2_QDNA_readCounts/")

load("bam_by_patient.RData")

qdna_analysis <- function(p){
  #patient number
  p_num <- p
  print(paste0("Starting analysis number: ", p_num))
  
  #load patient readCounts (from 2a_qdna_process_bam.R)
  load(paste0("./2_QDNA_readCounts/",all_readCounts[p_num]))
  
  #Plot1
  #Plot a raw copy number profile (read counts across the genome)
  plot(readCounts, logTransform=FALSE, ylim=c(-50, 2000))
  #Highlight filtered bins
  highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
  
  #filter
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE,
                                     mappability = 75)
  #Plot2
  #plot median read counts as a function of GC content and mappability
  isobarPlot(readCountsFiltered)
  
  #Estimate the correction for GC content and mappability
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  
  #Plot 3
  #a plot for the relationship between the observed standard deviation 
  #in the data and its read depth
  #Samples with low-quality DNA will be noisier than expected 
  #and appear further above the line than good-quality samples.
  noisePlot(readCountsFiltered)
  
  #apply the correction for GC content and mappability
  copyNumbers <- correctBins(readCountsFiltered) 
  copyNumbersNormalized <- normalizeBins(copyNumbers) #normalise
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized) #smooth outliers
  print("CN filtered")
  
  #Plot 4
  #Plot the copy number profile
  plot(copyNumbersSmooth)
  
  
  ##################
  ###SEGMENTATION###
  ##################
  #Segmentation with the CBS algorithm from DNAcopy
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="log2")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  print("Segmentation with DNAcopy")
  
  #Plot 5
  plot(copyNumbersSegmented)
  
  
  #call aberrations
  #EM algorithm
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  plot(copyNumbersCalled)
  
  #COPY NUMBER OF EACH BIN
  exportBins(copyNumbersCalled, type = "copynumber", 
             paste0("raw_cn_",patient_ids[p_num],".txt"), logTransform = FALSE)
  #COPY NUMBER OF SEGMENTS
  exportBins(copyNumbersCalled, type = "segments", 
             paste0("segment_",patient_ids[p_num],".txt"), 
             logTransform = FALSE)
}


p_vec <- 1:length(bam_by_patient)
sapply(p_vec, function(x) qdna_analysis(x))


# #Plot1
# #Plot a raw copy number profile (read counts across the genome)
# plot(readCounts, logTransform=FALSE, ylim=c(-50, 2000))
# #Highlight filtered bins
# highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
# 
# #filter
# readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE,
#                                    mappability = 75)
# #Plot2
# #plot median read counts as a function of GC content and mappability
# isobarPlot(readCountsFiltered)
# 
# #Estimate the correction for GC content and mappability
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# 
# #Plot 3
# #a plot for the relationship between the observed standard deviation 
# #in the data and its read depth
# #Samples with low-quality DNA will be noisier than expected 
# #and appear further above the line than good-quality samples.
# noisePlot(readCountsFiltered)
# 
# #apply the correction for GC content and mappability
# copyNumbers <- correctBins(readCountsFiltered) 
# copyNumbersNormalized <- normalizeBins(copyNumbers) #normalise
# copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized) #smooth outliers
# print("CN filtered")
# 
# #Plot 4
# #Plot the copy number profile
# plot(copyNumbersSmooth)
# 
# 
# ##################
# ###SEGMENTATION###
# ##################
# #Segmentation with the CBS algorithm from DNAcopy
# copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="log2")
# copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
# print("Segmentation with DNAcopy")
# 
# #Plot 5
# plot(copyNumbersSegmented)
# 
# 
# #call aberrations
# #EM algorithm
# copyNumbersCalled <- callBins(copyNumbersSegmented)
# plot(copyNumbersCalled)
# 
# #COPY NUMBER OF EACH BIN
# exportBins(copyNumbersCalled, type = "copynumber", 
#            "raw_cn.txt", logTransform = FALSE)
# #COPY NUMBER OF SEGMENTS
# exportBins(copyNumbersCalled, type = "segments", "segment.txt", 
#            logTransform = FALSE)


