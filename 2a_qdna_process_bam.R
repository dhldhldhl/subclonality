#setwd
wd <- "/rds/project/rds-csoP2nj6Y6Y/ctDNA/dl"
setwd(wd)
bam_dir <- "/rds/project/rds-csoP2nj6Y6Y/ctDNA/data/Antwerp/bam_files/mapped_hg38/"
###############################################################################

#load packages
library(future)
library(QDNAseq)
library(QDNAseq.hg38)

###############################################################################

#parallelise
future::plan("multisession")

###############################################################################

#500kb bin annotation of genome
bins <- getBinAnnotations(binSize=500, genome="hg38")
print("Bin annotations obtained")

#load bam_by_patient from 1_sort_bam.R
load("bam_by_patient.RData")

process_bam <- function(p){
  #patient number
  p_num <- p
  
  #Processing BAM files
  setwd(bam_dir)
  test_bams <- bam_by_patient[[p_num]]
  readCounts <- binReadCounts(bins, bamfiles=test_bams)
  print("Bam files processed")
  
  #save readCounts as RData file
  setwd(wd)
  filename <- paste0("2_QDNA_readCounts/readCounts_patient", 
                     patient_ids[p_num], ".RData")
  save(readCounts, file = filename)
  print("readCounts saved")
}

p_vec <- 1:length(bam_by_patient)
sapply(p_vec, function(x) process_bam(x))

# #patient number
# p_num <- 1
# 
# #Processing BAM files
# setwd(bam_dir)
# test_bams <- bam_by_patient[[p_num]]
# readCounts <- binReadCounts(bins, bamfiles=test_bams)
# print("Bam files processed")
# 
# #save readCounts as RData file
# setwd(wd)
# filename <- paste0("2_QDNA_readCounts/readCounts_patient", 
#                    patient_ids[p_num], ".RData")
# save(readCounts, file = filename)
# print("readCounts saved")