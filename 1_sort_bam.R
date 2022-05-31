##############################################################################

load("antwerpData.RData")
bam_files <- read.csv("bam_file_names.txt")
bam_files <- bam_files[-grep(pattern = ".bai", x = bam_files[,1]),]

patient_ids <- unique(ichorCNA$Patient_ID)
getPatientBam <- function(p){
  patient <- ichorCNA[ichorCNA$Patient_ID == patient_ids[p],]
  patient <- patient[order(patient$Sample_ID, decreasing = F),] #order by date
  
  patient$SLX
  patient_barcode <- strsplit(patient$Barcode, "-")
  barcodes <- NULL
  for(i in 1:length(patient_barcode)){
    split <- strsplit(patient$Barcode, "-")[[i]]
    barcode <- paste0(split[1], "_", split[2])
    barcodes <- c(barcodes, barcode)
  }
  
  identifier <- paste0(patient$SLX, ".", barcodes)
  patient_bam <- sapply(1:length(identifier), 
                        function(x) bam_files[grep(pattern = identifier[x], x = bam_files)])
  return(patient_bam)
}

bam_by_patient <- sapply(1:length(patient_ids), function(x) getPatientBam(x))
(bam_by_patient[[2]])
##############################################################################

# library(DNAcopy)
# data(coriell)
# 
# 
# length(bam_files) == nrow(ichorCNA)
# 
# ichorCNA[ichorCNA$Barcode == "D703tp-D508tp",]
# identifier[7]
# 
# raw_cn <- read.csv("test_samples/raw_cn.txt", sep = "\t")
# segment <- read.csv("test_samples/segment.txt", sep = "\t")
# 
# 
# head(raw_cn)
# head(segment)
# 
# seg.df <- segment[,-(1:4)]
# colnames(seg.df) <- paste0("Sample", 1:8)
# cn.df <- raw_cn[,-(1:4)]
# colnames(cn.df) <- paste0("Sample", 1:8)
# 
# length(unique(seg.df[,5]))
# length(unique(cn.df[,5]))
# 
# plot(density(seg.df[,5]))