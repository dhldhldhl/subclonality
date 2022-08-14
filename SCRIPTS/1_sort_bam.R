##############################################################################
load("../DATA/antwerpData.RData") #This imports RECIST and ichorCNA
bam_files <- read.csv("../DATA/bam_file_names.txt") #This imports all the file names in the HPC data folder
bam_files <- bam_files[-grep(pattern = ".bai", x = bam_files[,1]),] #removes files of .bam.bai, leaving only .bam

patient_ids <- unique(ichorCNA$Patient_ID)
patient_ids <- sort(patient_ids)

#a function to link each patient w/ their bam files
getPatientBam <- function(patient_id){
  patient <- ichorCNA[ichorCNA$Patient_ID == patient_id,] #get ichorCNA rows for one patient
  patient <- patient[order(patient$Date, decreasing = F),] #order by date 
  
  #changing barcode string so it matches bam_file_names.txt; 
  #replaces "-" to "_"
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

bam_by_patient <- sapply(patient_ids, function(x) getPatientBam(x))
names(bam_by_patient) <- paste0("patient_", patient_ids)
save(bam_by_patient, patient_ids, file = "../DATA/bam_by_patient.RData")


