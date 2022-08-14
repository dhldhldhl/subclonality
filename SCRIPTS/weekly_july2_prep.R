##############################################################################
#Categorise patients

#Patients with only two samples
sample_num_patient <- sapply(1:80, function(x) length(bam_by_patient[[x]]))
two_samp_patients <- which(sample_num_patient == 2)

#Patients with over 6 samples
timely_patients <- which(sample_num_patient > 6)

go_patients <- 1:80
go_patients <- go_patients[!(go_patients %in% two_samp_patients)]
go_patients <- go_patients[!(go_patients %in% timely_patients)]

##############################################################################
#800
#initialise results vector
go_results_800 <- vector(mode = "list", length = length(go_patients))

#Run for go_patients (i.e., non 2 or timely sample patients)
for(patient_x in 1:length(go_patients)){
  tryCatch({
    cat("Starting patient :", patient_x, "\n")
    go_results_800[[patient_x]] <- run_liquidCNA(go_patients[patient_x])
  }, error=function(e){cat("ERROR")})
}

names(go_results_800) <- paste0("patient_",patient_ids[go_patients])

##########################
#1500
#initialise results vector
go_results_1500 <- vector(mode = "list", length = length(go_patients))

#Run for go_patients (i.e., non 2 or timely sample patients)
for(patient_x in 1:length(go_patients)){
  tryCatch({
    cat("Starting patient :", patient_x, "\n")
    go_results_1500[[patient_x]] <- run_liquidCNA(go_patients[patient_x])
  }, error=function(e){cat("ERROR")})
}

names(go_results_1500) <- paste0("patient_",patient_ids[go_patients])

##########################
#3000
go_results_3000 <- vector(mode = "list", length = length(go_patients))

#Run for go_patients (i.e., non 2 or timely sample patients)
for(patient_x in 1:length(go_patients)){
  tryCatch({
    cat("Starting patient :", patient_x, "\n")
    go_results_3000[[patient_x]] <- run_liquidCNA(go_patients[patient_x])
  }, error=function(e){cat("ERROR")})
}

names(go_results_3000) <- paste0("patient_",patient_ids[go_patients])

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# 7. stitch timely

#number of BAM files our timely patients have
timely_bam_num <- sample_num_patient[timely_patients]

#number of bam files to batch each timely patient by:
nbatch1 <- ceiling(timely_bam_num/2)
nbatch2 <- timely_bam_num - nbatch1

#option 1) make topSample of batch1 baseSample of batch2
unlist(sapply(1:length(timely_res1), function(x) tail(timely_res1[[x]]$time, 2)[1]))

#initialise results vector for option 1
timely_res2.opt1 <- vector(mode = "list", length = length(timely_patients))

###############################################################################
p_num = timely_patients[1]
seg.df <- read.delim(paste0("../DATA/2_QDNA_CNout/segment_", patient_ids[p_num], ".txt"))
cn.df <- read.delim(paste0("../DATA/2_QDNA_CNout/raw_cn_", patient_ids[p_num], ".txt"))

seg.df <- seg.df[,-(1:4)]
cn.df <- cn.df[,-(1:4)]

colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))

###############################################################################
timely.n = 1
timely = TRUE
batch1 = FALSE
nbatch = nbatch2[x]

####################
#option 1: make topSample of batch1 baseSample of batch2
timely_res2.opt1 <- vector(mode = "list", length = length(timely_patients))

for(x in 1:length(timely_patients)){
  tryCatch({
    cat("\n Starting", x, "\n")
    batch1.top <- tail(timely_res1[[x]]$time, 2)[1]
    base.id <- as.numeric(strsplit(batch1.top, "e")[[1]][2])
    
    timely_res2.opt1[[x]] <- run_liquidCNA(timely_patients[x], 
                                           timely = TRUE,
                                           batch1 = FALSE,
                                           nbatch = nbatch2[x])
    
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

names(timely_res2.opt1) <- paste0("patient_",patient_ids[timely_patients],"_batch2")



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


####################
#option 2:  time sample with highest ratio as baseSamp for batch2
timely_res2.opt2 <- vector(mode = "list", length = length(timely_patients))

for(x in 1:length(timely_patients)){
  tryCatch({
    cat("\n Starting", x, "\n")
    
    batch1.rats <- as.numeric(timely_res1[[x]]$rat)
    batch1.max.rat <- timely_res1[[x]]$time[which(batch1.rats == max(batch1.rats))]
    base.id <- as.numeric(strsplit(batch1.max.rat, "e")[[1]][2])
    
    timely_res2.opt2[[x]] <- run_liquidCNA(timely_patients[x], 
                                           timely = TRUE,
                                           batch1 = FALSE,
                                           nbatch = nbatch2[x])
    
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

names(timely_res2.opt2) <- paste0("patient_",patient_ids[timely_patients],"_batch2")



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

####################
#option3: Sample1 is baseSample for both batch
timely_res2.opt3 <- vector(mode = "list", length = length(timely_patients))

for(x in 1:length(timely_patients)){
  tryCatch({
    cat("\n Starting", x, "\n")
    timely_res2.opt3[[x]] <- run_liquidCNA(timely_patients[x], 
                                           timely = TRUE,
                                           batch1 = FALSE,
                                           nbatch = nbatch2[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

names(timely_res2.opt3) <- paste0("patient_",patient_ids[timely_patients],"_batch2")



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


##############################################################################
#8. APPEND TO RECIST
RECIST$rat <- NA
RECIST$purity_mean <- NA

#get one patient
p.id <- patient_ids[1]
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

#Get patients' liquidCNA results... just need to figure this out
p.res <- go_results$patient_339
#order the samples so it matches RECIST ordering
p.res.ord <- p.res[order(p.res$time, decreasing = F),]
p.res.ord.rat <- as.numeric(p.res.ord$rat)
p.res.ord.purity <- as.numeric(p.res.ord$purity_mean)

#rat and purity for samples with RECIST
rats <- p.res.ord.rat[samp.w.RECIST]
purities <- p.res.ord.purity[samp.w.RECIST]

#append
RECIST.rownum <- (rownames(patient.RECIST))
RECIST[RECIST.rownum,]$purity_mean <- purities
RECIST[RECIST.rownum,]$rat <- rats

##############################################################################

go.purity.001 <- vector(mode = "list", length = length(go_patients))

for(patient_x in go_patients){
  tryCatch({
    cat("Starting patient :", patient_x, "\n")
    go.purity.001[[patient_x]] <- run_liquidCNA(patient_x)
  }, error=function(e){cat("ERROR at:", patient_x, "\n")})
}

##############################################################################

null_patients <- which(sapply(1:length(liquidCNA_results), function(x) is.null(liquidCNA_results[[x]])))
go_patients %in% null_patients
null_patients %in% timely_patients



run_liquidCNA(12)











