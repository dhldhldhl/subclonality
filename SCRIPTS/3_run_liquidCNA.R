library(pracma); library(ggplot2); library(ggpubr); library(reshape2); 
library(mixtools); library(fitdistrplus); library(dplyr); library(QDNAseq); 
require(gtools); library(gridExtra); library(devtools)
source_url("https://raw.githubusercontent.com/elakatos/liquidCNA/main/mixture_estimation_functions.R")

################################################################################

#setwd & load 
setwd("/Users/dhl/laheen Dropbox/DHL DHL/Cambridge/Internship/subclonality/")

#setwd("/rds/project/rds-csoP2nj6Y6Y/ctDNA/dl/subclonality")
load("../DATA/bam_by_patient.RData")

#function run_liquidCNA
source("SCRIPTS/3_liquidCNA.R")

################################################################################

#initialise list to store results
liquidCNA_results <- vector(mode = "list", length = 80)

################################################################################

#Categorise patients

#Patients with only two samples
sample_num_patient <- sapply(1:80, function(x) length(bam_by_patient[[x]]))
two_samp_patients <- which(sample_num_patient == 2)

#Patients with over 6 samples
timely_patients <- which(sample_num_patient > 6)

go_patients <- 1:80
go_patients <- go_patients[!(go_patients %in% two_samp_patients)]
go_patients <- go_patients[!(go_patients %in% timely_patients)]

################################################################################

#Run for go_patients (i.e., non 2 or timely sample patients)
#EMsteps should be set to 1500 to not get NAs for absolute ratios

for(patient_x in go_patients){
  tryCatch({
    cat("Starting patient :", patient_x, "\n")
    liquidCNA_results[[patient_x]] <- run_liquidCNA(patient_x)
  }, error=function(e){cat("ERROR at:", patient_x, "\n")})
}

################################################################################
#Run for timely patients

#number of BAM files our timely patients have
timely_bam_num <- sample_num_patient[timely_patients]

#number of bam files to batch each timely patient by:
nbatch1 <- ceiling(timely_bam_num/2)
nbatch2 <- timely_bam_num - nbatch1

#Batch 1
#initialise result vector
timely_res1 <- vector(mode = "list", length = length(timely_patients))

#RUN!
for(x in 1:length(timely_patients)){
  tryCatch({
    cat("\n Starting", x, "\n")
    timely_res1[[x]] <- run_liquidCNA(timely_patients[x], 
                                      timely = TRUE,
                                      batch1 = TRUE,
                                      nbatch = nbatch1[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

names(timely_res1) <- paste0("patient_",patient_ids[timely_patients],"_batch1")

#option3: Sample1 is baseSample for both batch
# timely_res2 <- vector(mode = "list", length = length(timely_patients))
# 
# for(x in 1:length(timely_patients)){
#   tryCatch({
#     cat("\n Starting", x, "\n")
#     base.id <- tail(timely_res1[[x]]$time,1)
#     base.id <- as.numeric(strsplit(base.id, "e")[[1]][2])
#     timely_res2[[x]] <- run_liquidCNA(timely_patients[x], 
#                                            timely = TRUE,
#                                            batch1 = FALSE,
#                                            nbatch = nbatch2[x])
#   }, error=function(e){cat("ERROR at:", x, "\n")})
# }
# 
# names(timely_res2) <- paste0("patient_",patient_ids[timely_patients],"_batch2")
# 

timely_res2.opt3 <- vector(mode = "list", length = length(timely_patients))

for(x in 1:length(timely_patients)){
  tryCatch({
    base.id <- tail(timely_res1[[x]]$time,1)
    base.id <- as.numeric(strsplit(base.id, "e")[[1]][2])
    cat("\n Starting", x, "\n")
    timely_res2.opt3[[x]] <- run_liquidCNA(timely_patients[x], 
                                           timely = TRUE,
                                           batch1 = FALSE,
                                           nbatch = nbatch2[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

names(timely_res2.opt3) <- paste0("patient_",patient_ids[timely_patients],"_batch2")


#Stitch

timely_res <- sapply(1:length(timely_patients), function(x) rbind(timely_res1[[x]][-which(timely_res1[[x]]$relratio == 0),], timely_res2.opt3[[x]]),
                     simplify = FALSE)
names(timely_res) <- paste0("patient_",patient_ids[timely_patients])

for(x in 1:length(timely_patients)){
  liquidCNA_results[[timely_patients[x]]] <- timely_res[[x]]
}

################################################################################
mean(unlist(sapply(1:length(liquidCNA_results), function(x) liquidCNA_results[[x]]$cutOff)))
#mean cutOff is 0.17
median(unlist(sapply(1:length(liquidCNA_results), function(x) liquidCNA_results[[x]]$cutOff)))
#median cutOff is 0.16

#cut chosen as 0.16

#Run for 2 time patients
#Run!
for(patient_x in two_samp_patients){
  tryCatch({
    cat("Starting patient :", patient_x, "\n")
    liquidCNA_results[[patient_x]] <- run_2samp_liquidCNA(patient_x)
  }, error=function(e){cat("ERROR at:", patient_x, "\n")})
}


################################################################################

names(liquidCNA_results) <- paste0("patient_", patient_ids)
save(liquidCNA_results, file = "../DATA/liquidCNA_results_aug2.RData")

load(file = "../DATA/liquidCNA_results_aug1.RData")
################################################################################

timely_res1

#RUN!
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


# #RUN!
# timely_res2.opt2 <- vector(mode = "list", length = length(timely_patients))
# 
# for(x in 1:length(timely_patients)){
#   tryCatch({
#     cat("\n Starting", x, "\n")
#     
#     batch1.rats <- as.numeric(timely_res1[[x]]$rat)
#     batch1.max.rat <- timely_res1[[x]]$time[which(batch1.rats == max(batch1.rats))]
#     base.id <- as.numeric(strsplit(batch1.max.rat, "e")[[1]][2])
#     
#     timely_res2.opt2[[x]] <- run_liquidCNA(timely_patients[x], 
#                                            timely = TRUE,
#                                            batch1 = FALSE,
#                                            nbatch = nbatch2[x])
#     
#   }, error=function(e){cat("ERROR at:", x, "\n")})
# }
# 
# names(timely_res2.opt2) <- paste0("patient_",patient_ids[timely_patients],"_batch2")

#RUN!
timely_res2.opt3 <- vector(mode = "list", length = length(timely_patients))

for(x in 1:length(timely_patients)){
  tryCatch({
    base.id <- tail(timely_res1[[x]]$time,1)
    base.id <- as.numeric(strsplit(base.id, "e")[[1]][2])
    cat("\n Starting", x, "\n")
    timely_res2.opt3[[x]] <- run_liquidCNA(timely_patients[x], 
                                      timely = TRUE,
                                      batch1 = FALSE,
                                      nbatch = nbatch2[x])
  }, error=function(e){cat("ERROR at:", x, "\n")})
}

names(timely_res2.opt3) <- paste0("patient_",patient_ids[timely_patients],"_batch2")



timely_res.opt1 <- sapply(1:length(timely_patients), function(x) rbind(timely_res1[[x]][-which(timely_res1[[x]]$relratio == 0),], timely_res2.opt1[[x]]), simplify = F)
names(timely_res.opt1) <- paste0("patient_",patient_ids[timely_patients])

# timely_res.opt2 <- sapply(1:length(timely_patients), function(x) rbind(timely_res1[[x]][-which(timely_res1[[x]]$relratio == 0),], timely_res2.opt2[[x]]), simplify = F)
# names(timely_res.opt2) <- paste0("patient_",patient_ids[timely_patients])

timely_res.opt3 <- sapply(1:length(timely_patients), function(x) rbind(timely_res1[[x]][-which(timely_res1[[x]]$relratio == 0),], timely_res2.opt3[[x]]), simplify = F)
names(timely_res.opt3) <- paste0("patient_",patient_ids[timely_patients])


x=6
base.to.first <- function(rats){
  numeric.rat <- as.numeric(rats)
  base <- tail(numeric.rat,1)
  non.base <- numeric.rat[1:(length(numeric.rat)-1)]
  ordered <- c(base, non.base)
  return(ordered)
}
##
opt1_ordering <- timely_res.opt1[[x]]$time
non.base <- opt1_ordering[1:(length(opt1_ordering)-1)]
opt1_ordering <- c("Sample2", non.base)

opt3_ordering <- timely_res.opt3[[x]]$time
base <- tail(opt3_ordering,1)
non.base <- opt3_ordering[1:(length(opt3_ordering)-1)]
opt3_ordering <- c(base, non.base)

##
pdf(file="Fig_stitch_opt.pdf", width = 7, height = 4.5)
par(mfrow = c(1,2))

plot(base.to.first(timely_res.opt1[[x]]$rat), type = "l", xaxt = "n", 
     main = "\n \n Option 1 & 2",
     ylab = "Sub-clonal ratio estimates", xlab = "Sample Order")
axis(1, at=1:8, labels=opt1_ordering)
abline(v=nbatch1[x]+0.5, lty = 2, col = "red")
plot(base.to.first(timely_res.opt3[[x]]$rat), xaxt = "n", type = "l",
     main = "\n \n Option 3", 
     ylab = "Sub-clonal ratio estimates", xlab = "Sample Order")
abline(v=nbatch1[x]+0.5, lty = 2, col = "red")
axis(1, at=1:8, labels=opt3_ordering)
title("\n \n Comparison of stitching options | Patient 3068", outer = T)
dev.off()
