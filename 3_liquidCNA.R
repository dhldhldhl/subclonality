###############################################################################
library(pracma); library(ggplot2); library(ggpubr); library(reshape2); 
library(mixtools); library(fitdistrplus); library(dplyr); library(QDNAseq); 
require(gtools); library(gridExtra); library(devtools)
source_url("https://raw.githubusercontent.com/elakatos/liquidCNA/main/mixture_estimation_functions.R")
###############################################################################

load("bam_by_patient.RData")
p_vec <- 1:length(patient_ids)
p <- p_vec[1]
#read outputs of QDNAseq

seg.df <- read.delim(paste0("2_QDNA_CNout/segment_", patient_ids[p], ".txt"))
cn.df <- read.delim(paste0("2_QDNA_CNout/raw_cn_", patient_ids[p], ".txt"))

seg.df <- seg.df[,-(1:4)]
cn.df <- cn.df[,-(1:4)]

colnames(seg.df) <- paste0("Sample", 1:ncol(seg.df))
colnames(cn.df) <- paste0("Sample", 1:ncol(cn.df))

