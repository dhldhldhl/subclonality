#setwd
wd <- "/rds/project/rds-csoP2nj6Y6Y/ctDNA/dl"
setwd(wd)

#package directory
pkg_dir <- "/rds/project/rds-csoP2nj6Y6Y/ctDNA/dl/packages"

#install packages
if (!require("future", lib.loc = pkg_dir, quietly = TRUE)) {
  install.packages("future", lib = pkg_dir,repos = "http://cran.us.r-project.org")
}
if (!require("BiocManager", lib.loc = pkg_dir, quietly = TRUE)) {
  install.packages("BiocManager", lib = pkg_dir, repos = "http://cran.us.r-project.org")
}
if (!require("devtools", lib.loc = pkg_dir, quietly = TRUE)) {
  install.packages("devtools", lib = pkg_dir, repos = "http://cran.us.r-project.org")
}
BiocManager::install("QDNAseq")
devtools::install_github("asntech/QDNAseq.hg38@main")

#For liquidCNA
packages <- c("ggplot2", "pracma", "ggpubr", "reshape2", "mixtools", 
              "fitdistrplus", "dplyr", "gtools", "gridExtra")
install.packages(setdiff(packages, rownames(installed.packages())))  

library(pracma); library(ggplot2); library(ggpubr); library(reshape2); 
library(mixtools); library(fitdistrplus); library(dplyr); library(QDNAseq); 
require(gtools); library(gridExtra); library(devtools)

