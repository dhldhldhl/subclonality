load("../DATA/weekly_july3.RData")
save(pVec.46, final.results.46, seg.cns.46, seg.cns.corr.46, 
     seg.dcn.46,seg.dcn.toUse.46, cn.df.58, seg.df.58, seg.sub.58,
     seg.plot.58, seg.dcn.toUse.58, pVec.58, c.df,
     file = "../DATA/weekly_july3.RData")

p_num <- match(3080, patient_ids)
pVec.46 <- pVec
final.results.46 <- final.results
seg.cns.46 <- seg.cns
seg.cns.corr.46 <- seg.cns.corr
seg.dcn.46 <- seg.dcn
seg.dcn.toUse.46 <- seg.dcn.toUse
samp <- "Sample3" #"Sample2"

p_num <- match(3301, patient_ids)
pVec.58 <- pVec
cn.df.58 <- cn.df
seg.df.58 <- seg.df
seg.sub.58 <- seg.sub
seg.plot.58 <- seg.plot
seg.dcn.toUse.58 <- seg.dcn.toUse

#########################
rats <- unlist(sapply(1:length(liquidCNA_results), function(x) as.numeric(liquidCNA_results[[x]]$rat)))
purities <- unlist(sapply(1:length(liquidCNA_results), function(x) as.numeric(liquidCNA_results[[x]]$purity_mean)))

# 31/283 estimations have sub-clonality ratio outside of [0,1]
out.rats<- which(rats < 0 | rats > 1)
length(out.rats)

is.out.rat <- function(x){
  rats <- as.numeric(liquidCNA_results[[x]]$rat)
  return(any(rats < 0 | rats > 1))
}

# 22/80 patients have sub-clonality ratio outside of [0,1]
out.rat.patients <- which(sapply(1:length(liquidCNA_results), function(x) is.out.rat(x)))
length((out.rat.patients))
sample_num_patient[out.rat.patients]

#purities of out.rats
rats[out.rats]
purities[out.rats]

plot(purities[out.rats], rats[out.rats],
     ylab = "Sub-ratio estimate outside [0,1]",
     xlab = "Corresponding purity estimates")

out.results <- sapply(out.rat.patients, function(x) liquidCNA_results[[x]], simplify = F)
names(out.results) <- paste0("patient_", patient_ids[out.rat.patients])

##############################################################################

##############################################################################

###############################################################################

##############################################################################
boxplot(RECIST$rat, 
        main = "Boxplot of all subclonality estimations",
        las = 1,
        col = c("lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE)

##############################################################################
#SUBCLONALITY X PROGRESSION!
no.prog <- which(RECIST$Progression == "NO")
yes.prog <- which(RECIST$Progression == "YES")

library("dplyr")
library("ggpubr")
no.prog.rat <- RECIST[no.prog,]$rat
yes.prog.rat <- RECIST[yes.prog,]$rat

no.prog.rat <- no.prog.rat[!(is.na(no.prog.rat))]
yes.prog.rat <- yes.prog.rat[!(is.na(yes.prog.rat))]

no.prog.rat <- no.prog.rat[no.prog.rat<=1]
yes.prog.rat <- yes.prog.rat[yes.prog.rat<=1]

#as normality cannot be assumed,
#do Unpaired Two-Samples Wilcoxon Test 
wilcox.rat <- wilcox.test(no.prog.rat, yes.prog.rat, alternative = "two.sided")
wilcox.rat$p.value

wilcox.purity <- wilcox.test(RECIST[no.prog,]$purity_mean, RECIST[yes.prog,]$purity_mean, alternative = "two.sided")
wilcox.purity$p.value


pdf(file="Fig_2_PROG.pdf", width = 7, height = 5)
par(mfrow=c(1,2))
#
boxplot(no.prog.rat,
        yes.prog.rat, 
        ylim = c(0, 1.15),
        main = "(A) Subclonality estimation \n against progression",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink", "lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "Subclonal ratio estimation")
text(1.95, 1.1, paste0("Wilcox p-value: \n ", round(wilcox.rat$p.value, 4)))

#PURITY X PROGRESSION!
boxplot(RECIST[no.prog,]$purity_mean,
        RECIST[yes.prog,]$purity_mean, 
        ylim = c(0, 0.6),
        main = "(B) Purity estimation \n against progression",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "liquidCNA purity estimation")
text(1.95, 0.57, paste0("Wilcox p-value: \n ", round(wilcox.purity$p.value, 4)))

dev.off()
###########################################################
no.prog <- which(RECIST$Progression == "NO")
yes.prog <- which(RECIST$Progression == "YES")

yes.prog.rat <- RECIST[yes.prog, ]$rat
no.prog.rat <- RECIST[no.prog, ]$rat

yes.prog.pur <- RECIST[yes.prog, ]$purity_mean
no.prog.pur <- RECIST[no.prog, ]$purity_mean

yes.prog.rat <- yes.prog.rat[yes.prog.pur>=0.13]
no.prog.rat <- no.prog.rat[no.prog.pur>=0.13]

yes.prog.rat <- yes.prog.rat[!is.na(yes.prog.rat)]

#over 1 is uninterpretable
yes.prog.rat <-yes.prog.rat[yes.prog.rat<=1]
no.prog.rat <- no.prog.rat[no.prog.rat<=1]





pdf(file="Fig_S1.pdf", width = 4.3, height = 5.3)
boxplot(yes.prog.rat,
        no.prog.rat, 
        ylim = c(0, 1),
        main = "Subclonality estimation \n against progression \n (samples with purity  >= 0.13)",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = F,
        notch = F,
        ylab = "liquidCNA purity estimation")
text(1.95, 0.95, paste0("Wilcox p-value: \n ", round(wilcox.test(yes.prog.rat, no.prog.rat, alternative = "two.sided")$p.value, 4)))

dev.off()

p.met.res <- wilcox.test(yes.met.rat, no.met.rat, alternative = "two.sided")
##############################################################################
#SUBCLONALITY X METASTASIS!
patient.metastasis <- function(p.id){
  p.RECIST <- RECIST[which(RECIST$Patient_ID == p.id),]
  if(p.RECIST$time[1] <= 10){
    p.RECIST <- p.RECIST[-1,]
  }
  m <- "y" %in% unname(unlist(p.RECIST[,5:14]))
  return(m)
}

p.metastasis <- sapply((RECIST$Patient_ID), function(x) patient.metastasis(x))

yes.met.rat <- RECIST[p.metastasis, ]$rat
no.met.rat <- RECIST[p.metastasis==FALSE, ]$rat

yes.met.rat <- yes.met.rat[!(is.na(yes.met.rat))]
no.met.rat <- no.met.rat[!(is.na(no.met.rat))]

#over 1 is uninterpretable
yes.met.rat <-yes.met.rat[yes.met.rat<=1]
no.met.rat <- no.met.rat[no.met.rat<=1]

p.met.res <- wilcox.test(yes.met.rat, no.met.rat, alternative = "two.sided")

#PURITY X METASTASIS!
p.met.purity.res <- wilcox.test(RECIST[p.metastasis,]$purity_mean,
                                RECIST[p.metastasis==FALSE,]$purity_mean,
                                alternative = "two.sided")



###!draw
pdf(file="Fig_3_METASTAT.pdf", width = 7, height = 5)

par(mfrow=c(1,2))
#draw
boxplot(no.met.rat, 
        yes.met.rat,
        ylim = c(0, 1.15),
        main = "(A) Subclonality against \n new metastasis",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "Subclonal ratio estimation")
text(1.9, 1.1, paste0("Wilcox test: \np-value = ", round(p.met.res$p.value, 4)))

boxplot(RECIST[p.metastasis,]$purity_mean,
        RECIST[p.metastasis==FALSE,]$purity_mean, 
        ylim = c(0, 0.6),
        main = "(B) Purity against \n new metastasis",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "liquidCNA purity estimation")
text(1.9, 0.58, paste0("Wilcox test: \np-value = ", round(p.met.purity.res$p.value, 4)))

dev.off()
#################################


##############################################################################
pdf(file="Fig_24a.pdf", width = 6.1, height = 5)
plot.estimation.w.RECIST(match(2734,patient_ids))
dev.off()

pdf(file="Fig_24b.pdf", width = 6.1, height = 5)
plot.estimation.w.RECIST(match(2977,patient_ids))
dev.off()


 shapiro.test(no.prog.rat)
shapiro.test(yes.prog.rat)
# p-value < 0.05 implying that the distribution of the data 
# are significantly different from normal distribution
#Distribution of estimated subclonalities are skewed
# normality should not be assumed

 ##############################################################################
#get purity estimations from new ichorCNA sent by solon
new.ichor <- read.csv("../DATA/ichor_estimates.csv")

mean(new.ichor$read_count)
median(new.ichor$read_count)

# bam_by_patient.id <- sapply(1:length(bam_by_patient), function(x) unlist(strsplit(bam_by_patient[[x]], ".bam")))
# names(bam_by_patient.id) <- paste0("patient_", patient_ids)
# get.ichor.purities <- function(x){
#   sample.row <- match(bam_by_patient.id[[x]], new.ichor$id)
#   sample.ichor.purity <- new.ichor[sample.row,]$ichorCNA
#   return(sample.ichor.purity)
# }
# ichor.purities <- sapply(1:80, function(x) get.ichor.purities(x))
# names(ichor.purities) <- paste0("patient_", patient_ids)

#for each patient get 'key'
getKey <- function(p_num){
  split <- strsplit(bam_by_patient[[p_num]], ".H")
  keys <- sapply(1:length(split), function(x) split[[x]][1])
  return(keys)
}

patient.keys <- sapply(1:80, function(x) getKey(x))
names(patient.keys) <- paste0("patient_", patient_ids)

use.key <- function(x){
  sample.row <- match(patient.keys[[x]], new.ichor$key)
  sample.ichor.purity <- new.ichor[sample.row,]$ichorCNA
  return(sample.ichor.purity)
}

ichor.purities <- sapply(1:80, function(x) use.key(x))
names(ichor.purities) <- paste0("patient_", patient_ids)

#liquidCNA purities
purities <- (sapply(1:length(liquidCNA_results), 
                          function(x) as.numeric(liquidCNA_results[[x]]$purity_mean)))
names(purities) <- paste0("patient_", patient_ids)

pdf(file="Fig_6.pdf", width = 5.43, height = 6)
plot(unlist(purities), unlist(ichor.purities), xlab = "liquidCNA estimation", ylab = "ichorCNA estimation",
     main = "Purity Estimation | liquidCNA vs ichorCNA",
     ylim = c(0,0.8), xlim = c(0,0.8), col = "darkgreen", pch = 3)
rsqrd <- rsq(unlist(purities), unlist(ichor.purities))
text(0.725, 0.78, paste0("r2 = ",round(rsqrd, 5)))
dev.off()
##############################################################################
#ichorCNA PURITY X PROGRESSION!
yes.prog <- which(RECIST$Progression == "YES")

par(mfrow=c(1,2))
#draw
boxplot(RECIST[no.prog,]$purity_mean,
        RECIST[yes.prog,]$purity_mean, 
        ylim = c(0, 0.6),
        main = "Purity estimation against progression: \n liquidCNA",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "liquidCNA purity estimation")

#draw
boxplot(RECIST[-yes.prog,]$ichorCNA,
        RECIST[yes.prog,]$ichorCNA, 
        ylim = c(0, 0.8),
        main = "\n ichorCNA",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "ichorCNA purity estimation")

#########################################################
rats <- as.numeric(unlist(sapply(1:80, function(x) liquidCNA_results[[x]]$rat)))
purities <- as.numeric(unlist(sapply(1:80, function(x) liquidCNA_results[[x]]$purity_median)))

na.rats <- which(is.na(rats))
rats <- rats[-na.rats]
purities <- purities[-na.rats]

#calculate coefficient of correlation
rsq <- function (x, y) cor(x, y) ^ 2
subclone.purity.rsqrd <- rsq(rats, purities)

#plot
pdf(file="Fig_14.pdf", width = 5.43, height = 5.43)
plot(rats, purities, main = "Purity estimates against subclonality",
     ylab = "Purity estimates", xlab = "Subclonality estimates", 
     col = "darkred", 
     #xlim = c(0,3), 
     pch = 3)
text(2.35, 0.5, paste0("r2 = ",round(subclone.purity.rsqrd, 5)))
dev.off()
#############################################################
pdf(file="Fig_15b.pdf", width = 5.43, height = 7)
par(mfrow=c(3,1))
plot.3C(1, cn.df, seg.df, seg.sub, seg.plot, patient_id = 1005)
abline(v=370, lwd = 2.5, lty = 2, col = "red")
plot.3C(2, cn.df, seg.df, seg.sub, seg.plot, patient_id = 1005)
abline(v=370, lwd = 2.5, lty = 2, col = "red")
plot.3C(3, cn.df, seg.df, seg.sub, seg.plot, patient_id = 1005)
abline(v=370, lwd = 2.5, lty = 2, col = "red")
dev.off()

#############################################################

load("../DATA/liquidCNA_results_july4.RData")
load("../DATA/liquidCNA_results_aug2.RData")
rats <- unlist(sapply(1:length(liquidCNA_results), 
                      function(x) as.numeric(liquidCNA_results[[x]]$rat)))
purities <- unlist(sapply(1:length(liquidCNA_results), 
                          function(x) as.numeric(liquidCNA_results[[x]]$purity_mean)))
out.rats<- which(rats < 0 | rats > 1)
length(out.rats)
#before: 25/283 estimates have sub-clonality ratio outside of [0,1]
#now: 16

is.out.rat <- function(x){
  rats <- as.numeric(liquidCNA_results[[x]]$rat)
  return(any(rats < 0 | rats > 1))
}

# this is 19/80 patients have sub-clonality ratio outside of [0,1]
out.rat.patients <- which(sapply(1:length(liquidCNA_results), function(x) is.out.rat(x)))
length((out.rat.patients))

#purities of out.rats
rats[out.rats]
purities[out.rats]

out.results <- sapply(out.rat.patients, function(x) liquidCNA_results[[x]], simplify = F)
names(out.results) <- paste0("patient_", patient_ids[out.rat.patients])

pdf(file="Fig_18_bigSubclones.pdf", width = 5.43, height = 5.5)
plot(purities[out.rats], rats[out.rats],
     ylab = "Estimated subclonal ratio",
     xlab = "Corresponding purity estimates",
     main = "Subclonality estimates outside [0,1] \n against their respective purity",
     col = "slateblue", pch = 17)
dev.off()

##############
load("../DATA/liquidCNA_results_july4.RData")
df3080 <- data.frame(liquidCNA_results$patient_3080$time,
                 round(as.numeric(liquidCNA_results$patient_3080$rat),4),
                 liquidCNA_results$patient_3080$purity_mean)
df3080 <- rbind(df3080[3,], df3080[1:2,])
colnames(df3080) <- c("Sample order",
                      "Subclonality",
                      "Purity")
rownames(df3080) <- c("(baseline)",
                      " " ,
                      "  ")


library(grid)
library(gridExtra)
library(gtable)

t1 <- tableGrob(df3080)
title <- textGrob("Patient 3080 estimates",gp=gpar(fontsize=14))
padding <- unit(5,"mm")

table <- gtable_add_rows(
  t1, 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 2, 1, ncol(table))

grid.newpage()
grid.draw(table)

pdf(file="Tab1.pdf", width = 5.43, height = 2)
grid.newpage()
grid.draw(table)
dev.off()

###################
load("../DATA/weekly_july3.RData")
t1 <- tableGrob(round(head(seg.cns.46),3), rows = NULL)
title1 <- textGrob("(A) Segment CN values",gp=gpar(fontsize=14))
padding <- unit(5,"mm")

table1 <- gtable_add_rows(
  t1, 
  heights = grobHeight(title1) + padding,
  pos = 0)
table1 <- gtable_add_grob(
  table1, 
  title1, 
  1, 1, 1, ncol(table1))

missed1 <- convertWidth(sum(table1$widths), "in", valueOnly = TRUE) -
  convertWidth(grobWidth(title1), "in", valueOnly = TRUE)

if(missed1 < 0 ) # need to do something about it
  table1$widths <- table1$widths + unit(abs(missed1)/ncol(table1), "in")

grid.newpage()
grid.draw(table1)

###

t2 <- tableGrob(round(head(seg.cns.corr.46),3), rows = NULL)
title2 <- textGrob("(B) Purity-corrected segment CN",gp=gpar(fontsize=14))
padding <- unit(5,"mm")

table2 <- gtable_add_rows(
  t2, 
  heights = grobHeight(title2) + padding,
  pos = 0)
table2 <- gtable_add_grob(
  table2, 
  title2, 
  1, 1, 1, ncol(table2))

missed2 <- convertWidth(sum(table2$widths), "in", valueOnly = TRUE) -
  convertWidth(grobWidth(title2), "in", valueOnly = TRUE)

if(missed2 < 0 ) # need to do something about it
  table2$widths <- table2$widths + unit(abs(missed2)/ncol(table2), "in")

grid.newpage()
grid.draw(table2)

###

t1 <- tableGrob(round(head(seg.dcn.46),3), rows = NULL)
title <- textGrob("(C) Segment dCN values",gp=gpar(fontsize=14))
padding <- unit(5,"mm")

table3 <- gtable_add_rows(
  t1, 
  heights = grobHeight(title) + padding,
  pos = 0)
table3 <- gtable_add_grob(
  table3, 
  title, 
  1, 1, 1, ncol(table3))

missed <- convertWidth(sum(table3$widths), "in", valueOnly = TRUE) -
  convertWidth(grobWidth(title), "in", valueOnly = TRUE)

if(missed < 0 ) # need to do something about it
  table3$widths <- table$widths + unit(abs(missed)/ncol(table3), "in")

grid.newpage()
grid.draw(table3)

pdf(file="Fig_19_lowPurity.pdf", width = 9, height = 4)
grid.arrange(table1, table2, table3, nrow=1, ncol=3)
dev.off()
   #########################################################################

plot.estimation.w.RECIST(match(3188, patient_ids))

###################################################################
p_num = match(3269, patient_ids)
p.id = 3269
pdf(file="Fig_20.pdf", width = 9, height = 4.5)
ts = 3
par(mfrow=c(1,2))
plot(cn.df[,ts], col = "gray", xlab = "Bin", ylab = "CN states",
     main= paste0("\n \n  Sample", ts," with low purity (pi = 0.13)"))
ts = 1
plot(cn.df[,ts], col = "gray", xlab = "Bin", ylab = "CN states", 
     ylim = c(1.5, 10),
     main= paste0("\n \n Sample", ts," with high purity (pi = 0.48)"))
title("\n \n CN profile of time samples from Patient 3269", outer = T)
dev.off()

###################################################################
p_num = match(2734, patient_ids)




