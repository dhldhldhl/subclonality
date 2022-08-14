##################################################################
#P-value

#T-test for between subclonality of YES or NO progression
no.prog <- which(RECIST$Progression == "NO")
yes.prog <- which(RECIST$Progression == "YES")

no.prog.rat <- RECIST[no.prog,]$rat
yes.prog.rat <- RECIST[yes.prog,]$rat

no.prog.rat <- no.progression[!(is.na(no.progression))]
yes.prog.rat <- yes.progression[!(is.na(yes.progression))]

library("dplyr")
library("ggpubr")

shapiro.test(no.prog.rat) # p-value = 2.928e-09
shapiro.test(yes.prog.rat) # p-value < 2.2e-16
# p-value < 0.05 implying that the distribution of the data 
# are significantly different from normal distribution
# Distribution of estimated subclonalities are skewed
# normality should not be assumed

#so t-test shouldn't be done, but let's see what happens
t.test(no.prog.rat, yes.prog.rat, var.equal = FALSE) #p-value = 0.1091

#as normality cannot be assumed,
#do Unpaired Two-Samples Wilcoxon Test 
wilcox.rat <- wilcox.test(no.prog.rat, yes.prog.rat, alternative = "two.sided")
wilcox.rat$p.value

#############################
shapiro.test(RECIST[no.prog,]$purity_mean) #p-value = 7.904e-14
shapiro.test(RECIST[yes.prog,]$purity_mean) #p-value = 1.398e-12
#for purity as well normality cannot be assumed

t.test(RECIST[no.prog,]$purity_mean, RECIST[yes.prog,]$purity_mean, var.equal = FALSE)
#p-value = 0.003662
wilcox.purity <- wilcox.test(RECIST[no.prog,]$purity_mean, RECIST[yes.prog,]$purity_mean, alternative = "two.sided")
wilcox.purity$p.value

#######################################################################################
#Subclonality for estimating Metastasis

########################### by patient...
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
p.met.res <- wilcox.test(yes.met.rat, no.met.rat, alternative = "two.sided")

#draw
boxplot(no.met.rat, 
        yes.met.rat,
        ylim = c(0, 3),
        main = "Subclonality against patients with new metastasis",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "Subclonal ratio estimation")
text(2.1, 2.9, paste0("Wilcox test: \np-value = ",round(p.met.res$p.value, 5)))



########################### by time sample...
ts.met <- function(row_num){
  p.RECIST <- RECIST[row_num,]
  if(p.RECIST$time[1] <= 10){
    m = FALSE
  } else {
    m <- "y" %in% unname(unlist(p.RECIST[,5:14]))
  }
  return(m)
}

ts.metastasis <- sapply(1:nrow(RECIST), function(x) ts.met(x))

yes.metastasis <- RECIST[ts.metastasis, ]$rat
no.metastasis <- RECIST[ts.metastasis==FALSE, ]$rat

ts.met.res <- wilcox.test(yes.metastasis, no.metastasis, alternative = "two.sided")

#draw
boxplot(no.metastasis, 
        yes.metastasis,
        ylim = c(0, 3),
        main = "Subclonality against new metastasis (by time sample)",
        at = c(1,2),
        names = c("NO", "YES"),
        las = 1,
        col = c("pink","lightblue"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE,
        ylab = "Subclonal ratio estimation")
text(2.1, 2.9, paste0("Wilcox test: \np-value = ",round(ts.met.res$p.value, 5)))


############### with purity...

table(is.na(as.numeric(unlist(sapply(1:80, function(x) liquidCNA_results[[x]]$rat)))))


#################################################
ts = 3
plot.3C(ts, cn.df, seg.df, seg.sub, seg.plot,  patient_id = 1005)
abline(v=365, lty = 2, col = "red"); abline(v=370, lty = 2, col = "red"); abline(v=375, lty = 2, col = "red")


