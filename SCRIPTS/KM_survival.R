#threshold <- (mean(yes.met.rat[!(is.na(yes.met.rat))]) - mean(no.met.rat[!(is.na(no.met.rat))]))/2 + mean(no.met.rat[!(is.na(no.met.rat))])
threshold = 0.3656956

survival.df <- data.frame(Patient_ID = RECIST$Patient_ID,
                          First.Meta = RECIST$Date.Meta,
                          Scan.date = RECIST$Date,
                          time = NA,
                          R.Censor = NA,
                          subclone = RECIST$rat,
                          Group = NA,
                          Event = NA)

for(x in 1:nrow(survival.df)){
  ifelse(survival.df[x,]$subclone > threshold, 
         survival.df[x,]$Group <- 2, 
         survival.df[x,]$Group <- 1)
  #Group: 2 = high subclonality group; 1 = low subclonality group
}

survival.df$First.Meta
First.Meta.new <- structure(numeric(length(survival.df$First.Meta)), class="Date")

for(x in 1:length(First.Meta.new)){
  new.date <- as.Date(survival.df$First.Meta[x], format = "%m/%d/%y")
  First.Meta.new[x] <- new.date
}
survival.df$First.Meta <- First.Meta.new 

time <- sapply(1:nrow(survival.df), function(x) survival.df[x,]$Scan.date - survival.df[x,]$First.Meta)
time[which(time < 0)] <- 0
#why are there negative times, what exactly is Date.Meta?
survival.df$time <- time
save(time, file = "../DATA/time.RData")

patient.metastasis <- function(p.id){
  p.RECIST <- RECIST[which(RECIST$Patient_ID == p.id),]
  if(p.RECIST$time[1] <= 10){
    p.RECIST <- p.RECIST[-1,]
  }
  m <- "y" %in% unname(unlist(p.RECIST[,5:14]))
  return(m)
}


ts.metastasis <- function(p.id){
  p.RECIST <- RECIST[which(RECIST$Patient_ID == p.id),]
  if(p.RECIST$time[1] <= 10){
    p.RECIST <- p.RECIST[-1,]
  }
  m <- "y" %in% unname(unlist(p.RECIST[,5:14]))
  return(m)
}

ts.each.metastasis <- sapply(RECIST$Patient_ID, function(x) ts.metastasis(x))
survival.df[which(ts.each.metastasis),]$Event <- 2 #yes metastatsis event
survival.df[which(ts.each.metastasis==F),]$Event <- 1 #no metastatsis event

censor <- rep(NA, nrow(survival.df))
for(i in 1:length(id)){
  ifelse(id[i+1]==id[i], censor[i] <- 1, censor[i] <- 0)
}
tail(censor, 1) <- 1
survival.df$R.Censor <- censor

Surv(survival.df$time, survival.df$R.Censor)[1:10]




model.1 = coxph(Surv(tstart,tstop,status) ~ var1 + var2+ … + vark, 
                method="breslow", 
                robust=TRUE, 
                data = example1)
summary(model.1)

########################################################################
######RECIST wrong date

RECIST[RECIST$Patient_ID==3357,]

#########################################
load("../DATA/liquidCNA_results_aug2.RData")
load("../DATA/antwerpData.RData")
RECIST <- append.RECIST()
survival.df2 <- data.frame(id = RECIST$Patient_ID,
                          tstart=0,
                          tstop=RECIST$time,
                          status=NA,
                          event=NA,
                          sub.group=NA,
                          subclone = RECIST$rat,
                          purity.group=NA,
                          purity = RECIST$purity_mean)
############################
#filling tstart
p.ids <- unique(survival.df2$id)
fill.time <- function(patient_id){
  p.rows <- which(survival.df2$id == patient_id)
  if(length(p.rows) >= 2){
    for(samp in 2:length(p.rows)){
      survival.df2[p.rows[samp],]$tstart <<- survival.df2[p.rows[samp]-1,]$tstop
    }
  } else {
  }
}

sapply(p.ids, function(x) fill.time(x))
############################
#filling status
# 1 is recurrent event
# 0 is censored
id <- survival.df2$id
status <- rep(NA, nrow(survival.df2))
for(i in 1:length(id)){
  ifelse(id[i+1]==id[i], status[i] <- 1, status[i] <- 0)
}
status[length(status)] <- 0
survival.df2$status <- status

############################
#filling event

# ts.metastasis <- function(p.id){
#   p.RECIST <- RECIST[which(RECIST$Patient_ID == p.id),]
#   if(p.RECIST$time[1] <= 10){
#     p.RECIST <- p.RECIST[-1,]
#   }
#   m <- "y" %in% unname(unlist(p.RECIST[,5:14]))
#   return(m)
# }

survival.df2[which(RECIST$New.Met),]$event <- 2 #yes metastatsis event
survival.df2[which(RECIST$New.Met==F),]$event <- 1 #no metastatsis event

############################
#filling sub.group
threshold = 0.3656956
for(x in 1:nrow(survival.df2)){
  ifelse(survival.df2[x,]$subclone > threshold, 
         survival.df2[x,]$sub.group <- 2, 
         survival.df2[x,]$sub.group <- 1)
  #Group: 2 = high subclonality group; 1 = low subclonality group
}

############################
#filling purity.group
for(x in 1:nrow(survival.df2)){
  ifelse(survival.df2[x,]$purity > 0.13, 
         survival.df2[x,]$purity.group <- 2, 
         survival.df2[x,]$purity.group <- 1)
  #Group: 2 = high purity group; 1 = low purity group
}


############################

survival.df2 <- survival.df2[survival.df2$subclone<=1,]

####################################
# model.1 = coxph(Surv(tstart,tstop,status) ~ (subclone + purity + cluster(id)), 
#                 method="breslow", 
#                 robust=TRUE, 
#                 data = survival.df2)
# summary(model.1)
# 
# plot(survfit(model.1), 
#      xlab = "Days", 
#      ylab = "Overall probability of event (metastasis)")
# 
# summary(survfit(model.1), times = 365.25)
# 
# 
# sfit <- survfit(model.1)
# cumhaz.upper <- -log(sfit$upper)
# cumhaz.lower <- -log(sfit$lower)
# cumhaz <- sfit$cumhaz # same as -log(sfit$surv)
# 
# plot(cumhaz, xlab="Days ahead", ylab="Cumulative hazard",
#      ylim=c(min(cumhaz.lower), max(cumhaz.upper)))
# lines(cumhaz.lower)
# lines(cumhaz.upper)

############################

model.2 = coxph(Surv(tstart,tstop,status) ~ (sub.group + purity.group + cluster(id)), 
                method="breslow", 
                robust=TRUE, 
                data = survival.df2)
summary(model.2)

sfit <- survfit(model.2)
cumhaz.upper <- -log(sfit$upper)
cumhaz.lower <- -log(sfit$lower)
cumhaz <- sfit$cumhaz # same as -log(sfit$surv)

pdf(file="Fig_26.pdf", width = 5.43, height = 5.1)
plot(cumhaz, xlab="Days ahead", ylab="Cumulative hazard",
     ylim=c(min(cumhaz.lower), max(cumhaz.upper)),
     pch=21, main = "Cumulative hazard of recurrent metastasis")
lines(cumhaz.lower)
lines(cumhaz.upper)
dev.off()

model.df <- data.frame(Hazard_ratio = c(1.2434 , 0.8246),
                       lower_95 = c(0.7992, 0.5083),
                       upper_95 = c(1.935, 1.338),
                       p_value = c(0.334, 0.435))

rownames(model.df) <- c("Subclonality", "Purity")

pdf(file="Fig_26_tab.pdf", width = 5.43, height = 1.6)
ss <- tableGrob(model.df)
grid.arrange(ss)
dev.off()
 
##################################
model.3 = coxph(Surv(tstart,tstop,status) ~ (sub.group + cluster(id)), 
                method="breslow", 
                robust=TRUE, 
                data = survival.df2)
summary(model.3)

model.4 = coxph(Surv(tstart,tstop,status) ~ (purity.group + cluster(id)), 
                method="breslow", 
                robust=TRUE, 
                data = survival.df2)
summary(model.4)
