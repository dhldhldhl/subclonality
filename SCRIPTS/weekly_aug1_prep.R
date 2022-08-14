#Calculating time

###################
#First meta
First.Meta <- RECIST$Date.Meta
First.Meta.new <- structure(numeric(length(First.Meta)), class="Date")

for(x in 1:length(First.Meta.new)){
  new.date <- as.Date(survival.df$First.Meta[x], format = "%m/%d/%y")
  First.Meta.new[x] <- new.date
}

###################
#Scan date
Scan.date <- RECIST$Date

time <- Scan.date - First.Meta.new

###############################################################################
#How many recurrent metastasis?
#There were 39 metastasis
table(RECIST$New.Met)

#a function which for each patient counts new metastasis == TRUE
recurrent.met <- function(p.id){
  p.RECIST <- RECIST[RECIST$Patient_ID==p.id,]
  number <- sum(p.RECIST$New.Met, na.rm = TRUE)
  return(number)
}


recur <- sapply(unique(RECIST$Patient_ID), function(x) recurrent.met(x))

table(recur)
#of the 39 metastasis
#26 are secondary
#5 are tertiary
###############################################################################

