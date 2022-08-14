#patient1.seg.data <- seg.data
save(copyNumbersSegmented, patient1.seg.data, file = "../DATA/weekly_july1.RData")

save(go_results,two_samp_results, timely_res1, timely_res2, 
     file = "../DATA/liquidCNA_results_july1.RData")
RECIST

final.results$purity_mean <- as.character(pHat.df[1,match(final.results$time, 
                                                          names(pHat.df))])

RECIST$time

ids <- RECIST$Patient_ID

splt <- split(ids, cumsum(c(1, diff(ids) != 0)))

length(splt)

table( ids )
ichorCNA$