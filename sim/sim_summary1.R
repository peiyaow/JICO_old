load("/Users/peiyaow/Documents/GitHub/continuum/result_PCR.RData")
summary = Reduce("+", RESULT)/50
row.names(summary) = c("2step", "iter", "iter.orthIndiv", "PCR.global", "PLS.global", "PCR.group", "PLS.group")
summary
RANK

load("/Users/peiyaow/Documents/GitHub/continuum/result_PCR_tune.RData")
summary = (Reduce("+", RESULT)/10)[c(10, 20, 30, 31:34),]
row.names(summary) = c("2step", "iter", "iter.orthIndiv", "PCR.global", "PLS.global", "PCR.group", "PLS.group")
summary
RANK

load("/Users/peiyaow/Documents/GitHub/continuum/result_sim.RData")
summary = (Reduce("+", RESULT)/10)
row.names(summary) = c("2step", "PCR.global", "PLS.global", "PCR.group", "PLS.group")
RANK
