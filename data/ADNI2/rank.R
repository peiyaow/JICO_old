library(RcppAlgos)
# install.packages("RcppAlgos")
hehe = list(list(c(0,0,0)), list(c(1,0,0)), list(c(2,0,0), c(1,1,0)), list(c(3,0,0), c(2,1,0), c(1,1,1)),
            list(c(4,0,0), c(3,1,0), c(2,2,0), c(2,1,1)))

RANK = list()
for (s in 0:1){ # all possible rank up to sum = 1..4
  for (rankJ in 0:s){
    for (cc in hehe[[s-rankJ+1]]){
      RANK = list.append(RANK, t(apply(unique(permuteGeneral(cc)), 1, function(x) c(rankJ, x))))
    }
  }
}
RANK = do.call(rbind, RANK)

# save(RANK, file = "rank1.RData")
nrow(RANK)
parameter.set = list()
for (i in 1:nrow(RANK)){
  for (gam in gam.list){
    rankJ = RANK[i,1]
    rankA = RANK[i,-1]
    parameter = list(gam = gam, rankJ = rankJ, rankA = rankA)
    parameter.set = list.append(parameter.set, parameter)
  }
}
