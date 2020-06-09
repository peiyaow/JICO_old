L = 10
a = seq(0, 1, length.out = L+1)
gam.list = a/(1-a)
gam.list[L+1] = 1e10
parameter.set = list()
for (gam in gam.list){
  m = 5
  if (gam == 0){
    m = 1
  }
  ml.joint = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m = m,
                                # center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                criteria = "1se", plot = T)
  rankJ = ml.joint$rankJ
  print(rankJ)
  parameter = list(gam = gam, rankJ = rankJ, rankA = rep(0, G))
  parameter.set = list.append(parameter.set, parameter)
  
  ml.separate = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = gam, 
                                                     nfolds = 10, m = m,
                                                     # center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                                     criteria = "min", plot = T))
  rankA = sapply(ml.separate, function(ml) ml$rankA)
  print(rankA)
  rankJ.max = min(rankA)
  for (rankJ in 0:rankJ.max){
    parameter = list(gam = gam, rankJ = rankJ, rankA = rankA - rankJ)
    parameter.set = list.append(parameter.set, parameter)
  }
}
