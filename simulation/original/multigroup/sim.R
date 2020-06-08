library(pracma)
library(MASS)

path = "~/Documents/GitHub/continuum/"
setwd(path)
source("./function/jive_continuum.R")
source("./function/cv_multigroup.R")

r = 2
r1 = 1
r2 = 1
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 200
G = 2
LL = 10

s = 10
t = 1
t1 = 1 #.5
t2 = 1 #.5

beta = c(rep(1/25, 50), rep(0, p-50))*s
beta1 = c(rep(0, p-25), rep(1/50, 25))*s
beta2 = c(rep(0, p-50), rep(1/50, 25), rep(0, 25))*s

RESULT = list()
RANK = list()
for (ii in 1:10){
  set.seed(myseeds[ii])
  S1 = mvrnorm(n1, rep(0, r), diag(r)*t)
  S2 = mvrnorm(n2, rep(0, r), diag(r)*t)
  U = mvrnorm(r, rep(0, p), diag(p))#%*%V%*%t(V)
  
  T1 = mvrnorm(n1, rep(0, r1), diag(r1)*t1)
  T2 = mvrnorm(n2, rep(0, r2), diag(r2)*t2)
  U1 = mvrnorm(r1, rep(0, p), diag(p))#%*%V1%*%t(V1)
  U2 = mvrnorm(r2, rep(0, p), diag(p))#%*%V2%*%t(V2)
  
  E1 = mvrnorm(n1, rep(0, p), diag(p))*.2
  E2 = mvrnorm(n2, rep(0, p), diag(p))*.2
  
  X1 = S1%*%U + T1%*%U1 + E1
  X2 = S2%*%U + T2%*%U2 + E2
  
  sum(E1^2)/sum(X1^2)
  sum(E2^2)/sum(X2^2)
  
  e1 = rnorm(n1)*1
  e2 = rnorm(n2)*1
  
  Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
  Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2
  
  sum(e1^2)/sum(Y1^2)
  sum(e2^2)/sum(Y2^2)
  
  while(!(sum(e1^2)/sum(Y1^2) <= 0.1)*(sum(e2^2)/sum(Y2^2) <= 0.1)){
    e1 = rnorm(n1)*1
    e2 = rnorm(n2)*1
    
    Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
    Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2
  }
  
  X = rbind(X1, X2)
  Y = rbind(Y1, Y2)
  
  X.list = list(X1, X2)
  Y.list = list(Y1, Y2)
  
  S1 = mvrnorm(n1, rep(0, r), diag(r)*t)
  S2 = mvrnorm(n2, rep(0, r), diag(r)*t)
  T1 = mvrnorm(n1, rep(0, r1), diag(r1)*t1)
  T2 = mvrnorm(n2, rep(0, r2), diag(r2)*t2)
  
  E1 = mvrnorm(n1, rep(0, p), diag(p))*.2
  E2 = mvrnorm(n2, rep(0, p), diag(p))*.2
  
  X1 = S1%*%U + T1%*%U1 + E1
  X2 = S2%*%U + T2%*%U2 + E2
  
  e1 = rnorm(n1)*1
  e2 = rnorm(n2)*1
  
  Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
  Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2
  
  while(!(sum(e1^2)/sum(Y1^2) <= 0.1)*(sum(e2^2)/sum(Y2^2) <= 0.1)){
    e1 = rnorm(n1)*1
    e2 = rnorm(n2)*1
    
    Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
    Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2
  }
  
  X.test = rbind(X1, X2)
  Y.test = rbind(Y1, Y2)
  
  X.test.list = list(X1, X2)
  Y.test.list = list(Y1, Y2)
  
  # train models
  ml.pls = plsr(Y~X, validation = "CV", center = F, scale = T)
  ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
  
  ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = T)
  ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
  
  # ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = T, intercept = F, nlambda = L-1)
  
  ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F, scale = T))
  ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))
  
  ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F, scale = T))
  ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))
  
  # ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, nlambda = L-1,
  #                                                  standardize = T, intercept = F))
  
  gam.list = c(exp(seq(log(0.5), log(1.5), length.out = LL-2)), 1e10)
  parameter.set = list()
  # for (gam in gam.list){
  #   ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m = 10, 
  #                                 center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", plot = T)
  #   rankJ = ml.step1$rankJ
  #   rankJ
  #   ml.step2 = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = gam, nfolds = 10, m = 10, 
  #                                                   center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", plot = T))
  #   rankA = sapply(ml.step2, function(ml) ml$rankA)
  #   rankA
  #   rankA = sapply(ml.step2, function(ml) max(ml$rankA - rankJ, 0))
  #   parameter = list(gam = gam, rankJ = rankJ, rankA = rankA)
  #   parameter.set = list.append(parameter.set, parameter)
  # }
  rMSE = list()
  for (gam in gam.list){
    ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m = 10,
                                  center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", plot = T)
    rMSE = list.append(rMSE, ml.step1$rMSE)
  }
  rMSE = do.call(rbind, rMSE)
  
  plot(1-rMSE[,4]/rMSE[,1], col = 4, ylim = c(0, 1), type = "l")
  points(1-rMSE[,2]/rMSE[,1], col = 2, type = "l")
  points(1-rMSE[,3]/rMSE[,1], col = 3, type = "l")
  points(1-rMSE[,5]/rMSE[,1], col = 5, type = "l")
  points(1-rMSE[,6]/rMSE[,1], col = 6, type = "l")
  
  rMSE = list()
  for (gam in gam.list){
    ml.step2 = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = gam, nfolds = 10, m = 10,
                                                    center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "1se", 
                                                    plot = T))
    
#    rMSE = list.append(rMSE, lapply(1:G, function(g) ml.step2[[g]]$rMSE))
    rMSE = list.append(rMSE, ml.step2[[1]]$rMSE)
  }
  rMSE = do.call(rbind, rMSE)
  
  plot(1-rMSE[,4]/rMSE[,1], col = 2, ylim = c(0, 1), type = "l")
  points(1-rMSE[,2]/rMSE[,1], type = "l")
  points(1-rMSE[,3]/rMSE[,1], col = 3, type = "l")
  
  points(1-rMSE[,5]/rMSE[,1], col = 4, type = "l")
  points(1-rMSE[,6]/rMSE[,1], col = 5, type = "l")
  
  #   rankA = sapply(ml.step2, function(ml) ml$rankA)
  #   rankA
  #   rankJ.max = min(rankA)
  #   for (rankJ in 0:rankJ.max){
  #     parameter = list(gam = gam, rankJ = rankJ, rankA = rankA - rankJ)
  #     parameter.set = list.append(parameter.set, parameter)
  #   }
  # }
  
  ml.2step.best = cv.continnum.2step(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, 
                                     center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min")
  ml.2step.list = list()
  for (parameter in parameter.set){
    print(parameter)
    ml.2step.list = list.append(ml.2step.list, 
                                continuum.2step(X.list, Y.list, lambda = 0, 
                                                gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                                center.X = F, scale.X = T, center.Y = F, scale.Y = T))
  }
  
  # ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, 
  #                                  center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", orthIndiv = F)
  # ml.iter.list = list()
  # for (parameter in parameter.set){
  #   print(parameter)
  #   ml.iter.list = list.append(ml.iter.list, 
  #                              continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
  #                                                        gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
  #                                                        center.X = F, scale.X = T, center.Y = F, scale.Y = T, orthIndiv = F))
  # }
  # 
  # ml.iter.orthIndiv.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, 
  #                                            center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", orthIndiv = T)
  # ml.iter.orthIndiv.list = list()
  # for (parameter in parameter.set){
  #   print(parameter)
  #   ml.iter.orthIndiv.list = list.append(ml.iter.orthIndiv.list, 
  #                                        continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
  #                                                                  gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
  #                                                                  center.X = F, scale.X = T, center.Y = F, scale.Y = T, orthIndiv = T))
  # }
  
  
  # testing
  MSE = list()
  
  MSE.2step = list()
  for (ml in ml.2step.list){
    MSE.2step = list.append(MSE.2step, sapply(1:G, function(g) 
      mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  }
#  MSE.2step = list.append(MSE.2step, MSE.2step[[ml.2step.best$ix]])
  MSE = list.append(MSE, MSE.2step[[ml.2step.best$ix]])
  
  # MSE.iter = list()
  # for (ml in ml.iter.list){
  #   MSE.iter = list.append(MSE.iter, sapply(1:G, function(g) 
  #     mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  # }
  # MSE.iter = list.append(MSE.iter, MSE.iter[[ml.iter.best$ix]])
  # MSE[LL+(1:LL)] = MSE.iter
  # 
  # MSE.orthIndiv.iter = list()
  # for (ml in ml.iter.orthIndiv.list){
  #   MSE.orthIndiv.iter = list.append(MSE.orthIndiv.iter, sapply(1:G, function(g) 
  #     mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  # }
  # MSE.orthIndiv.iter = list.append(MSE.orthIndiv.iter, MSE.orthIndiv.iter[[ml.iter.orthIndiv.best$ix]])
  # MSE[2*LL+(1:LL)] = MSE.orthIndiv.iter
  
  # global models
  ml = ml.pls
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pls)[,,1] - Y.test.list[[g]])^2)))
  
  ml = ml.pcr
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pcr)[,,1] - Y.test.list[[g]])^2)))
  
  # ml = ml.ridge
  # MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newx = X.test.list[[g]], s = ml.ridge$lambda.min) - Y.test.list[[g]])^2)))
  
  # group models
  ml = ml.pls.list
  MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
  )^2)))
  
  ml = ml.pcr.list
  MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]])[,,1]
  )^2)))
  
  
  MSE = do.call(rbind, MSE)
  print(ii)
  print(MSE)
  print(do.call(c, ml.2step.best$parameter))
  RESULT = list.append(RESULT, MSE)
  RANK = list.append(RANK, do.call(c, ml.2step.best$parameter))
                     # rbind(do.call(c, ml.2step.best$parameter), 
                     #             do.call(c, ml.iter.best$parameter), 
                     #             do.call(c, ml.iter.orthIndiv.best$parameter)))
}

save(RESULT, RANK, file = "result_sim.RData")

apply(Reduce("+", RESULT)/10, 1, mean)

a = seq(0, .99, by = 0.01)
plot(a/(1-a), a)

a = 0.2

