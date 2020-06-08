library(MASS)
library(pls)
library(pracma)
library(glmnet)

r = 2
r1 = 1
r2 = 1
n1 = 50
n2 = 50
n = n1 + n2
p = 200
G = 2
L = (c(12, 5))*10
L1 = (1.5)*10
L2 = (1.5)*10
alpha = c(1, 1)
alpha1 = rep(1, r1)
alpha2 = rep(1, r2)
LL = 10

RESULT = list()
RANK = list()
for (ii in 1:10){
  set.seed(myseeds[ii])
  Q = randortho(p, type = "orthonormal")
  U = t(matrix(Q[,1:r], ncol = r))
  U1 = t(matrix(Q[, r + (1:r1) ], ncol = r1))
  U2 = t(matrix(Q[, r+r1 + (1:r2) ], ncol = r2))
  
  # training data
  # S1 = mvrnorm(n, rep(0, r), diag(L))
  # S2 = mvrnorm(n, rep(0, r), diag(L))
  # T1 = mvrnorm(n, rep(0, r1), DIAG(L1))
  # T2 = mvrnorm(n, rep(0, r2), DIAG(L2))
  
  Q1 = matrix(randortho(n1, type = "orthonormal")[,1:r], ncol = r)
  Q2 = matrix(randortho(n2, type = "orthonormal")[,1:r], ncol = r)
  V1 = matrix(randortho(n1, type = "orthonormal")[,1:r1], ncol = r1)
  V2 = matrix(randortho(n2, type = "orthonormal")[,1:r2], ncol = r2)
  
  # S = V%*%diag(L)
  # S1 = S[1:n1,]
  # S2 = S[-(1:n1),]
  
  S1 = Q1%*%DIAG(L)
  S2 = Q2%*%DIAG(L)
  T1 = V1%*%DIAG(L1)
  T2 = V2%*%DIAG(L2)
  
  E1 = mvrnorm(n1, rep(0, p), diag(p))*.1
  E2 = mvrnorm(n2, rep(0, p), diag(p))*.1
  
  X1 = S1%*%U + T1%*%U1 + E1
  X2 = S2%*%U + T2%*%U2 + E2
  
  sum(E1^2)/sum(X1^2)
  sum(E2^2)/sum(X2^2)
  
  e1 = rnorm(n1)*2
  e2 = rnorm(n2)*2
  
  Y1 = S1%*%alpha + T1%*%alpha1 + e1
  Y2 = S2%*%alpha + T2%*%alpha2 + e2
  
  sum(e1^2)/sum(Y1^2)
  sum(e2^2)/sum(Y2^2)
  
  X = rbind(X1, X2)
  Y = rbind(Y1, Y2)
  
  X.list = list(X1, X2)
  Y.list = list(Y1, Y2)
  
  # testing data
  Q1 = matrix(randortho(n1, type = "orthonormal")[,1:r], ncol = r)
  Q2 = matrix(randortho(n2, type = "orthonormal")[,1:r], ncol = r)
  V1 = matrix(randortho(n1, type = "orthonormal")[,1:r1], ncol = r1)
  V2 = matrix(randortho(n2, type = "orthonormal")[,1:r2], ncol = r2)
  
  # S = V%*%diag(L)
  # S1 = S[1:n1,]
  # S2 = S[-(1:n1),]
  
  S1 = Q1%*%DIAG(L)
  S2 = Q2%*%DIAG(L)
  T1 = V1%*%DIAG(L1)
  T2 = V2%*%DIAG(L2)
  
  E1 = mvrnorm(n1, rep(0, p), diag(p))*.1 #.01 .1
  E2 = mvrnorm(n2, rep(0, p), diag(p))*.1 #.01 .1
  
  X1 = S1%*%U + T1%*%U1 + E1
  X2 = S2%*%U + T2%*%U2 + E2
  
  e1 = rnorm(n1)*2 #.1
  e2 = rnorm(n2)*2 #.1
  
  Y1 = S1%*%alpha + T1%*%alpha1 + e1
  Y2 = S2%*%alpha + T2%*%alpha2 + e2
  
  X.test = rbind(X1, X2)
  Y.test = rbind(Y1, Y2)
  
  X.test.list = list(X1, X2)
  Y.test.list = list(Y1, Y2)
  
  # train model
  ml.pls = plsr(Y~X, validation = "CV", center = F, scale = T)
  ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
  
  ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = T)
  ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
  
  # ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = T, intercept = T, nlambda = LL-1)
  
  ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F, scale = T))
  ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))
  
  ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F, scale = T))
  ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))
  
  # ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, nlambda = LL-1,
  #                                                  standardize = T, intercept = T))
  
  gam.list = c(exp(seq(log(0.5), log(1.5), length.out = LL-2)), 1e10)
  parameter.set = list()
  for (gam in gam.list){
    ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m = 5, 
                                  center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "scree")
    rankJ = ml.step1$rankJ
    ml.step2 = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = gam, nfolds = 10, m = 5, 
                                                    center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "1se"))
    rankA = sapply(ml.step2, function(ml) max(ml$rankA - rankJ, 0))
    parameter = list(gam = gam, rankJ = rankJ, rankA = rankA)
    parameter.set = list.append(parameter.set, parameter)
  }
  
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
  
  ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, 
                                   center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", orthIndiv = F)
  ml.iter.list = list()
  for (parameter in parameter.set){
    print(parameter)
    ml.iter.list = list.append(ml.iter.list, 
                               continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                                         gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                                         center.X = F, scale.X = T, center.Y = F, scale.Y = T, orthIndiv = F))
  }
  
  ml.iter.orthIndiv.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, 
                                             center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "min", orthIndiv = T)
  ml.iter.orthIndiv.list = list()
  for (parameter in parameter.set){
    print(parameter)
    ml.iter.orthIndiv.list = list.append(ml.iter.orthIndiv.list, 
                                         continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                                                   gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                                                   center.X = F, scale.X = T, center.Y = F, scale.Y = T, orthIndiv = T))
  }
  
  # testing
  MSE = list()
  
  MSE.2step = list()
  for (ml in ml.2step.list){
    MSE.2step = list.append(MSE.2step, sapply(1:G, function(g) 
      mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  }
  MSE.2step = list.append(MSE.2step, MSE.2step[[ml.2step.best$ix]])
  MSE[1:LL] = MSE.2step
  
  MSE.iter = list()
  for (ml in ml.iter.list){
    MSE.iter = list.append(MSE.iter, sapply(1:G, function(g) 
      mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  }
  MSE.iter = list.append(MSE.iter, MSE.iter[[ml.iter.best$ix]])
  MSE[LL+(1:LL)] = MSE.iter
  
  MSE.orthIndiv.iter = list()
  for (ml in ml.iter.orthIndiv.list){
    MSE.orthIndiv.iter = list.append(MSE.orthIndiv.iter, sapply(1:G, function(g) 
      mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  }
  MSE.orthIndiv.iter = list.append(MSE.orthIndiv.iter, MSE.orthIndiv.iter[[ml.iter.orthIndiv.best$ix]])
  MSE[2*LL+(1:LL)] = MSE.orthIndiv.iter
  
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
  
  # ml = ml.ridge.list
  # MSE = list.append(MSE, sapply(1:G, function(g) 
  #   mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))
  
  MSE = do.call(rbind, MSE)
  print(ii)
  print(MSE)
  RESULT = list.append(RESULT, MSE)
  RANK = list.append(RANK, rbind(do.call(c, ml.2step.best$parameter), 
                                 do.call(c, ml.iter.best$parameter), 
                                 do.call(c, ml.iter.orthIndiv.best$parameter)))
}

save(RESULT, RANK, file = "result_PCR_tune.RData")
