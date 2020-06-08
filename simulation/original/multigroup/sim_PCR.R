library(MASS)
library(pls)
library(pracma)

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

RESULT = list()
RANK = list()
myseeds = floor(runif(50)*1e4)
for (ii in 1:50){
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
  ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = 1e10, nfolds = 10, m = 10, 
                                center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "scree", plot = T)
  rankJ = ml.step1$rankJ
  # rankJ
  ml.step2 = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = 1e10, nfolds = 10, m = 10, 
                                                  center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "1se", plot = T))
  rankA = sapply(ml.step2, function(ml) max(ml$rankA - rankJ, 0))
  # rankA
  
  ml.2step = continuum.2step(X.list, Y.list, lambda = 0, 
                             gam = 1e10, rankJ = rankJ, rankA = rankA, 
                             center.X = F, scale.X = T, center.Y = F, scale.Y = T)
  
  ml.iter = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200, gam = 1e10, rankJ = rankJ, rankA = rankA,
                                      center.X = F, scale.X = T, center.Y = F, scale.Y = T, orthIndiv = F)
  
  ml.iter.orthIndiv = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200, gam = 1e10, rankJ = rankJ, rankA = rankA,
                                                center.X = F, scale.X = T, center.Y = F, scale.Y = T, orthIndiv = T)
  
  ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = T)
  ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
  
  ml.pls = plsr(Y~X, validation = "CV", center = T, scale = T)
  ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
  
  ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F, scale = T))
  ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))
  
  ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = T, scale = T))
  ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))
  
  # testing
  MSE = list()
  ml = ml.2step
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  
  ml = ml.iter
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  
  ml = ml.iter.orthIndiv
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  
  ml = ml.pcr
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pcr)[,,1] - Y.test.list[[g]])^2)))
  
  ml = ml.pls
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pcr)[,,1] - Y.test.list[[g]])^2)))
  
  ml = ml.pcr.list
  MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]])[,,1]
  )^2)))
  
  ml = ml.pls.list
  MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
  )^2)))
  
  MSE = do.call(rbind, MSE)
  print(ii)
  print(MSE)
  RESULT = list.append(RESULT, MSE)
  RANK = list.append(RANK, c(rankJ, rankA))
}

save(RESULT, RANK, file = "result_PCR.RData")
save(myseeds, file = "seeds.RData")
