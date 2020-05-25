library(MASS)
library(pls)
library(pracma)

r = 2
r1 = 2
r2 = 2
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 100 # 40
G = 2

# s = 10 #6
# beta = rep(1/5, p)
# beta1 = rep(c(1/100, 0, -1/100, 0), p/4)*s
# beta2 = -rep(c(1/100, 0, -1/100, 0), p/4)*s

# Q = randortho(p, type = "orthonormal")
# s = 0.5
# alpha = rep(1, r)
# alpha1 = c(0, 1)*s
# alpha2 = c(1, 0)*s
# beta = Q[,1:r]%*%alpha
# beta1 = Q[,(r+1):(r+r1)]%*%alpha1
# beta2 = Q[,(r+r1+1):(r+r1+r2)]%*%alpha2

# beta2 = rep(c(0, -1/100, 0, 1/100), p/4)*s
# beta1 = rep(c(1/100, -1/100), p/2)*s
# beta2 = -rep(c(1/100, -1/100), p/2)*s

s = 0.5
beta = c(rep(1/25, 50), rep(0, p-50))
beta1 = c(rep(1/50, 25), rep(0, p-25))*s
beta2 = c(rep(-1/50, 25), rep(0, p-25))*s

MSE.list = list()
for (ii in 1:10){
  # training data
  S1 = mvrnorm(n1, rep(0, r), diag(r))
  S2 = mvrnorm(n2, rep(0, r), diag(r))
  
  U = mvrnorm(r, rep(0, p), diag(p))
  
  T1 = mvrnorm(n1, rep(0, r1), diag(r1))
  T2 = mvrnorm(n2, rep(0, r2), diag(r2))
  
  U1 = mvrnorm(r1, rep(0, p), diag(p))
  U2 = mvrnorm(r2, rep(0, p), diag(p))
  
  E1 = mvrnorm(n1, rep(0, p), diag(p))*0.1
  E2 = mvrnorm(n2, rep(0, p), diag(p))*0.1
  
  X1 = S1%*%U + T1%*%U1 + E1
  X2 = S2%*%U + T2%*%U2 + E2
  
  e1 = rnorm(n1)*0.1
  e2 = rnorm(n2)*0.1
  
  Y1 = (S1%*%U + T1%*%U1)%*%beta + (T1%*%U1)%*%beta1 + e1
  Y2 = (S2%*%U + T2%*%U2)%*%beta + (T2%*%U2)%*%beta2 + e2
  
  X = rbind(X1, X2)
  Y = rbind(Y1, Y2)
  
  X.list = list(X1, X2)
  Y.list = list(Y1, Y2)
  
  
  # testing data
  S1 = mvrnorm(n1, rep(0, r), diag(r))
  S2 = mvrnorm(n2, rep(0, r), diag(r))
  
  #U = mvrnorm(r, rep(0, p), diag(p))
  
  T1 = mvrnorm(n1, rep(0, r1), diag(r1))
  T2 = mvrnorm(n2, rep(0, r2), diag(r2))
  
  #U1 = mvrnorm(r1, rep(0, p), diag(p))
  #U2 = mvrnorm(r2, rep(0, p), diag(p))
  
  E1 = mvrnorm(n1, rep(0, p), diag(p))*0.1
  E2 = mvrnorm(n2, rep(0, p), diag(p))*0.1
  
  X1.test = S1%*%U + T1%*%U1 + E1
  X2.test = S2%*%U + T2%*%U2 + E2
  
  e1 = rnorm(n1)*0.1
  e2 = rnorm(n2)*0.1
  
  Y1.test = (S1%*%U + T1%*%U1)%*%beta + (T1%*%U1)%*%beta1 + e1
  Y2.test = (S2%*%U + T2%*%U2)%*%beta + (T2%*%U2)%*%beta2 + e2
  
  X.test.list = list(X1.test, X2.test)
  Y.test.list = list(Y1.test, Y2.test)
  
  # global method
  ml.pls = plsr(Y~X, validation = "CV", center = F)
  ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
  
  ml.pcr = pcr(Y~X, validation = "CV", center = F)
  ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
  
  ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = F, intercept = F)
  
  # group method
  ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F))
  ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))
  
  ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F))
  ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))
  
  ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, standardize = F, intercept = F))
  
  # my method
  gam.list = c(log(c(seq(exp(0.75), exp(1), length.out = 4), seq(exp(1), exp(2), length.out = 4)[-1])), 1e10)
  ml.list = list()
  
  for (gam in gam.list){
    print(gam)
    ml.list = list.append(ml.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = gam, rankJ = r, rankA = c(r1, r2),
                                                             center.X = F, scale.X = F, center.Y = F, scale.Y = F))
  }
  
  # testing
  MSE = list()
  for (ml in ml.list){
    MSE = list.append(MSE, sapply(1:G, function(g) 
      mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  }
  
  # comparing methods
  # global
  ml = ml.pls
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = r)[,,1] - Y.test.list[[g]])^2)))
  
  ml = ml.pcr
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = r)[,,1] - Y.test.list[[g]])^2)))
  
  ml = ml.ridge
  MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newx = X.test.list[[g]], s = ml.ridge$lambda.min) - Y.test.list[[g]])^2)))
  
  # group
  # pls
  ml = ml.pls.list
  MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = r + r.list[[g]])[,,1]
  )^2)))
  
  # pcr
  ml = ml.pcr.list
  MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = r + r.list[[g]])[,,1]
  )^2)))
  
  # ridge
  ml = ml.ridge.list
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))
  
  MSE = do.call(rbind, MSE)
  row.names(MSE) = c(paste0("C_", round(gam.list, 2)[1:3]), "C_pls", paste0("C_", round(gam.list, 2)[5:7]), "C_pcr", 
                     "global_pls", "global_pcr", "global_ridge",
                     "group_pls", "group_pcr", "group_ridge")
  MSE.list = list.append(MSE.list, MSE)
  # t(t(apply(MSE, 1, mean)))
}

C_min = mean(sapply(MSE.list, function(MSE) min(apply(MSE[1:8,], 1, mean))))
MSE.res = Reduce("+", MSE.list)/10

rbind(C_min,t(t(c(apply(MSE.res, 1, mean)))))
