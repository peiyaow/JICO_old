library(MASS)
library(pls)
library(purrr)
library(pracma)

r = 2
r1 = 2
r2 = 2
n = 50
p = 100
G = 2
L = c(10, 10)
L1 = c(10, 10)
L2 = c(10, 10)
alpha = c(100, 50)
alpha1 = rep(10, r1)
alpha2 = rep(10, r2)

# beta = rep(...)
# beta1 = rep(1, p)
# beta2 = rep(1, p)

Q = randortho(p, type = "orthonormal")
DIFF = list()
CONV = list()
NRUN = list()
RESULT = list()
NCOMP = list()

for (ii in 1:10){
  S1 = mvrnorm(n, rep(0, r), diag(L))
#  S1 = matrix(rep(L, n), ncol = r)
  S2 = mvrnorm(n, rep(0, r), diag(L))
#  S2 = matrix(rep(L, n), ncol = r)
  #  S = rbind(S1, S2)
#  V = eigen(t(S)%*%S)$vectors
  
#  U = mvrnorm(r, rep(0, p), diag(p))
  U = Q[1:r,]

  T1 = mvrnorm(n, rep(0, r1), DIAG(L1))
#  T1 = matrix(rep(L1, n), ncol = r1)
  T2 = mvrnorm(n, rep(0, r2), DIAG(L2))
#  T2 = matrix(rep(L2, n), ncol = r2)
#  V1 = eigen(t(T1)%*%T1)$vectors
#  V2 = eigen(t(T2)%*%T2)$vectors

#  W1 = mvrnorm(r1, rep(0, p), diag(p))
  W1 = matrix(Q[(r+1):(r+r1),], ncol = p, byrow = F)
#  W2 = mvrnorm(r2, rep(0, p), diag(p))
  W2 = matrix(Q[(r+r1+1):(r+r1+r2),], ncol = p, byrow = F)

  E1 = mvrnorm(n, rep(0,p), diag(p)*.01)
  E2 = mvrnorm(n, rep(0,p), diag(p)*.01)
  
#  X1 = S1%*%V%*%U + T1%*%V1%*%W1 + E1
#  X2 = S2%*%V%*%U + T2%*%V2%*%W2 + E2
  
  X1 = S1%*%U + T1%*%W1 + E1
  X2 = S2%*%U + T2%*%W2 + E2
  
  e1 = (rnorm(n))
  e2 = (rnorm(n))
  # Y1 = S1%*%U%*%beta + T1%*%W1%*%beta1 + e1
  # Y2 = S2%*%U%*%beta + T2%*%W2%*%beta2 + e2
  # Y1 = S1%*%V%*%alpha + T1%*%V1%*%alpha1 + e1
  # Y2 = S2%*%V%*%alpha + T2%*%V2%*%alpha2 + e2

  Y1 = S1%*%alpha + T1%*%alpha1 + e1
  Y2 = S2%*%alpha + T2%*%alpha2 + e2
  
  # cor(S1, Y1)
  # cor(T1, Y1)
  # t(S1)%*%S1
  # t(S2)%*%S2
  # t(T1)%*%T1
  # eigen(t(T2)%*%T2)$values
  # eigen(t(E1)%*%E1)$values
  
  X = rbind(X1, X2)
  Y = rbind(Y1, Y2)
#  t(X%*%t(U))%*%(X)%*%t(U)
  cor(rbind(S1, S2), Y)
  cor(rbind(T1, T2), Y)
  cor(T1, Y1)
  cor(T2, Y2)
  
  X.list = list(X1, X2)
  Y.list = list(Y1, Y2)
  
  # testing data
  S1 = mvrnorm(n, rep(0, r), diag(L))
  S2 = mvrnorm(n, rep(0, r), diag(L))
#  S = rbind(S1, S2)
#  V = eigen(t(S)%*%S)$vectors
  T1 = mvrnorm(n, rep(0, r1), DIAG(L1))
  T2 = mvrnorm(n, rep(0, r2), DIAG(L2))
#  V1 = eigen(t(T1)%*%T1)$vectors
#  V2 = eigen(t(T2)%*%T2)$vectors
  E1 = mvrnorm(n, rep(0,p), diag(p)*.01)
  E2 = mvrnorm(n, rep(0,p), diag(p)*.01)
  
  X1 = S1%*%U + T1%*%W1 + E1
  X2 = S2%*%U + T2%*%W2 + E2
  
  e1 = (rnorm(n))
  e2 = (rnorm(n))
  # Y1 = S1%*%U%*%beta + T1%*%W1%*%beta1 + e1
  # Y2 = S2%*%U%*%beta + T2%*%W2%*%beta2 + e2
  Y1 = S1%*%alpha + T1%*%alpha1 + e1
  Y2 = S2%*%alpha + T2%*%alpha2 + e2
  
  X.test = rbind(X1, X2)
  Y.test = rbind(Y1, Y2)
  
  X.test.list = list(X1, X2)
  Y.test.list = list(Y1, Y2)
  
  result = list()
  conv = list()
  nrun = list()
  ncomp = list()
  # ml.continuum = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2), maxiter = 300)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  # 
  # ml.continuum = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1.2, rankJ = r, rankA = c(r1, r2), maxiter = 300)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  # 
  # ml.continuum = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2), maxiter = 300)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  # 
  # ml.continuum1 = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = .8, rankJ = r, rankA = c(r1, r2), maxiter = 300)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  # 
  # ml.continuum2 = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = .5, rankJ = r, rankA = c(r1, r2), maxiter = 300)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

  ml.jive = jive.multigroup(X.list, rankJ = r, rankA = c(r1, r2), method = "given", orthIndiv = FALSE)
  beta.jive = C2beta.jive(X.list, Y.list, ml.jive$U[[ml.jive$nrun]], ml.jive$W[[ml.jive$nrun]])
  Yhat.homo.jive.list = lapply(1:G, function(g) X.test.list[[g]]%*%beta.jive$beta.C[[g]])
  Yhat.heter.jive.list = lapply(1:G, function(g) X.test.list[[g]]%*%beta.jive$beta.Cind[[g]])
  MSE.intercept.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n))^2))
  MSE.homo.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n) - Yhat.homo.jive.list[[g]])^2))
  MSE.heter.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n) - Yhat.heter.jive.list[[g]])^2))
  MSE.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n) - Yhat.homo.jive.list[[g]] - Yhat.heter.jive.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.jive, MSE.homo.jive, MSE.heter.jive, MSE.jive))
  nrun = list.append(nrun, ml.jive$nrun)
  
  ml.continuum.pcr = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2), maxiter = 300, scale.X = T)
  ml.continuum = ml.continuum.pcr
  Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  conv = list.append(conv, ml.continuum$converged)
  nrun = list.append(nrun, ml.continuum$nrun)
  
  diff = mean((do.call(rbind, ml.continuum$J) - do.call(rbind, ml.jive$joint))^2) + mean((do.call(rbind, ml.continuum$I) - do.call(rbind, ml.jive$individual))^2)
  DIFF = list.append(DIFF, diff)
  
  ml.continuum = continuum.2step.v1(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2))
  Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  
  ml.pcr = mycvpcr(X.list, Y.list)
  Yhat.homo.pcr.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pcr$beta.C[[g]])
  Yhat.heter.pcr.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pcr$beta.Cind[[g]])
  MSE.intercept.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n))^2))
  MSE.homo.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n) - Yhat.homo.pcr.list[[g]])^2))
  MSE.heter.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n) - Yhat.heter.pcr.list[[g]])^2))
  MSE.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n) - Yhat.homo.pcr.list[[g]] - Yhat.heter.pcr.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.pcr, MSE.homo.pcr, MSE.heter.pcr, MSE.pcr))
  ncomp = list.append(ncomp, c(ml.pcr$rankA, ml.pcr$rankJ))
  # ml.continuum = continuum.2step.v2(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2))
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  
  ml.continuum.pls = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2), maxiter = 300, scale.X = T)
  ml.continuum = ml.continuum.pls
  Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  conv = list.append(conv, ml.continuum$converged)
  nrun = list.append(nrun, ml.continuum$nrun)
  CONV = list.append(CONV, do.call(c, conv))
  NRUN = list.append(NRUN, do.call(c,nrun))
  
  ml.continuum = continuum.2step.v1(X.list, Y.list, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2))
  Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  
  ml.pls = mycvplsr(X.list, Y.list)
  Yhat.homo.pls.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pls$beta.C[[g]])
  Yhat.heter.pls.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pls$beta.Cind[[g]])
  MSE.intercept.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n))^2))
  MSE.homo.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n) - Yhat.homo.pls.list[[g]])^2))
  MSE.heter.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n) - Yhat.heter.pls.list[[g]])^2))
  MSE.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n) - Yhat.homo.pls.list[[g]] - Yhat.heter.pls.list[[g]])^2))
  result = list.append(result, rbind(MSE.intercept.pls, MSE.homo.pls, MSE.heter.pls, MSE.pls))
  ncomp = list.append(ncomp, c(ml.pls$rankA, ml.pls$rankJ))
  
  # ml.continuum = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1, rankJ = ml.pls$rankJ, rankA = ml.pls$rankA)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  
  # ml.continuum = continuum.2step.v1(X.list, Y.list, lambda = 0, gam = 1, rankJ = ml.pls$rankJ, rankA = ml.pls$rankA)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  
  # ml.continuum = continuum.2step.v1(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = ml.pcr$rankJ, rankA = ml.pcr$rankA)
  # Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
  # Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
  # MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
  # MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  # MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  # MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  # result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))
  RESULT = list.append(RESULT, result)
  NCOMP = list.append(NCOMP, ncomp)
}

RESULT.overall = lapply(RESULT, function(list) t(sapply(list, function(X) X[4,])))
RESULT.overall = reduce(RESULT.overall, `+`)/10
row.names(RESULT.overall) = c("JIVE", "C.PCR", "C_2step.PCR", "PCR", "C.PLS", "C_2step.PLS", "PLS")
NCOMP

save.image("PLS.RData")
c