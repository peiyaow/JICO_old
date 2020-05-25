library(MASS)
library(pls)
library(glmnet)
source("~/Documents/GitHub/continuum/function/jive.R")
source("~/Documents/GitHub/continuum/function/jive_continuum.R")

r = 2
r1 = 1
r2 = 1
n = 50
p = 100
G = 2
beta = rep(1/25, p)
beta1 = rep(-1/50, p)
beta2 = rep(1/50, p)

#MSE = list()
#flag = c()
#for (k in 1:50){
  # S1 = mvrnorm(n, rep(0, r), diag(r))
  # S2 = matrix(rbinom(n*r, 1, 1/2), n, r)
  # U = mvrnorm(r, rep(0, p), diag(p))
  # T1 = matrix(rbinom(n*r1, 1, 1/2), nrow = n, ncol = r1)
  # T2 = matrix(runif(n*r2), nrow = n, ncol = r2)
  # W1 = matrix(runif(r1*p), nrow = r1, ncol = p)
  # W2 = matrix(rbinom(r2*p, 1, 1/2), nrow = r2, ncol = p)
  
  
  S1 = mvrnorm(n, rep(0, r), diag(r))
  S2 = mvrnorm(n, rep(0, r), diag(r))
  U = mvrnorm(r, rep(0, p), diag(p))
  T1 = mvrnorm(n, rep(0, r1), diag(r1))
  T2 = mvrnorm(n, rep(0, r2), diag(r2))
  W1 = mvrnorm(r1, rep(0, p), diag(p))
  W2 = mvrnorm(r2, rep(0, p), diag(p))
#  E1 = mvrnorm(n, rep(0,p), diag(p))
#  E2 = mvrnorm(n, rep(0,p), diag(p))
  S1 = scale(S1)
  S2 = scale(S2)
  T1 = scale(T1)
  T2 = scale(T2)
  
  
  X1 = S1%*%U + T1%*%W1
  X2 = S2%*%U + T2%*%W2
  
  e1 = scale(rnorm(n))*2
  e2 = scale(rnorm(n))*2
#  Y1 = 1+X1%*%beta + e1
#  Y2 = 2+X2%*%beta + e2
  Y1 = 1+S1%*%U%*%beta + T1%*%W1%*%beta1 + e1
  Y2 = 2+S2%*%U%*%beta + T2%*%W2%*%beta2 + e2
  
  X = rbind(X1, X2)
  Y = rbind(Y1, Y2)
  
  X.list = list(X1, X2)
  Y.list = list(Y1, Y2)

  # testing data
  S1 = mvrnorm(n, rep(0, r), diag(r))
  S2 = mvrnorm(n, rep(0, r), diag(r))
  #U = mvrnorm(r, rep(0, p), diag(p))
  T1 = mvrnorm(n, rep(0, r1), diag(r1))
  T2 = mvrnorm(n, rep(0, r2), diag(r2))
  #W1 = mvrnorm(r1, rep(0, p), diag(p))
  #W2 = mvrnorm(r2, rep(0, p), diag(p))
  S1 = scale(S1)
  S2 = scale(S2)
  T1 = scale(T1)
  T2 = scale(T2)
  
  X1 = S1%*%U + T1%*%W1
  X2 = S2%*%U + T2%*%W2
  
  e1 = scale(rnorm(n))*2
  e2 = scale(rnorm(n))*2
#  Y1 = 1+X1%*%beta + e1
#  Y2 = 2+X2%*%beta + e2
  Y1 = 1+S1%*%U%*%beta + T1%*%W1%*%beta1 + e1
  Y2 = 2+S2%*%U%*%beta + T2%*%W2%*%beta2 + e2
  
  X.test = rbind(X1, X2)
  Y.test = rbind(Y1, Y2)
  
  X.test.list = list(X1, X2)
  Y.test.list = list(Y1, Y2)
  
  ml.jive = jive(X.list, rankJ = r, rankA = c(r1, r2), method = "given", orthIndiv = FALSE)
  
#  MSE[[k]] = list()
  
  # method 1 jive
  beta.jive = C2beta.jive(X.list, Y.list, ml.jive$U[[ml.jive$nrun]], ml.jive$W[[ml.jive$nrun]])
  Yhat.homo.jive.list = lapply(1:G, function(g) X.test.list[[g]]%*%beta.jive$beta.C[[g]])
  Yhat.heter.jive.list = lapply(1:G, function(g) X.test.list[[g]]%*%beta.jive$beta.Cind[[g]])
  MSE.intercept.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n))^2))
  MSE.homo.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n) - Yhat.homo.jive.list[[g]])^2))
  MSE.heter.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n) - Yhat.heter.jive.list[[g]])^2))
  MSE.jive = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(beta.jive$intercept[[g]], n) - Yhat.homo.jive.list[[g]] - Yhat.heter.jive.list[[g]])^2))
#  MSE[[k]][[1]] = rbind(MSE.intercept.jive, MSE.homo.jive, MSE.heter.jive, MSE.jive)
  
  # continuum
  ml.continuum.int = continuum.int.2step(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2))
#  flag[k] = ml.continuum.int$converged
  Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum.int$beta.C[[g]])
  Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum.int$beta.Cind[[g]])
  MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum.int$intercept[[g]], n))^2))
  MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum.int$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
  MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum.int$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
  MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum.int$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
  rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum)
#  MSE[[k]][[2]] = rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum)
#}


rbind(MSE.jive, MSE.continuum)

