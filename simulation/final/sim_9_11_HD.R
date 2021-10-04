# ---------------------- reading shell command --------------------- 
# args = (commandArgs(TRUE))
# cat(args, "\n")
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# ------------------------------------------------------------------ 
myseed_I = 1
myseed_X = 1234

method = "PCR2"
source('~/GitHub/continuum/function/jive_continuum.R')
# setwd("D:/git-project/multi-source/contiuum from peiyao/simulation/final")
# source('~/function/jive_continuum.R')
library(MASS)
library(rlist)
library(glmnet)
library(methods)
library(pls)
library(doParallel)
# current = getwd()
# setwd("/nas/longleaf/home/peiyao/continuum/")
# setwd(current)
set.seed(myseed_X)
# cl = makeCluster(10)
# registerDoParallel(cl)
# getDoParWorkers()

L = 50
a = seq(0, 1, length.out = L+1)
gam.list = a/(1-a)
gam.list[L+1] = 1e10
CT_NULL = list()
CT_random = list()
CT_far = list()
nrun_NULL = rep(NULL,length(gam.list))
nrun_random = rep(NULL,length(gam.list))
nrun_far = rep(NULL,length(gam.list))
CT_homo_NULL = list()
CT_homo_random = list()
CT_homo_far = list()
CT_heter_NULL = list()
CT_heter_random = list()
CT_heter_far = list()
MSE_NULL = rep(NULL,length(gam.list))
MSE_random = rep(NULL,length(gam.list))
MSE_far = rep(NULL,length(gam.list))
G = 2
n1 = 50
n2 = 50
n = n1 + n2
p = 200
r = 1
r1 = 1
r2 = 1
r.list = list(r1, r2)
L = 10

alpha = rep(1, r)
alpha1 = rep(1, r1) #OLS: 0
alpha2 = rep(1, r2) #OLS: 0 


if (method == "PCR1"){
  X1 = mvrnorm(n1, rep(0, p), diag(p))
  X2 = mvrnorm(n2, rep(0, p), diag(p))
  X = rbind(X1, X2)
  q = r
  q1 = r1
  q2 = r2
  V = matrix(svd(X)$v[,1:r], ncol = r)
  V1 = matrix(svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:r1], ncol = r1)
  V2 = matrix(svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:r2], ncol = r2)
  
  
  X.list = list(X1, X2)
  
  X1 = mvrnorm(n1, rep(0, p), diag(p))
  X2 = mvrnorm(n2, rep(0, p), diag(p))
  
  X.test = rbind(X1, X2)
  
  X.test.list = list(X1, X2)
  
}
if (method == "PLS"){
  X1 = mvrnorm(n1, rep(0, p), diag(p))
  X2 = mvrnorm(n2, rep(0, p), diag(p))
  X = rbind(X1, X2)
  
  q = min(n, p)/2
  q1 = min(n1, p)/2
  q2 = min(n2, p)/2
  V = matrix(svd(X)$v[,1:q], ncol = q)%*%rep(1/sqrt(q), q)
  V1 = matrix(svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:q1], ncol = q1)%*%rep(1/sqrt(q1), q1)
  V2 = matrix(svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:q2], ncol = q2)%*%rep(1/sqrt(q2), q2)
  
  X.list = list(X1, X2)
  X1 = mvrnorm(n1, rep(0, p), diag(p))
  X2 = mvrnorm(n2, rep(0, p), diag(p))
  
  X.test = rbind(X1, X2)
  
  X.test.list = list(X1, X2)
  
}
if (method == "PCR2"){
  
  library(pracma)
  
  Q = randortho(n1)
  U = matrix(Q[,1:r], ncol = r)
  U1 = matrix(Q[,r+(1:r1)], ncol = r1)
  Q = randortho(n2)
  U = rbind(U,matrix(Q[,1:r], ncol = r))
  U2 = matrix(Q[,r+(1:r2)], ncol = r2)
  Q = randortho(p)
  V = matrix(Q[,1:r], ncol = r)
  V1 = matrix(Q[,r+(1:r1)], ncol = r1)
  V2 = matrix(Q[,r+r1+(1:r2)], ncol = r2)
  D = (rnorm(r,1,0))
  D1 = (rnorm(r1,1,0))
  D2 = (rnorm(r2,1,0))
  if (r>1)
    D = diag(D)
  if (r1>1)
    D1 = diag(D1)
  if (r2>1)
    D2 = diag(D2)
  J = U%*%D%*%t(V)
  J1 = J[1:n1,]
  J2 = J[(n1+1):n,]
  
  I1 = U1%*%(D1)%*%t(V1)
  I2 = U2%*%(D2)%*%t(V2)
  E1 = mvrnorm(n1,rep(0,p),diag(0.4,p))
  E2 = mvrnorm(n2,rep(0,p),diag(0.4,p))
  X1=J1+I1
  X2=J2+I2
  X = rbind(X1,X2)
  
  X.list = list(X1, X2)
  
  
  Q = randortho(n1)
  U = matrix(Q[,1:r], ncol = r)
  U1 = matrix(Q[,r+(1:r1)], ncol = r1)
  Q = randortho(n2)
  U = rbind(U,matrix(Q[,1:r], ncol = r))
  U2 = matrix(Q[,r+(1:r2)], ncol = r2)
  # D = (rnorm(r,1,0.1))
  # D1 = (rnorm(r1,1,0.4))
  # D2 = (rnorm(r2,1,0.4))
  if (r>1)
    D = diag(D)
  if (r1>1)
    D1 = diag(D1)
  if (r2>1)
    D2 = diag(D2)
  J = U%*%D%*%t(V)
  J1 = J[1:n1,]
  J2 = J[(n1+1):n,]
  
  I1 = U1%*%(D1)%*%t(V1)
  I2 = U2%*%(D2)%*%t(V2)
  E1 = mvrnorm(n1,rep(0,p),diag(0.4,p))
  E2 = mvrnorm(n2,rep(0,p),diag(0.4,p))
  X1=J1+I1
  X2=J2+I2
  
  X.test = rbind(X1, X2)
  
  X.test.list = list(X1, X2)
}

a = seq(0, 1, length.out = L+1)
gam.list = a/(1-a)
gam.list[L+1] = 1e10
# print(i)

e1 = rnorm(n1)*.2
Y1 = X.list[[1]]%*%V%*%alpha + X.list[[1]]%*%V1%*%alpha1 + e1

e2 = rnorm(n2)*.2
Y2 = X.list[[2]]%*%V%*%alpha + X.list[[2]]%*%V2%*%alpha2 + e2

Y = rbind(Y1, Y2)
Y.list = list(Y1, Y2)

e1 = rnorm(n1)*.2
Y1 = X.test.list[[1]]%*%V%*%alpha + X.test.list[[1]]%*%V1%*%alpha1 + e1

e2 = rnorm(n2)*.2
Y2 = X.test.list[[2]]%*%V%*%alpha + X.test.list[[2]]%*%V2%*%alpha2 + e2

# X.test = rbind(X1, X2)
Y.test = rbind(Y1, Y2)

Y.test.list = list(Y1, Y2)

rep=1
for(i in 1:rep)
{
  set.seed(myseed_I)
  MSE_temp = rep(NULL,length(gam.list))
  j=1
  set.seed(i*j*123)
  I.random = list(matrix(rnorm(n1*p),nrow = n1,ncol = p),matrix(rnorm(n2*p),nrow = n2,ncol = p))
  
  ml.111.list = list()
  for (gam in gam.list){
    # print(gam)
    ml.111.list = list.append(ml.111.list, continuum.multigroup.iter(X.list, Y.list, lambda = 1e-10, maxiter = 300, gam = gam, rankJ = 1, rankA = c(1, 1),center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv =F, I.initial = I.random))
  }
  nrun_random = rbind(nrun_random,unlist(lapply(ml.111.list, function(x) x$nrun)))
  CT_random = list.append(CT_random,lapply(ml.111.list, function(x) x$CT.list))
  CT_homo_random = list.append(CT_homo_random,lapply(ml.111.list, function(x) x$CThomo.list))
  CT_heter_random = list.append(CT_heter_random,lapply(ml.111.list, function(x) x$CTheter.list))
  
  MSE = list()
  for (ml in ml.111.list){
    MSE = list.append(MSE, sapply(1:G, function(g)
      mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
  }
  MSE_random =  rbind(MSE_random,rowSums(do.call(rbind, MSE)))
}

# save(MSE_random, file = "MSE_initial_PCR_111_rep300_9_10_new_1e9.RData")
plot(x = a, y = MSE_random[1,],type = "l",ylab = "MSE", xlab = expression(paste("a = ", gamma, "/(1+", gamma, ")")), main = "Overall MSE of JICO with different initial value")
# for(k in 2:10){
#   lines(x=a, y= MSE_random[k,])
# }

# lines(x=a, y= MSE_NULL,col = "red")

# 
# 
# ml.111.list =  foreach(j=1:length(gam.list), .packages=c('MASS', 'rlist', 'foreach')) %dopar% {
#   source('~/GitHub/continuum/function/jive_continuum.R')
#   continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 300, gam = gam.list[j], rankJ = 1, rankA = c(1, 1),
#                             center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F)
# }  
# nrun_NULL = rbind(nrun_NULL,unlist(lapply(ml.111.list, function(x) x$nrun)))
# CT_NULL = list.append(CT_NULL,lapply(ml.111.list, function(x) x$CT.list))
# CT_homo_NULL = list.append(CT_homo_NULL,lapply(ml.111.list, function(x) x$CThomo.list))
# CT_heter_NULL = list.append(CT_heter_NULL,lapply(ml.111.list, function(x) x$CTheter.list))
# 
# MSE = list()
# for (ml in ml.111.list){
#   MSE = list.append(MSE, sapply(1:G, function(g) 
#     mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
# }
# MSE_NULL =  rbind(MSE_NULL,rowSums(do.call(rbind, MSE)))
# 
# plot(x = a, y = MSE_random[1,],type = "l",ylab = "MSE", xlab = expression(paste("a = ", gamma, "/(1+", gamma, ")")), main = "Overall MSE of JICO with different initial value")
# for(k in 2:30){
#   lines(x=a, y= MSE_random[k,])
# }
# 
# lines(x=a, y= MSE_NULL,col = "red")
# 
# 
# source("D:/git-project/multi-source/contiuum from peiyao/function/jive_continuum.R")
# MSE_random_PCR = NULL
# for (i in 1:30){
#   print(i)
#   MSE_temp = rep(NULL,length(gam.list))
#   MSE = 1e10
#   for (j in 1:10){
#     set.seed(i*j*123)
#     I.random = list(matrix(rnorm(n1*p),nrow = n1,ncol = p),matrix(rnorm(n2*p),nrow = n2,ncol = p))
#     
#     ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 300, gam = 1e10, rankJ = 1, rankA = c(1, 1),
#                                    center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F,I.initial = I.random)
#     
#     MSE = min(MSE,sapply(1:G, function(g) 
#       mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
#   }
#   MSE_random_PCR =  c(MSE_random_PCR,MSE)
# }
# 
# ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 300, gam = 1e10, rankJ = 1, rankA = c(1, 1),
#                                center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F)
# 
# MSE = min(MSE,sapply(1:G, function(g) 
#   mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
