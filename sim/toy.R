library(pls)
library(caret)
library(glmnet)
library(methods)

source("~/Documents/GitHub/continuum/function/JICO.R")

G = 2
n1 = 50
n2 = 50
n = n1 + n2
p = 200
r = 1
r1 = 1
r2 = 1
r.list = list(r1, r2)
L = 50

alpha = rep(1, r)
alpha1 = rep(1, r1) 
alpha2 = rep(1, r2) 

X1 = mvrnorm(n1, rep(0, p), diag(p))
X2 = mvrnorm(n2, rep(0, p), diag(p))
X = rbind(X1, X2)

q = r
q1 = r1
q2 = r2
V = matrix(svd(X)$v[,1:q], ncol = q)%*%rep(1/sqrt(q), q)
V1 = matrix(svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:q1], ncol = q1)%*%rep(1/sqrt(q1), q1)
V2 = matrix(svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:q2], ncol = q2)%*%rep(1/sqrt(q2), q2)

e1 = rnorm(n1)*.2
Y1 = X1%*%V%*%alpha + X1%*%V1%*%alpha1 + e1

e2 = rnorm(n2)*.2
Y2 = X2%*%V%*%alpha + X2%*%V2%*%alpha2 + e2

Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

gam = 1e12

# fitting oracle
ml.JICO = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 300, gam = gam, rankJ = 1, rankA = c(1, 1),
                          center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F)

Yhat.list = lapply(1:G, function(g) as.numeric(ml.JICO$intercept[[g]])+ X.list[[g]]%*%ml.JICO$beta.C[[g]] + X.list[[g]]%*%ml.JICO$beta.Cind[[g]])






