# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
cat(args, "\n")
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------ 

library(MASS)
library(rlist)
library(glmnet)
library(methods)
library(pls)

current = getwd()
setwd("/nas/longleaf/home/peiyao/continuum/")
source("./function/jive_continuum.R")
setwd("/nas/longleaf/home/peiyao/alpha/")
source("./function/main_function.R")
setwd(current)
set.seed(myseed)

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
alpha1 = rep(0, r1) #OLS: 0
alpha2 = rep(0, r2) #OLS: 0 

# Q = randortho(p)
# V = matrix(Q[,1:r], ncol = r)
# V1 = matrix(Q[,r+(1:r1)], ncol = r1)
# V2 = matrix(Q[,r+r1+(1:r2)], ncol = r2)

X1 = mvrnorm(n1, rep(0, p), diag(p))
X2 = mvrnorm(n2, rep(0, p), diag(p))
X = rbind(X1, X2)

q = min(n, p)
q1 = min(n1, p)
q2 = min(n2, p)
V = svd(X)$v[,1:q]%*%rep(1/sqrt(q), q)
V1 = svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:q1]%*%rep(1/sqrt(q1), q1)
V2 = svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:q2]%*%rep(1/sqrt(q2), q2)

e1 = rnorm(n1)*.2
Y1 = X1%*%V%*%alpha + X1%*%V1%*%alpha1 + e1

e2 = rnorm(n2)*.2
Y2 = X2%*%V%*%alpha + X2%*%V2%*%alpha2 + e2

sum(e1^2)/sum(Y1^2)
sum(e2^2)/sum(Y2^2)

Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

X1 = mvrnorm(n1, rep(0, p), diag(p))
X2 = mvrnorm(n2, rep(0, p), diag(p))

e1 = rnorm(n1)*.2
Y1 = X1%*%V%*%alpha + X1%*%V1%*%alpha1 + e1

e2 = rnorm(n2)*.2
Y2 = X2%*%V%*%alpha + X2%*%V2%*%alpha2 + e2

X.test = rbind(X1, X2)
Y.test = rbind(Y1, Y2)

X.test.list = list(X1, X2)
Y.test.list = list(Y1, Y2)

# ------------------------------- ALPHA -----------------------------------------
X2U.list = lapply(1:G, function(ix) X2U1(X.list[[ix]], K = 10, plot = F))
H.list = lapply(X2U.list, function(list) list$H)
K.list = lapply(X2U.list, function(list) list$K)
P.list = lapply(X2U.list, function(list) list$P)
L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 

F.train.list = lapply(1:G, function(ix) X2U.list[[ix]]$F_)
U.train.list = lapply(1:G, function(ix) X2U.list[[ix]]$U)

FnU.test.list = lapply(1:G, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
F.test.list = lapply(FnU.test.list, function(list) list$F_) 
U.test.list = lapply(FnU.test.list, function(list) list$U)

# OLS.F
data.F.train.list = lapply(1:G, function(ix) data.frame(Y = Y.list[[ix]], 
                                                        F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:G, function(l) lm(Y~., data = data.F.train.list[[l]]))

# OLS.U
U.train = do.call(rbind, U.train.list)
HY.train.list = lapply(1:G, function(ix) H.list[[ix]]%*%Y.list[[ix]])
HY.train = do.call(c, HY.train.list)

ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0, standardize = T)
EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5, standardize = T)
lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = T)

HYhat.test.ridge.OLS.list = lapply(1:G, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.EN.OLS.list = lapply(1:G, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.lasso.OLS.list = lapply(1:G, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
PYhat.test.list = lapply(1:G, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)

# Yhat.test.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
Yhat.test.ridge.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
Yhat.test.EN.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
Yhat.test.lasso.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])

mse.ridge.OLS.list = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
mse.EN.OLS.list = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
mse.lasso.OLS.list = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)

#------------------------------------------ALPHA-0----------------------------------------------------
X2U.list = lapply(1:G, function(ix) X2U2(X.list[[ix]], K = 0, plot = F))
H.list = lapply(X2U.list, function(list) list$H)
K.list = lapply(X2U.list, function(list) list$K)
P.list = lapply(X2U.list, function(list) list$P)
L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 

F.train.list = lapply(1:G, function(ix) X2U.list[[ix]]$F_)
U.train.list = lapply(1:G, function(ix) X2U.list[[ix]]$U)

FnU.test.list = lapply(1:G, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
F.test.list = lapply(FnU.test.list, function(list) list$F_) 
U.test.list = lapply(FnU.test.list, function(list) list$U) 

# OLS.F
data.F.train.list = lapply(1:G, function(ix) data.frame(Y = Y.list[[ix]], 
                                                        F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:G, function(l) lm(Y~., data = data.F.train.list[[l]]))

# OLS.U
U.train = do.call(rbind, U.train.list)
HY.train.list = lapply(1:G, function(ix) H.list[[ix]]%*%Y.list[[ix]])
HY.train = do.call(c, HY.train.list)

ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0, standardize = T)
EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5, standardize = T)
lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = T)

# HYhat.test.OLS.list = lapply(1:G, function(ix) U.test.list[[ix]]%*%beta.OLS.U)
HYhat.test.ridge.OLS.list = lapply(1:G, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.EN.OLS.list = lapply(1:G, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.lasso.OLS.list = lapply(1:G, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
PYhat.test.list = lapply(1:G, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)

# Yhat.test.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
Yhat.test.ridge.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
Yhat.test.EN.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
Yhat.test.lasso.OLS.list = lapply(1:G, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])

# mse.OLS.list = compute.mse(Y.test.list, Yhat.test.OLS.list)
mse.ridge.OLS.list0 = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
mse.EN.OLS.list0 = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
mse.lasso.OLS.list0 = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)

MSE = rbind(mse.ridge.OLS.list$mse.vec, mse.EN.OLS.list$mse.vec, mse.lasso.OLS.list$mse.vec,
            mse.ridge.OLS.list0$mse.vec, mse.EN.OLS.list0$mse.vec, mse.lasso.OLS.list0$mse.vec)

file.name = "result.csv"
write.table(t(c(myseed, as.vector(MSE))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)
