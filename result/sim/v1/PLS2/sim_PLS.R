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
alpha1 = rep(.6, r1) #OLS: 0
alpha2 = rep(.6, r2) #OLS: 0 

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

a = seq(0, 1, length.out = L+1)
gam.list = a/(1-a)
gam.list[L+1] = 1e10

ml.100.list = list()
for (gam in gam.list){
  print(gam)
  ml.100.list = list.append(ml.100.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 1000, gam = gam, rankJ = 1, rankA = c(0, 0),
                                                           center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F))
}

ml.200.list = list()
for (gam in gam.list){
  print(gam)
  ml.200.list = list.append(ml.200.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 1000, gam = gam, rankJ = 2, rankA = c(0, 0),
                                                                   center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F))
}

ml.111.list = list()
for (gam in gam.list){
  ml.111.list = list.append(ml.111.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 1000, gam = gam, rankJ = 1, rankA = c(1, 1),
                                     center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F))
}

ml.011.list = list()
for (gam in gam.list){
  print(gam)
  ml.011.list = list.append(ml.011.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 1000, gam = gam, rankJ = 0, rankA = c(1, 1),
                                                           center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F))
}

ml.022.list = list()
for (gam in gam.list){
  ml.022.list = list.append(ml.022.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 1000, gam = gam, rankJ = 0, rankA = c(2, 2),
                                                           center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = F))
}

ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = F, intercept = F)

ml.pls = plsr(Y~X, validation = "CV", center = F, scale = F)
ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)

ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = F)
ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)

ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, 
                                                  standardize = F, intercept = F))

ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F, scale = F))
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))

ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F, scale = F))
ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))


MSE = list()
for (ml in ml.100.list){
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}

for (ml in ml.200.list){
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}

for (ml in ml.111.list){
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}

for (ml in ml.011.list){
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}

for (ml in ml.022.list){
  MSE = list.append(MSE, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}

ml = ml.ridge
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newx = X.test.list[[g]], s = ml.ridge$lambda.min) - Y.test.list[[g]])^2)))

ml = ml.pls
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = r)[,,1] - Y.test.list[[g]])^2)))

ml = ml.pcr
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = r)[,,1] - Y.test.list[[g]])^2)))

# group models
ml = ml.ridge.list
MSE = list.append(MSE, sapply(1:G, function(g)
  mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))

ml = ml.pls.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = r + r.list[[g]])[,,1]
)^2)))

ml = ml.pcr.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = r + r.list[[g]])[,,1]
)^2)))

MSE = do.call(rbind, MSE)
# row.names(MSE) = c("iter.OLS", "iter.PLS", "iter.PCR", "global.ridge", "global.PLS", 
#                    "global.PCR", "group.ridge", "group.PLS", 
#                    "group.PCR")
file.name = "result.csv"
write.table(t(c(myseed, as.vector(MSE))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

