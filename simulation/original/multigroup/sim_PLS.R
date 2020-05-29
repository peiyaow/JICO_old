# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
cat(args, "\n")
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------ 

library(pls)
library(caret)
library(glmnet)
library(methods)
library(pracma)
library(MASS)

current = getwd()
setwd("/nas/longleaf/home/peiyao/continuum/")
source("./data/ADNI2/loaddata.R")
source("./function/jive_continuum.R")
source("./function/cv_multigroup.R")
load("./simulation/original/multigroup/pls.RData")
setwd(current)

set.seed(myseed)

L = 50

# training data
Q1 = randortho(n1, type = "orthonormal")
Q2 = randortho(n2, type = "orthonormal")

J = U%*%D%*%t(V)
J1 = Q1%*%J[1:n1,]
J2 = Q2%*%J[n1+(1:n2),]
I1 = Q1%*%U1%*%D1%*%t(V1)
I2 = Q2%*%U2%*%D2%*%t(V2)

E1 = mvrnorm(n1, rep(0, p), diag(p))*.2
E2 = mvrnorm(n2, rep(0, p), diag(p))*.2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

sum(E1^2)/sum(X1^2)
sum(E2^2)/sum(X2^2)

e1 = rnorm(n1)*1
e2 = rnorm(n2)*1

Y1 = J1%*%beta + I1%*%beta1 + e1
Y2 = J2%*%beta + I2%*%beta2 + e2

sum(e1^2)/sum(Y1^2)
sum(e2^2)/sum(Y2^2)

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

# testing data
Q1 = randortho(n1, type = "orthonormal")
Q2 = randortho(n2, type = "orthonormal")

J = U%*%D%*%t(V)
J1 = Q1%*%J[1:n1,]
J2 = Q2%*%J[n1+(1:n2),]
I1 = Q1%*%U1%*%D1%*%t(V1)
I2 = Q2%*%U2%*%D2%*%t(V2)

E1 = mvrnorm(n1, rep(0, p), diag(p))*.2
E2 = mvrnorm(n2, rep(0, p), diag(p))*.2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e1 = rnorm(n1)*1
e2 = rnorm(n2)*1

Y1 = J1%*%beta + I1%*%beta1 + e1
Y2 = J2%*%beta + I2%*%beta2 + e2

X.test = rbind(X1, X2)
Y.test = rbind(Y1, Y2)

X.test.list = list(X1, X2)
Y.test.list = list(Y1, Y2)

# train models
ml.pls = plsr(Y~X, validation = "CV", center = F, scale = T)
ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)

ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = T)
ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)

ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = T, intercept = F, nlambda = L-1)

ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F, scale = T))
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))

ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F, scale = T))
ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))

ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, nlambda = L-1,
                                                  standardize = T, intercept = F))

gam.list = exp(seq(log(0.5), log(1.5), length.out = L-1))
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
                                   center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                   criteria = "min")
ml.2step.list = list()
for (parameter in parameter.set){
  print(parameter)
  ml.2step.list = list.append(ml.2step.list, 
                              continuum.2step(X.list, Y.list, lambda = 0, 
                                              gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                              center.X = F, scale.X = T, center.Y = F, scale.Y = T))
}

ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, 
                                 center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                 maxiter = 200, criteria = "min")
ml.iter.list = list()
for (parameter in parameter.set){
  print(parameter)
  ml.iter.list = list.append(ml.iter.list, 
                             continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                                       gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                                       center.X = F, scale.X = T, center.Y = F, scale.Y = T))
}


# testing 
MSE = list()
MSE.2step = list()
for (ml in ml.2step.list){
  MSE.2step = list.append(MSE.2step, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}
MSE.2step = list.append(MSE.2step, MSE.2step[[ml.2step.best$ix]])
MSE[1:L] = MSE.2step

MSE.iter = list()
for (ml in ml.iter.list){
  MSE.iter = list.append(MSE.iter, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}
MSE.iter = list.append(MSE.iter, MSE.iter[[ml.iter.best$ix]])
MSE[L+(1:L)] = MSE.iter

# global models
ml = ml.pls
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pls)[,,1] - Y.test.list[[g]])^2)))

ml = ml.pcr
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pcr)[,,1] - Y.test.list[[g]])^2)))

ml = ml.ridge
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newx = X.test.list[[g]], s = ml.ridge$lambda.min) - Y.test.list[[g]])^2)))

# group models
ml = ml.pls.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
)^2)))

ml = ml.pcr.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]])[,,1]
)^2)))

ml = ml.ridge.list
MSE = list.append(MSE, sapply(1:G, function(g) 
  mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))

MSE = do.call(rbind, MSE)
print(MSE)
RANK = t(sapply(parameter.set, function(parameter) c(parameter$rankJ, parameter$rankA)))

file.name = "result.csv"
write.table(t(c(myseed, as.vector(MSE))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "rank.csv"
write.table(t(c(myseed, as.vector(RANK))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)


