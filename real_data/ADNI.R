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

# source("~/Documents/GitHub/continuum/data/ADNI2/loaddata.R")
# source("~/Documents/GitHub/continuum/function/jive_continuum.R")
# source("~/Documents/GitHub/continuum/function/cv_multigroup.R")

current = getwd()
setwd("/nas/longleaf/home/peiyao/continuum/")
source("./data/ADNI2/loaddata.R")
load("./data/ADNI2/rank2.RData")
source("./function/jive_continuum.R")
source("./function/cv_multigroup.R")
source("./function/PLS.R")
setwd(current)

set.seed(myseed)

n = dim(X)[1]
p = dim(X)[2]
L = 50 # number of gam

idtrain = unlist(createDataPartition(label, times = 1, p = 4/5))
idtest = (1:n)[-idtrain]

X.train = X[idtrain,]
X.test = X[-idtrain,]

Y.train = matrix(Y[idtrain])
Y.test = matrix(Y[-idtrain])

label.train = matrix(label[idtrain])
label.test = matrix(label[-idtrain])

label.level = levels(label)
G = length(label.level)

X.list = lapply(label.level, function(l) X.train[label.train == l,])
X.test.list = lapply(label.level, function(l) X.test[label.test == l,])
Y.list = lapply(label.level, function(l) matrix(Y.train[label.train == l,]))
Y.test.list = lapply(label.level, function(l) matrix(Y.test[label.test == l,]))
label.list = lapply(label.level, function(l) matrix(label.train[label.train == l,]))
label.test.list = lapply(label.level, function(l) matrix(label.test[label.test == l,]))

X = do.call(rbind, X.list)
Y = do.call(rbind, Y.list)
label = do.call(rbind, label.list)

# -------------------------------- train models -------------------------------- 
ml.pls = plsr(Y~X, validation = "CV", center = T, scale = T)
ncomp.pls = selectNcomp(ml.pls, method = "onesigma", plot = F)

ml.pcr = pcr(Y~X, validation = "CV", center = T, scale = T)
ncomp.pcr = selectNcomp(ml.pcr, method = "onesigma", plot = F)

ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = T, intercept = T)

ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = T, scale = T))
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "onesigma", plot = F))

ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = T, scale = T))
ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "onesigma", plot = F))

ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0,
                                                  standardize = T, intercept = T))

# my models
# parameters
a = seq(0, 1, length.out = L+1)
gam.list = a/(1-a)
gam.list[L+1] = 1e10

parameter.set = list()
for (i in 1:nrow(RANK)){
  for (gam in gam.list){
    rankJ = RANK[i,1]
    rankA = RANK[i,-1]
    parameter = list(gam = gam, rankJ = rankJ, rankA = rankA)
    parameter.set = list.append(parameter.set, parameter)
  }
}

# tune best 2step model
ml.2step.best = cv.continnum.2step(X.list, Y.list, lambda = 0, parameter.set, 
                                   # center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                   nfolds = 10, criteria = "min")
ml.2step.list = list()
for (parameter in parameter.set){
  print(parameter)
  ml.2step.list = list.append(ml.2step.list, 
                              continuum.2step(X.list, Y.list, lambda = 0, 
                                              # center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                              gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA))
}

# tune best iterate model not orthogonal
ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, criteria = "min", orthIndiv = F)
ml.iter.list = list()
for (parameter in parameter.set){
  print(parameter)
  ml.iter.list = list.append(ml.iter.list, 
                             continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                                       gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                                       orthIndiv = F))
}

# -------------------------------- testing --------------------------------
MSE.2step = list()
for (ml in ml.2step.list){
  MSE.2step = list.append(MSE.2step, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}
MSE.2step = list.append(MSE.2step, MSE.2step[[ml.2step.best$ix]])

file.name = "rank_2step.csv"
write.table(t(c(myseed, do.call(c, ml.2step.best$parameter))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_2step.csv"
write.table(t(c(myseed, as.vector(do.call(rbind, MSE.2step)))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

MSE.iter = list()
for (ml in ml.iter.list){
  MSE.iter = list.append(MSE.iter, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}
MSE.iter = list.append(MSE.iter, MSE.iter[[ml.iter.best$ix]])
# MSE[(L+1)*nrow(RANK)+1:((L+1)*nrow(RANK))] = MSE.iter

file.name = "rank_iter.csv"
write.table(t(c(myseed, do.call(c, ml.iter.best$parameter))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_iter.csv"
write.table(t(c(myseed, as.vector(do.call(rbind, MSE.iter)))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

MSE = list()
MSE = list.append(MSE, MSE.2step[[ml.2step.best$ix]])
MSE = list.append(MSE, MSE.iter[[ml.iter.best$ix]])

# global models
ml = ml.pls
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict.wrapper(ml, X.test.list[[g]], ncomp.pls) - Y.test.list[[g]])^2)))

ml = ml.pcr
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict.wrapper(ml, X.test.list[[g]], ncomp.pcr) - Y.test.list[[g]])^2)))

ml = ml.ridge
MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newx = X.test.list[[g]], s = ml.ridge$lambda.min) - Y.test.list[[g]])^2)))

# group models
ml = ml.pls.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict.wrapper(ml[[g]], X.test.list[[g]], ncomp.pls.list[[g]]))^2)))

ml = ml.pcr.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict.wrapper(ml[[g]], X.test.list[[g]], ncomp.pcr.list[[g]]))^2)))

ml = ml.ridge.list
MSE = list.append(MSE, sapply(1:G, function(g) 
  mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))

print(do.call(c, ml.2step.best$parameter))
print(do.call(c, ml.iter.best$parameter))
print(do.call(rbind, MSE))

file.name = "result.csv"
write.table(t(c(myseed, as.vector(do.call(rbind, MSE)))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)


