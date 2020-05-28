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
source("./function/jive_continuum.R")
source("./function/cv_multigroup.R")
setwd(current)

set.seed(myseed)

n = dim(X)[1]
p = dim(X)[2]

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

# X.train.mean = apply(X.train, 2, mean)
# X.train.sd = apply(X.train, 2, sd)
# X.train = sweep(X.train, 2, X.train.mean)
# X.train = sweep(X.train, 2, X.train.sd, "/")
# X.test = sweep(X.test, 2, X.train.mean)
# X.test = sweep(X.test, 2, X.train.sd, "/")

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
ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)

ml.pcr = pcr(Y~X, validation = "CV", center = T, scale = T)
ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)

ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = T, intercept = T)

ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = T, scale = T))
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))

ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = T, scale = T))
ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))

ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, 
                                                  standardize = T, intercept = T))

# my models

# c(log(c(seq(exp(0.5), exp(1), length.out = 4), seq(exp(1), exp(1.5), length.out = 4)[-1])))
L = 20
gam.list = c(log(c(seq(exp(0), exp(1), length.out = L/2), seq(exp(1), exp(2), length.out = L/2)[-1])))
#c(log(c(seq(exp(0.2), exp(1), length.out = 7))))

# greedy selection of rankJ and rankA using 2 separate procedure
ml.2step.list = lapply(gam.list, function(gam) cv.continnum.2step.separate(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m1 = 10, m2 = 5))
parameter.set = lapply(ml.2step.list, function(ml) ml$parameter) # parameter set

# tune best 2step model
ml.2step.best = cv.continnum.2step(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, criteria = "min")

# tune best iterate model
ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, criteria = "min")

ml.iter.list = list()
for (parameter in parameter.set){
  print(parameter)
  ml.iter.list = list.append(ml.iter.list, 
                        continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
                                                  gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA))
}

# -------------------------------- testing --------------------------------
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
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]]+1)[,,1]
)^2)))

ml = ml.pcr.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]]+1)[,,1]
)^2)))

ml = ml.ridge.list
MSE = list.append(MSE, sapply(1:G, function(g) 
  mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))

print(do.call(rbind, MSE))

file.name = "result.csv"
write.table(t(c(myseed, as.vector(do.call(rbind, MSE)))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)
# do.call(rbind, MSE)%*%sapply(X.test.list, function(X) nrow(X))/sum(sapply(X.test.list, function(X) nrow(X)))



