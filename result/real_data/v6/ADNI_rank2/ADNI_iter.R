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
load("./data/ADNI2/rank3.RData")
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

ml.iter.list = list()
for (parameter in parameter.set){
  print(parameter)
  ml.iter.list = list.append(ml.iter.list,
                             continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 500,
                                                       gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
                                                       orthIndiv = F))
}

# testing
MSE.iter = list()
for (ml in ml.iter.list){
  MSE.iter = list.append(MSE.iter, sapply(1:G, function(g)
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}

file.name = "result_iter.csv"
write.table(t(c(myseed, as.vector(do.call(rbind, MSE.iter)))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)






