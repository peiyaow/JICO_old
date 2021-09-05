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

# setwd("~/Documents/GitHub/continuum/")
current = getwd()
setwd("/nas/longleaf/home/peiyao/continuum/")
load("./data/ADNI2/rank.RData")
# source("./data/ADNI2/loaddata_mac.R")
source("./data/ADNI2/loaddata.R")
source("./function/jive_continuum.R")
source("./function/cv_multigroup.R")
source("./function/PLS.R")
setwd(current)

set.seed(myseed)

label.level = levels(label)
X.list = lapply(label.level, function(l) X[label == l,])
Y.list = lapply(label.level, function(l) matrix(Y[label == l]))

L = 12
a = seq(0.2, 0.75, length.out = L)
gam.list = a/(1-a)
gam.list
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

ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, 
                                 nfolds = 10, maxiter = 200, criteria = "min", 
                                 orthIndiv = F)

print("Best iter parameter is:")
print(do.call(c, ml.iter.best$parameter))

file.name = "rank_iter.csv"
write.table(t(c(myseed, do.call(c, ml.iter.best$parameter))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)


