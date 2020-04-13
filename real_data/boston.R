library(mlbench)
library(caret)
library(glmnet)
library(pls)
source("~/Documents/GitHub/continuum/function/function.R")

data(BostonHousing2)
num = c("lon", "lat", "crim", "zn", "indus", "chas", "nox", "rm", "age", "dis", "tax", "ptratio", "b", "lstat")  
cat = c("town", "chas", "rad")
X.num = BostonHousing2[num]
X.cat = apply(BostonHousing2[cat], 2, as.factor)

Y = BostonHousing2["cmedv"]
X = apply(as.matrix(X.num[, -(1:2)]), 2, as.numeric)
label = X.cat[, 3]
summary(X.cat)

X = X[label %in% c("24", " 5", " 4"),]
Y = as.matrix(Y[label %in% c("24", " 5", " 4"),])
label = label[label %in% c("24", " 5", " 4")]

X = X[Y<=46,]
label = label[Y<=46]
Y = as.matrix(Y[Y<=46,])

n = dim(X)[1]
p = dim(X)[2]

ix.train = unlist(createDataPartition(label, times = 1, p = 3/4))
ix.test = (1:n)[-ix.train]
ix.list = list(ix.train, ix.test)
Y.list = lapply(1:2, function(x) Y[ix.list[[x]],]) # train,test
X.list = lapply(1:2, function(x) X[ix.list[[x]],])
label.list = lapply(1:2, function(x) label[ix.list[[x]]])

X.train = X.list[[1]]
X.test = X.list[[2]]
Y.train = Y.list[[1]]
Y.test = Y.list[[2]]
label.train = label.list[[1]]
label.test = label.list[[2]]

label.level = unique(label)[c(2,1,3)] # sort label
G = length(label.level)

X.train.list = lapply(label.level, function(l) X.train[label.train == l,])
X.test.list = lapply(label.level, function(l) X.test[label.test == l,])
Y.train.list = lapply(label.level, function(l) Y.train[label.train == l])
Y.test.list = lapply(label.level, function(l) Y.test[label.test == l])
label.train.list = lapply(label.level, function(l) label.train[label.train == l])
label.test.list = lapply(label.level, function(l) label.test[label.test == l])
n.train.vec = sapply(X.train.list, nrow)
n.test.vec = sapply(X.test.list, nrow)
n.vec = as.vector(table(label))
n.train = sum(n.train.vec)
n.test = sum(n.test.vec)

X.train = do.call(rbind, X.train.list)
X.test = do.call(rbind, X.test.list)
Y.train = do.call(c, Y.train.list)
Y.test = do.call(c, Y.test.list)
label.train = do.call(c, label.train.list)
label.test = do.call(c, label.test.list)
                     
# group-wise standardize X 
X.train.mean = lapply(X.train.list, colMeans)
X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))

X.train.list = lapply(1:G, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
X.train.list = lapply(1:G, function(ix) sweep(X.train.list[[ix]], 2, X.train.sd[[ix]], FUN = "/"))
X.test.list = lapply(1:G, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
X.test.list = lapply(1:G, function(ix) sweep(X.test.list[[ix]], 2, X.train.sd[[ix]], FUN = "/"))
X.train.list = lapply(1:G, function(g) matrix(X.train.list[[g]][!apply(X.train.list[[g]], 2, is.na)], nrow = n.train.vec[g]))
X.test.list = lapply(1:G, function(g) matrix(X.test.list[[g]][!apply(X.test.list[[g]], 2, is.na)], nrow = n.test.vec[g]))

# group-wise subtract Y mean
Y.train.mean = sapply(Y.train.list, mean)
Y.train.list = lapply(1:G, function(ix) Y.train.list[[ix]] - Y.train.mean[ix])

# global standardize X
X.train.mean = apply(X.train, 2, mean)
X.train.sd = apply(X.train, 2, sd)
X.train = sweep(X.train, 2, X.train.mean)
X.train = sweep(X.train, 2, X.train.sd, FUN = "/")

X.test = sweep(X.test, 2, X.train.mean)
X.test = sweep(X.test, 2, X.train.sd, FUN = "/")

ml.pls = plsr(Y.train ~ X.train, validation = "CV")
# ncomp.pls = selectNcomp(ml.pls, method = "onesigma", plot = F)
ncomp.pls = which.min(RMSEP(ml.pls)$val[1,,])-1

ml.pcr = pcr(Y.train~X.train, validation = "CV")
# ncomp.pcr = selectNcomp(ml.pcr, method = "onesigma", plot = F)
ncomp.pcr = which.min(RMSEP(ml.pcr)$val[1,,])-1

ml.ridge = cv.glmnet(x = X.train, y = Y.train, alpha = 0, lambda = exp(seq(log(1000), log(.1), length.out = 100)))
C.list = lapply(c(0, 0.2, 0.5, 1, 100), function(g) cv.continuum.exact(X.train, Y.train, gamma = g))

if (ncomp.pls){
  Yhat.pls = predict(ml.pls, newdata = X.test, ncomp = ncomp.pls)[,,1]
}else{
  Yhat.pls = rep(ml.pls$Ymeans, n*G)
}
if (ncomp.pcr){
  Yhat.pcr = predict(ml.pcr, newdata = X.test, ncomp = ncomp.pcr)[,,1]
}else{
  Yhat.pcr = rep(ml.pcr$Ymeans, n*G)
}
Yhat.ridge = predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.min)
# Yhat.ridge = predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.1se)

bhat.list = lapply(C.list, function(C) lm(Y.train~X.train%*%C)$coefficients)
Yhat.C.list = lapply(1:length(bhat.list), function(l) X.test%*%C.list[[l]]%*%bhat.list[[l]][-1]+bhat.list[[l]][1])


# mse.ols = mean((Yhat.ols - Y.test)^2)
mse.pls = mean((Yhat.pls - Y.test)^2)
mse.pcr = mean((Yhat.pcr - Y.test)^2)
mse.ridge = mean((Yhat.ridge - Y.test)^2)
mse.C = sapply(Yhat.C.list, function(Yhat) mean((Yhat - Y.test)^2))

print(c(mse.pls, mse.pcr, mse.ridge, mse.C))
# ----------------------------------------------------------------------

# -------------------------- group method ------------------------------
ml.pls.list = lapply(1:G, function(g) plsr(Y.train.list[[g]] ~ X.train.list[[g]], validation = "CV"))
# ncomp.pls.list = lapply(1:G, function(g) which.min(RMSEP(ml.pls.list[[g]])$val[1,,])-1)
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "onesigma", plot = F))

ml.pcr.list = lapply(1:G, function(g) pcr(Y.train.list[[g]]~X.train.list[[g]], validation = "CV"))
# ncomp.pcr.list = lapply(1:G, function(g) which.min(RMSEP(ml.pcr.list[[g]])$val[1,,])-1)
ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "onesigma", plot = F))


ml.C.list = lapply(1:G, function(g) lapply(c(0, 0.5, 1, 100), function(gam) cv.continuum.exact(X.train.list[[g]], Y.train.list[[g]], gamma = gam)))

ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.train.list[[g]], y = Y.train.list[[g]], alpha = 0, lambda = 
                                                    exp(seq(log(1000), log(.1), length.out = 100))))

Yhat.pls.list = lapply(1:G, function(g) if (ncomp.pls.list[[g]]){
  predict(ml.pls.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
} else{
  rep(0, n.test.vec[g])
})
Yhat.pcr.list = lapply(1:G, function(g) if (ncomp.pcr.list[[g]]){
  predict(ml.pcr.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]])[,,1]
} else{
  rep(0, n.test.vec[g])
})
Yhat.ridge.list = lapply(1:G, function(g) predict(ml.ridge.list[[g]], newx = X.test.list[[g]], s = ml.ridge.list[[g]]$lambda.min))

bhat.list.list = lapply(1:G, function(g) lapply(1:length(ml.C.list[[g]]), function(ix) lm(Y.train.list[[g]]~X.train.list[[g]]%*%ml.C.list[[g]][[ix]])$coefficients))
Yhat.C.list.list = lapply(1:G, function(g) lapply(1:length(bhat.list.list[[g]]), function(ix) X.test.list[[g]]%*%ml.C.list[[g]][[ix]]%*%bhat.list.list[[g]][[ix]][-1]+bhat.list.list[[g]][[ix]][1]))

mse.pls.list = lapply(1:G, function(g) mean((Yhat.pls.list[[g]]+Y.train.mean[g] - Y.test.list[[g]])^2))
mse.pcr.list = lapply(1:G, function(g) mean((Yhat.pcr.list[[g]]+Y.train.mean[g] - Y.test.list[[g]])^2))
mse.ridge.list = lapply(1:G, function(g) mean((Yhat.ridge.list[[g]]+Y.train.mean[g] - Y.test.list[[g]])^2))
mse.C.list = lapply(1:G, function(g) sapply(1:length(Yhat.C.list.list[[g]]), function(ix) mean((Yhat.C.list.list[[g]][[ix]]+Y.train.mean[g] - Y.test.list[[g]])^2)))

mse.group = rbind(do.call(c, mse.pls.list), do.call(c, mse.pcr.list), do.call(c, mse.ridge.list), do.call(cbind, mse.C.list))

mse.pls.group = sum(unlist(mse.pls.list)*n.test.vec)/n.test
mse.pcr.group = sum(unlist(mse.pcr.list)*n.test.vec)/n.test
mse.ridge.group = sum(unlist(mse.ridge.list)*n.test.vec)/n.test

print(mse.group)
# ----------------------------------------------------------------------

# pls int
cum.n.train.vec = c(0, cumsum(n.train.vec))
Yres.pls.list = lapply(1:G, function(g) if(ncomp.pls){
  ml.pls$residuals[(cum.n.train.vec[g]+1):cum.n.train.vec[g+1],, ncomp.pls]
}else{
  Y.train[label.train == label.level[g]]
})

Yres.mean.pls.list = lapply(1:G, function(g) mean(Yres.pls.list[[g]]))
Yres.pls.list = lapply(1:G, function(g) Yres.pls.list[[g]] - Yres.mean.pls.list[[g]])

ml.res.pls.list = lapply(1:G, function(g) plsr(Yres.pls.list[[g]] ~ X.train.list[[g]], validation = "CV"))
# ncomp.res.pls.list = lapply(1:G, function(g) which.min(RMSEP(ml.res.pls.list[[g]])$val[1,,])-1)
ncomp.res.pls.list = lapply(1:G, function(g) selectNcomp(ml.res.pls.list[[g]], method = "onesigma", plot = F))

Yhat.res.pls.list = lapply(1:G, function(g) if (ncomp.res.pls.list[[g]]){
  predict(ml.res.pls.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.res.pls.list[[g]])[,,1]
} else{
  rep(0, n.test.vec[g])
})
cum.n.test.vec = c(0, cumsum(n.test.vec))
Yhat.pls.int.list = lapply(1:G, function(g) Yhat.pls[(cum.n.test.vec[g]+1):cum.n.test.vec[g+1]]+Yhat.res.pls.list[[g]])
mse.pls.int.list = lapply(1:G, function(g) mean((Yhat.pls.int.list[[g]]+Yres.mean.pls.list[[g]] - Y.test.list[[g]])^2))
mse.pls.int = sum(unlist(mse.pls.int.list)*n.test.vec)/n.test

# ridge int
Yres.ridge.list = lapply(1:G, function(g) Y.train[label.train == label.level[g]] - predict(ml.ridge, s = ml.ridge$lambda.min, newx = X.train[label.train == label.level[g],]))
ml.res.ridge.list = lapply(1:G, function(g) cv.glmnet(y = Yres.ridge.list[[g]], x= X.train.list[[g]]))
Yhat.res.ridge.list = lapply(1:G, function(g) predict(ml.res.ridge.list[[g]], newx = X.test.list[[g]], s = ml.res.ridge.list[[g]]$lambda.min))
Yhat.ridge.int.list = lapply(1:G, function(g) Yhat.ridge[(cum.n.test.vec[g]+1):cum.n.test.vec[g+1]] + Yhat.res.ridge.list[[g]])
mse.ridge.int.list = lapply(1:G, function(g) mean((Yhat.ridge.int.list[[g]] - Y.test.list[[g]])^2))
mse.ridge.int = sum(unlist(mse.ridge.int.list)*n.test.vec)/n.test

print(c(mse.pls, mse.pcr, mse.ridge))
print(c(mse.pls.group, mse.pcr.group, mse.ridge.group))
print(c(mse.pls.int, mse.ridge.int))

