library(pls)
library(glmnet)
library(MASS)

n = 50
p = 100
G = 6

ksi = .1
b.list = list()
b.list[[1]] = matrix(rnorm(p, mean = 0, sd = 1))
for (j in 2:G){
  b.list[[j]] = diag(rbinom(p, size = 1, prob = ksi)*2-1)%*%b.list[[j-1]]
}
b = matrix(rep(.2, p))

s = seq(1, 3, length.out = 6)

# ------------------------- generate data ------------------------------
X.list = lapply(1:G, function(g) generateX(n, p, s[g]))
eps.list = lapply(1:G, function(l) scale(rnorm(n)))
Y.list = lapply(1:G, function(l) X.list[[l]]%*%(b.list[[l]]+b) + eps.list[[l]])

X.test.list = lapply(1:G, function(l) generateX(n, p, s[g]))
eps.test.list = lapply(1:G, function(l) scale(rnorm(n)))
Y.test.list = lapply(1:G, function(l) X.test.list[[l]]%*%(b.list[[l]]+b) + eps.test.list[[l]])

# -------------------------- global method -----------------------------
X = do.call(rbind, X.list)
Y = do.call(c, Y.list)

ml.pls = plsr(Y ~ X, validation = "CV")
# ncomp.pls = selectNcomp(ml.pls, method = "onesigma", plot = F)
ncomp.pls = which.min(RMSEP(ml.pls)$val[1,,])-1
ml.pcr = pcr(Y~X, validation = "CV")
# ncomp.pcr = selectNcomp(ml.pcr, method = "onesigma", plot = F)
ncomp.pcr = which.min(RMSEP(ml.pcr)$val[1,,])-1
ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0)

X.test = do.call(rbind, X.test.list)
Y.test = do.call(c, Y.test.list)

# Yhat.ols = cbind(1, X.test)%*%ml.ols$coefficients
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

# mse.ols = mean((Yhat.ols - Y.test)^2)
mse.pls = mean((Yhat.pls - Y.test)^2)
mse.pcr = mean((Yhat.pcr - Y.test)^2)
mse.ridge = mean((Yhat.ridge - Y.test)^2)

# -------------------------- group method -----------------------------
# ml.ols.list = lapply(1:G, function(g) lm(Y.list[[g]]~X.list[[g]]))
ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV"))
ncomp.pls.list = lapply(1:G, function(g) which.min(RMSEP(ml.pls.list[[g]])$val[1,,])-1)
# ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "onesigma", plot = F))

ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV"))
ncomp.pcr.list = lapply(1:G, function(g) which.min(RMSEP(ml.pcr.list[[g]])$val[1,,])-1)
# ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "onesigma", plot = F))

ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0, lambda = exp(seq(log(1000), log(.1), length.out = 100))))

# Yhat.ols.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%ml.ols.list[[g]]$coefficients)
Yhat.pls.list = lapply(1:G, function(g) if (ncomp.pls.list[[g]]){
  predict(ml.pls.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
} else{
  rep(ml.pls.list[[g]]$Ymeans, n)
})
Yhat.pcr.list = lapply(1:G, function(g) if (ncomp.pcr.list[[g]]){
  predict(ml.pcr.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]])[,,1]
} else{
  rep(ml.pcr.list[[g]]$Ymeans, n)
})
Yhat.ridge.list = lapply(1:G, function(g) predict(ml.ridge.list[[g]], newx = X.test.list[[g]], s = ml.ridge.list[[g]]$lambda.min))

# mse.ols.list = lapply(1:G, function(g) mean((Yhat.ols.list[[g]] - Y.test.list[[g]])^2))
mse.pls.list = lapply(1:G, function(g) mean((Yhat.pls.list[[g]] - Y.test.list[[g]])^2))
mse.pcr.list = lapply(1:G, function(g) mean((Yhat.pcr.list[[g]] - Y.test.list[[g]])^2))
mse.ridge.list = lapply(1:G, function(g) mean((Yhat.ridge.list[[g]] - Y.test.list[[g]])^2))

# mse.ols.group = mean(do.call(c, mse.ols.list))
mse.pls.group = do.call(mean, mse.pls.list)
mse.pcr.group = do.call(mean, mse.pcr.list)
mse.ridge.group = do.call(mean, mse.ridge.list)

# -------------------------- intergrated method -----------------------------
# pls int
Yres.pls.list = lapply(1:G, function(g) if(ncomp.pls){
  ml.pls$residuals[1:n+(g-1)*n,,ncomp.pls]
}else{
  Y.list[[g]] - ml.pls$Ymeans
})

ml.res.pls.list = lapply(1:G, function(g) plsr(Yres.pls.list[[g]] ~ X.list[[g]], validation = "CV"))
ncomp.res.pls.list = lapply(1:G, function(g) which.min(RMSEP(ml.res.pls.list[[g]])$val[1,,])-1)
# ncomp.res.pls.list = lapply(1:G, function(g) selectNcomp(ml.res.pls.list[[g]], method = "onesigma", plot = F))
Yhat.res.pls.list = lapply(1:G, function(g) if (ncomp.res.pls.list[[g]]){
  predict(ml.res.pls.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.res.pls.list[[g]])[,,1]
} else{
  rep(ml.res.pls.list[[g]]$Ymeans, n)
})

Yhat.pls.int.list = lapply(1:G, function(g) Yhat.pls[1:n+(g-1)*n]+Yhat.res.pls.list[[g]])
mse.pls.int.list = lapply(1:G, function(g) mean((Yhat.pls.int.list[[g]] - Y.test.list[[g]])^2))
mse.pls.int = do.call(mean, mse.pls.int.list)

# pls maxmin
hatb.pls.list = lapply(1:G, function(g) if(ncomp.pls.list[[g]]){
  ml.pls.list[[g]]$coefficients[,,ncomp.pls.list[[g]]]
}else{
  rep(0,p)
})
hatb.pls = do.call(cbind, hatb.pls.list)
w.pls.maxmin = maxmin(X, hatb.pls)$solution
hatb.pls.maxmin = hatb.pls%*%w.pls.maxmin
Yhat.pls.maxmin.list = lapply(1:G, function(g) ml.pls.list[[g]]$Ymeans + X.test.list[[g]]%*%hatb.pls.maxmin)
# Yhat.pls.maxmin.list = lapply(1:G, function(g) ml.pls$Ymeans + X.test.list[[g]]%*%hatb.pls.maxmin)

mse.pls.maxmin.list = lapply(1:G, function(g) mean((Yhat.pls.maxmin.list[[g]] - Y.test.list[[g]])^2))
mse.pls.maxmin = do.call(mean, mse.pls.maxmin.list)

c(mse.pls, mse.pcr, mse.ridge)
c(mse.pls.group, mse.pcr.group, mse.ridge.group)
c(mse.pls, mse.pls.group, mse.pls.int, mse.pls.maxmin)



