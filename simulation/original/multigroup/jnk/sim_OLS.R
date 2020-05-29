library(pracma)
# Ctot = randortho(p, type = "orthonormal")
# C = matrix(Ctot[,1:r], ncol = r)
# Ctot = Ctot[,-(1:r)]
# C1 = matrix(Ctot[, 1:r1])
# Ctot = Ctot[,-(1:r1)]
# C2 = matrix(Ctot[, 1:r2])
# 
# beta = C%*%c(1, -1)
# beta1 = C1*0.5
# beta2 = C2*0.5

beta = rep(1/25, p)
beta1 = c(rep(0, 50), rep(-1/50, 25), rep(0, 25))
beta2 = c(rep(0, 50), rep(1/50, 25), rep(0, 25))

X.list = lapply(1:G, function(g) (mvrnorm(n, rep(0, p), diag(p))))
e1 = (rnorm(n))
e2 = (rnorm(n))
Y1 = X.list[[1]]%*%(beta+beta1) + e1
Y2 = X.list[[2]]%*%(beta+beta2) + e2
Y.list = list(Y1, Y2)

X.test.list = lapply(1:G, function(g) (mvrnorm(n, rep(0, p), diag(p))))
e1 = (rnorm(n))
e2 = (rnorm(n))
Y1 = X.test.list[[1]]%*%(beta+beta1) + e1
Y2 = X.test.list[[2]]%*%(beta+beta2) + e2
Y.test.list = list(Y1, Y2)

result = list()

ml.continuum = continuum.int.2step(X.list, Y.list, lambda = 100, gam = 0, rankJ = 1, rankA = c(1, 1))
Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

ml.pls = mycvplsr(X.list, Y.list)
Yhat.homo.pls.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pls$beta.C[[g]])
Yhat.heter.pls.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pls$beta.Cind[[g]])
MSE.intercept.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n))^2))
MSE.homo.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n) - Yhat.homo.pls.list[[g]])^2))
MSE.heter.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n) - Yhat.heter.pls.list[[g]])^2))
MSE.pls = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pls$intercept[[g]], n) - Yhat.homo.pls.list[[g]] - Yhat.heter.pls.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.pls, MSE.homo.pls, MSE.heter.pls, MSE.pls))

ml.continuum = continuum.int.2step(X.list, Y.list, lambda = 0, gam = 1, rankJ = ml.pls$rankJ, rankA = ml.pls$rankA)
Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

ml.continuum = continuum.int.2step(X.list, Y.list, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2))
Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

ml.pcr = mycvpcr(X.list, Y.list)
Yhat.homo.pcr.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pcr$beta.C[[g]])
Yhat.heter.pcr.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.pcr$beta.Cind[[g]])
MSE.intercept.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n))^2))
MSE.homo.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n) - Yhat.homo.pcr.list[[g]])^2))
MSE.heter.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n) - Yhat.heter.pcr.list[[g]])^2))
MSE.pcr = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.pcr$intercept[[g]], n) - Yhat.homo.pcr.list[[g]] - Yhat.heter.pcr.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.pcr, MSE.homo.pcr, MSE.heter.pcr, MSE.pcr))

ml.continuum = continuum.int.2step(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = ml.pcr$rankJ, rankA = ml.pcr$rankA)
Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

ml.continuum = continuum.int.2step(X.list, Y.list, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2))
Yhat.homo.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.C[[g]])
Yhat.heter.continuum.list = lapply(1:G, function(g) X.test.list[[g]]%*%ml.continuum$beta.Cind[[g]])
MSE.intercept.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n))^2))
MSE.homo.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]])^2))
MSE.heter.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.heter.continuum.list[[g]])^2))
MSE.continuum = sapply(1:G, function(g) mean((Y.test.list[[g]] - rep(ml.continuum$intercept[[g]], n) - Yhat.homo.continuum.list[[g]] - Yhat.heter.continuum.list[[g]])^2))
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

result = t(sapply(result, function(X) X[4,]))
row.names(result) = c("OLS", "PLS", "PLS.C", "PLS.Cg", "PCR", "PCR.C", "PCR.Cg")



