library(MASS)
library(pls)
library(glmnet)

source("~/continuum/simulation/function/function.R")

p = 100
n = 100
H1 = cbind(matrix(3, ncol = 50, nrow = 50), matrix(4, ncol = 50, nrow = 50))
H2 = matrix(3.5, ncol = n, nrow = p-50)
H = rbind(H1, H2)
dim(H2) # p by n
beta = c(rep(1/25, 50), rep(1/25, p-50))
beta1 = c(rep(1/50, 25), rep(0, p-25))
beta2 = c(rep(-1/50, 25), rep(0, p-25))

for (i in 1:10){
    # -------------------------------------- generate data ---------------------------------------
    E = mvrnorm(p, rep(0, n), diag(n))
    X = H + E
    X = t(X)
    f1 = rnorm(n/2)*1.5
    f2 = rnorm(n/2)*1.5
    
    X1 = X[1:50,]
    X2 = X[51:100,]
    X.list = list(X1, X2)
    Y1 = X1%*%(beta+beta1) + f1
    Y2 = X2%*%(beta+beta2) + f2
    Y.list = list(Y1, Y2)
    Y = c(Y1, Y2)
    
    E.test = mvrnorm(p, rep(0, n), diag(n))
    X.test = H + E.test
    X.test = t(X.test)
    f1.test = rnorm(n/2)*1.5
    f2.test = rnorm(n/2)*1.5
    
    X1.test = X.test[1:50,]
    X2.test = X.test[51:100,]
    X.test.list = list(X1.test, X2.test)
    Y1.test = X1.test%*%(beta+beta1) + f1.test
    Y2.test = X2.test%*%(beta+beta2) + f2.test
    Y.test.list = list(Y1.test, Y2.test)
    Y.test = c(Y1.test, Y2.test)
    # ------------------------------------------------------------------------------------------
    
    # -------------------------------- global method -------------------------------------------
    ml.pls = plsr(Y ~ X, validation = "CV", scale = T)
    # ncomp.pls = selectNcomp(ml.pls, method = "onesigma", plot = F)
    ncomp.pls = which.min(RMSEP(ml.pls)$val[1,,])-1
    
    ml.pcr = pcr(Y~X, validation = "CV", scale = T)
    # ncomp.pcr = selectNcomp(ml.pcr, method = "onesigma", plot = F)
    ncomp.pcr = which.min(RMSEP(ml.pcr)$val[1,,])-1
    
    ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0)
    ml.lasso = cv.glmnet(x = X, y = Y, alpha = 1)
    
    ml.continuum = cv.continuum.ridge(X, Y, lambda = 0, gam = .5)
    C.pls = ml.continuum$C
    beta.C.pls = C2beta(X, Y, ml.continuum$C, ml.continuum$lam)$coef
    
    ml.continuum = cv.continuum.ridge(X, Y, lambda = 0, gam = 1e10)
    C.pcr = ml.continuum$C
    beta.C.pcr = C2beta(X, Y, ml.continuum$C, ml.continuum$lam)$coef
    
    
    ml.ridge.scale = cv.glmnet(x = X, y = scale(Y), alpha = 0)
    lambda = ml.ridge.scale$lambda.min
    ml.continuum = cv.continuum.ridge(X, Y, lambda = lambda, gam = 0)
    C.ridge0 = ml.continuum$C[,1]
    beta.C.ridge0 = C2beta(X, Y, C.ridge0, ml.continuum$lam)$coef
    C.ridge = ml.continuum$C
    beta.C.ridge = C2beta(X, Y, C.ridge, ml.continuum$lam)$coef
    
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
    Yhat.lasso = predict(ml.lasso, newx = X.test, s = ml.lasso$lambda.min)
    Yhat.continuum.ridge0 = cbind(1, X.test)%*%beta.C.ridge0
    Yhat.continuum.ridge = cbind(1, X.test)%*%beta.C.ridge
    Yhat.continuum.pls = cbind(1, X.test)%*%beta.C.pls
    Yhat.continuum.pcr = cbind(1, X.test)%*%beta.C.pcr
    
    mse.pls = mean((Yhat.pls - Y.test)^2)
    mse.pcr = mean((Yhat.pcr - Y.test)^2)
    mse.ridge = mean((Yhat.ridge - Y.test)^2)
    mse.lasso = mean((Yhat.lasso - Y.test)^2)
    mse.continuum.pls = mean((Yhat.continuum.pls - Y.test)^2)
    mse.continuum.pcr = mean((Yhat.continuum.pcr - Y.test)^2)
    mse.continuum.ridge = mean((Yhat.continuum.ridge - Y.test)^2)
    mse.continuum.ridge0 = mean((Yhat.continuum.ridge0 - Y.test)^2)
    
    c(mse.pls, mse.pcr, mse.ridge, mse.lasso, mse.continuum.pls, mse.continuum.pcr, mse.continuum.ridge, mse.continuum.ridge0)
    # ------------------------------------------------------------------------------------------
    
    # ---------------------------------- groupwise method --------------------------------------
    ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", scale = T))
    ncomp.pls.list = lapply(1:G, function(g) which.min(RMSEP(ml.pls.list[[g]])$val[1,,])-1)
    # ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "onesigma", plot = F))
    
    ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", scale = T))
    ncomp.pcr.list = lapply(1:G, function(g) which.min(RMSEP(ml.pcr.list[[g]])$val[1,,])-1)
    # ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "onesigma", plot = F))
    
    ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0))
    ml.lasso.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 1))
    
    ml.continuum.pls.list = lapply(1:G, function(g) cv.continuum.ridge(X.list[[g]], Y.list[[g]], 0, gam = 1/2))
    beta.C.pls.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.list[[g]], ml.continuum.pls.list[[g]]$C, ml.continuum.pls.list[[g]]$lam)$coef)
    
    ml.continuum.pcr.list = lapply(1:G, function(g) cv.continuum.ridge(X.list[[g]], Y.list[[g]], 0, gam = 1e10))
    beta.C.pcr.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.list[[g]], ml.continuum.pcr.list[[g]]$C, ml.continuum.pcr.list[[g]]$lam)$coef)
    
    ml.ridge.scale.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = scale(Y.list[[g]]), alpha = 0))
    lambda.list = lapply(1:G, function(g) ml.ridge.scale.list[[g]]$lambda.min)
    ml.continuum.ridge.list = lapply(1:G, function(g) cv.continuum.ridge(X.list[[g]], Y.list[[g]], lambda.list[[g]], gam = 0))
    beta.C.ridge.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.list[[g]], ml.continuum.ridge.list[[g]]$C, ml.continuum.ridge.list[[g]]$lam)$coef)
    beta.C.ridge0.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.list[[g]], ml.continuum.ridge.list[[g]]$C[,1], ml.continuum.ridge.list[[g]]$lam)$coef)
    
    Yhat.pls.list = lapply(1:G, function(g) if (ncomp.pls.list[[g]]){
      predict(ml.pls.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
    } else{
      rep(ml.pls.list[[g]]$Ymeans, n/2)
    })
    Yhat.pcr.list = lapply(1:G, function(g) if (ncomp.pcr.list[[g]]){
      predict(ml.pcr.list[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pcr.list[[g]])[,,1]
    } else{
      rep(ml.pcr.list[[g]]$Ymeans, n/2)
    })
    Yhat.ridge.list = lapply(1:G, function(g) predict(ml.ridge.list[[g]], newx = X.test.list[[g]], s = ml.ridge.list[[g]]$lambda.min))
    Yhat.lasso.list = lapply(1:G, function(g) predict(ml.lasso.list[[g]], newx = X.test.list[[g]], s = ml.lasso.list[[g]]$lambda.min))
    Yhat.continuum.ridge.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.ridge.list[[g]])
    Yhat.continuum.ridge0.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.ridge0.list[[g]])
    Yhat.continuum.pls.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.pls.list[[g]])
    Yhat.continuum.pcr.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.pcr.list[[g]])
    
    mse.pls.group.vec = sapply(1:G, function(g) mean((Yhat.pls.list[[g]] - Y.test.list[[g]])^2))
    mse.pcr.group.vec = sapply(1:G, function(g) mean((Yhat.pcr.list[[g]] - Y.test.list[[g]])^2))
    mse.ridge.group.vec = sapply(1:G, function(g) mean((Yhat.ridge.list[[g]] - Y.test.list[[g]])^2))
    mse.lasso.group.vec = sapply(1:G, function(g) mean((Yhat.lasso.list[[g]] - Y.test.list[[g]])^2))
    mse.continuum.ridge.group.vec = sapply(1:G, function(g) mean((Yhat.continuum.ridge.list[[g]] - Y.test.list[[g]])^2))
    mse.continuum.ridge0.group.vec = sapply(1:G, function(g) mean((Yhat.continuum.ridge0.list[[g]] - Y.test.list[[g]])^2))
    mse.continuum.pls.group.vec = sapply(1:G, function(g) mean((Yhat.continuum.pls.list[[g]] - Y.test.list[[g]])^2))
    mse.continuum.pcr.group.vec = sapply(1:G, function(g) mean((Yhat.continuum.pcr.list[[g]] - Y.test.list[[g]])^2))
    
    mse.pls.group = mean(mse.pls.group.vec)
    mse.pcr.group = mean(mse.pcr.group.vec)
    mse.ridge.group = mean(mse.ridge.group.vec)
    mse.lasso.group = mean(mse.lasso.group.vec)
    mse.continuum.ridge.group = mean(mse.continuum.ridge.group.vec)
    mse.continuum.ridge0.group = mean(mse.continuum.ridge0.group.vec)
    mse.continuum.pls.group = mean(mse.continuum.pls.group.vec)
    mse.continuum.pcr.group = mean(mse.continuum.pcr.group.vec)
    
    c(mse.pls.group, mse.pcr.group, mse.ridge.group, mse.lasso.group, mse.continuum.pls.group, mse.continuum.pcr.group, mse.continuum.ridge.group, mse.continuum.ridge0.group)
    # --------------------------------------------------------------------------------------------
    
    # ----------------------------------- integrated method --------------------------------------
    # C = ml.continuum$C
    # t(C)%*%t(scale(X))%*%scale(X)%*%C
    # t(ml.pcr$projection)%*%t(scale(X))%*%scale(X)%*%ml.pcr$projection
    
    # continuum ridge0
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.list[[g]])%*%beta.C.ridge0)
    Y.res.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.continuum.list[[g]])
    C.res.list = lapply(1:G, function(g) continuum.ridge.res(X.list[[g]], Y.res.list[[g]], C.ridge0, lambda, gam = 0))
    beta.C.res.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.res.list[[g]], C.res.list[[g]][,1], 0)$coef)
    Yhat.continuum.res.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.res.list[[g]])
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.ridge0)
    mse.continuum.ridge0.int.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]] + Yhat.continuum.res.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.ridge0.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.ridge0.int = mean(mse.continuum.ridge0.int.vec)
    
    g = 1
    t(C.res.list[[g]])%*%t(scale(X.list[[g]]))%*%scale(X.list[[g]])%*%C.ridge0
    
    
    # continuum ridge
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.list[[g]])%*%beta.C.ridge)
    Y.res.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.continuum.list[[g]])
    C.res.list = lapply(1:G, function(g) continuum.ridge.res(X.list[[g]], Y.res.list[[g]], C.ridge, lambda, gam = 0))
    beta.C.res.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.res.list[[g]], C.res.list[[g]][,1], 0)$coef)
    Yhat.continuum.res.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.res.list[[g]])
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.ridge)
    mse.continuum.ridge.int.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]] + Yhat.continuum.res.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.ridge.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.ridge.int = mean(mse.continuum.ridge.int.vec)
    
    g = 2
    t(C.res.list[[g]])%*%t(scale(X.list[[g]]))%*%scale(X.list[[g]])%*%C.ridge
    
    
    # continuum pls
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.list[[g]])%*%beta.C.pls)
    Y.res.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.continuum.list[[g]])
    C.res.list = lapply(1:G, function(g) continuum.ridge.res(X.list[[g]], Y.res.list[[g]], C.pls, 0, gam = 1/2))
    beta.C.res.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.res.list[[g]], C.res.list[[g]][,1], 0)$coef)
    Yhat.continuum.res.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.res.list[[g]])
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.pls)
    mse.continuum.pls.int.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]] + Yhat.continuum.res.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.pls.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.pls.int = mean(mse.continuum.pls.int.vec)
    
    g = 1
    t(C.res.list[[g]])%*%t(scale(X.list[[g]]))%*%scale(X.list[[g]])%*%C.pls
    
    # continuum pcr
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.list[[g]])%*%beta.C.pcr)
    Y.res.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.continuum.list[[g]])
    C.res.list = lapply(1:G, function(g) continuum.ridge.res(X.list[[g]], Y.res.list[[g]], C.pcr, 0, gam = 1e10))
    beta.C.res.list = lapply(1:G, function(g) C2beta(X.list[[g]], Y.res.list[[g]], C.res.list[[g]][,1], 0)$coef)
    Yhat.continuum.res.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.res.list[[g]])
    Yhat.continuum.list = lapply(1:G, function(g) cbind(1, X.test.list[[g]])%*%beta.C.pcr)
    mse.continuum.pcr.int.vec = sapply(1:G, function(g) mean((Yhat.continuum.list[[g]] + Yhat.continuum.res.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.pcr.vec = lapply(1:G, function(g) mean((Yhat.continuum.list[[g]]-Y.test.list[[g]])^2))
    mse.continuum.pcr.int = mean(mse.continuum.pcr.int.vec)
    
    g = 1
    t(C.res.list[[g]])%*%t(scale(X.list[[g]]))%*%scale(X.list[[g]])%*%C.pcr
    
    write.table(t(c(c(mse.pls, mse.pcr, mse.ridge, mse.lasso, mse.continuum.pls, mse.continuum.pcr, mse.continuum.ridge, mse.continuum.ridge0),
                    c(mse.pls.group, mse.pcr.group, mse.ridge.group, mse.lasso.group, mse.continuum.pls.group, mse.continuum.pcr.group, mse.continuum.ridge.group, mse.continuum.ridge0.group),
                    c(mse.continuum.pls.int, mse.continuum.pcr.int, mse.continuum.ridge.int, mse.continuum.ridge0.int))), file = "MSE.csv", sep = ',', append = T, col.names = F, row.names = F)
    
    write.table(t(c(c(mse.continuum.pls.vec, mse.continuum.pcr.vec, mse.continuum.ridge.vec, mse.continuum.ridge0.vec),
                    c(mse.pls.group.vec, mse.pcr.group.vec, mse.ridge.group.vec, mse.lasso.group.vec, mse.continuum.pls.group.vec, mse.continuum.pcr.group.vec, mse.continuum.ridge.group.vec, mse.continuum.ridge0.group.vec),
                    c(mse.continuum.pls.int.vec, mse.continuum.pcr.int.vec, mse.continuum.ridge.int.vec, mse.continuum.ridge0.int.vec))), file = "MSE_group.csv", sep = ',', append = T, col.names = F, row.names = F)
    
}
