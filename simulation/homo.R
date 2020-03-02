library(pls)
library(glmnet)
library(MASS)

n = 100
p = 50
s = 0.5
Sigma = diag(p)*s^2
e = 0.5
beta = rep(1, p)

X1 = mvrnorm(n, rep(0, p), Sigma)
X2 = mvrnorm(n, rep(0, p), Sigma)
eps = rnorm(n, mean = 0, sd = e)
Y = X1%*%beta + eps 
X = cbind(X1, X2)

X1.test = mvrnorm(n, rep(0, p), Sigma)
X2.test = mvrnorm(n, rep(0, p), Sigma)
eps.test = rnorm(n, mean = 0, sd = e)
Y.test = X1.test%*%beta + eps.test 
X.test = cbind(X1.test, X2.test)

ml.ols = lm(Y~X)
ml.pls = plsr(Y ~ X, validation = "CV")
ncomp.pls = selectNcomp(ml.pls, method = "onesigma", plot = T)
ml.pcr = pcr(Y~X, validation = "CV")
ncomp.pcr = selectNcomp(ml.pcr, method = "onesigma", plot = T)
ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0)

Yhat.ols = cbind(1, X.test)%*%ml.ols$coefficients
Yhat.pls = predict(ml.pls, newdata = X.test, ncomp = ncomp.pls)[,,1]
Yhat.pcr = predict(ml.pcr, newdata = X.test, ncomp = ncomp.pcr)[,,1]
Yhat.ridge = predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.1se)

mse.ols = mean((Yhat.ols - Y.test)^2)
mse.pls = mean((Yhat.pls - Y.test)^2)
mse.pcr = mean((Yhat.pcr - Y.test)^2)
mse.ridge = mean((Yhat.ridge - Y.test)^2)
