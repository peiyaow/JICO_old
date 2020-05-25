library(MASS)
library(pls)
library(purrr)
library(pracma)

r = 2
r1 = 2
r2 = 2
n = 50
p1 = 50
p2 = 50
L = c(50, 10)
L1 = c(1, 1)
L2 = c(1, 1)
p = p1 + p2
G = 2
alpha = rep(1, r)
alpha1 = rep(1, r1)
alpha2 = rep(1, r2)

Q = randortho(p, type = "orthonormal")
U = Q[,1:r]%*%diag(L)
Q1 = randortho(p1, type = "orthonormal")
Q2 = randortho(p2, type = "orthonormal")
W1 = Q1[,1:r1]%*%DIAG(L1)
W2 = Q2[,1:r2]%*%DIAG(L2)

P = randortho(n, type = "orthonormal")
S = matrix(P[1:r,], ncol = n)
S1 = matrix(P[(r+1):(r+r1),], ncol = n)
S2 = matrix(P[(r+r1+1):(r+r1+r2),], ncol = n)

E1 = t(mvrnorm(n, rep(0, p1), diag(p1))*0.01)
E2 = t(mvrnorm(n, rep(0, p2), diag(p2))*0.01)

J = U%*%S
J1 = J[1:p1,]
J2 = J[(p1+1):(p1+p2),]
I1 = W1%*%S1
I2 = W2%*%S2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e = rnorm(n)*0.01
Y = t(S)%*%alpha + t(S1)%*%alpha1 + t(S2)%*%alpha2 + e

X = rbind(X1, X2)
X.list = list(X1, X2)

# ----------------------------------------------- testing ----------------------------------------------------
P = randortho(n, type = "orthonormal")
S = matrix(P[1:r,], ncol = n)
S1 = matrix(P[(r+1):(r+r1),], ncol = n)
S2 = matrix(P[(r+r1+1):(r+r1+r2),], ncol = n)

E1 = t(mvrnorm(n, rep(0, p1), diag(p1))*0.01)
E2 = t(mvrnorm(n, rep(0, p2), diag(p2))*0.01)

J = U%*%S
J1 = J[1:p1,]
J2 = J[(p1+1):(p1+p2),]
I1 = W1%*%S1
I2 = W2%*%S2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e = rnorm(n)*0.01
Y.test = t(S)%*%alpha + t(S1)%*%alpha1 + t(S2)%*%alpha2 + e

X.test = rbind(X1, X2)
X.test.list = list(X1, X2)

# ----------------------------------------------- run model ----------------------------------------------------
result = list()
ml.jive = jive.multisource(X.list, rankJ = r, rankA = c(r1, r2), method = "given", orthIndiv = T)

ml.continuum.pcr = continuum.multisource.iter(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2))
ml.continuum = ml.continuum.pcr
ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
Yhat.homo = t(ml$J)%*%ml.continuum$J%*%ml.continuum$beta.C
Yhat.heter.list = lapply(1:G, function(g) t(ml$I[[g]])%*%ml.continuum$I[[g]]%*%ml.continuum$beta.Cind[[g]])
Yhat.heter = do.call("+", Yhat.heter.list)
MSE.intercept.continuum = mean((Y.test - ml.continuum$intercept)^2)
MSE.homo.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo)^2)
MSE.heter.continuum = mean((Y.test - ml.continuum$intercept - Yhat.heter)^2)
MSE.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo - Yhat.heter)^2)
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

diff = mean(ml.continuum$J - do.call(rbind, ml.jive$joint))^2 + mean((do.call(rbind, ml.continuum$I) - do.call(rbind, ml.jive$individual))^2)
diff

ml.continuum.pls = continuum.multisource.iter(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2))
ml.continuum = ml.continuum.pls
ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
Yhat.homo = t(ml$J)%*%ml.continuum$J%*%ml.continuum$beta.C
Yhat.heter.list = lapply(1:G, function(g) t(ml$I[[g]])%*%ml.continuum$I[[g]]%*%ml.continuum$beta.Cind[[g]])
Yhat.heter = do.call("+", Yhat.heter.list)
MSE.intercept.continuum = mean((Y.test - ml.continuum$intercept)^2)
MSE.homo.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo)^2)
MSE.heter.continuum = mean((Y.test - ml.continuum$intercept - Yhat.heter)^2)
MSE.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo - Yhat.heter)^2)
result = list.append(result, rbind(MSE.intercept.continuum, MSE.homo.continuum, MSE.heter.continuum, MSE.continuum))

result




