m = 10
r = 2
r1 = 2
r2 = 2
n = 50
p1 = 50
p2 = 50
p = p1 + p2

X = mvrnorm(n, rep(0, p), diag(p))
# X.mean = apply(X, 2, mean)
# X = sweep(X, 2, X.mean)
Q = randortho(n, type = "orthonormal")
V = Q[,1:r]
J = V%*%solve(t(V)%*%V)%*%t(V)%*%X
X1 = X[,1:p1]
X2 = X[,(p1+1):p]
V1 = Q[,(1+r):(r+r1)]
V2 = Q[,(1+r1+r):(r+r1+r2)]
I1 = V1%*%solve(t(V1)%*%V1)%*%t(V1)%*%X1
I2 = V2%*%solve(t(V2)%*%V2)%*%t(V2)%*%X2
# t(I1)%*%I2
# t(I1)%*%J
# t(I2)%*%J

beta = rep(1/5, n)*10
beta1 = rep(c(1/100, -1/100), n/2)*10
beta2 = rep(c(0, 1/50), n/2)*10

Y = J%*%t(J)%*%beta + I1%*%t(I1)%*%beta1 + I2%*%t(I2)%*%beta2 #+ e
ml = continuum.ridge.fixm(J, Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = FALSE)
# alpha = C2beta((J)%*%t(J), Y, ml$a, lambda = 0)$alpha
alpha = C2beta((J)%*%t(J), Y, ml$C, lambda = 0)$alpha
alpha
#ml$a%*%alpha

# ml_ = continuum.ridge.fixm(J%*%t(J), Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = T)
# alpha_ = C2beta((J)%*%t(J), Y, ml_$C, lambda = 0)$alpha
# ml_$C%*%alpha_
 
ml1 = continuum.ridge.fixm(I1, Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = FALSE)
alpha1 = C2beta((I1)%*%t(I1), Y, ml1$C, lambda = 0)$alpha
alpha1

ml2 = continuum.ridge.fixm(I2, Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = FALSE)
alpha2 = C2beta((I2)%*%t(I2), Y, ml2$C, lambda = 0)$alpha
alpha2


Q = randortho(p, type = "orthonormal")
U = Q[,1:r]%*%(ml$E[1:r, 1:r]^(1/2))
# S = t(ml$V[,1:r])
C = ml$C
# J = U%*%S
J = U%*%t(C)
# a = ml$V[,1:r]%*%solve(ml$E[1:r, 1:r])^(1/2)%*%ml$Z

# Y = t(J)%*%(J)%*%a%*%alpha #+ e
Y = t(J)%*%(J)%*%C%*%alpha #+ e
ml = continuum.ridge.fix(t(J), Y, 1, lambda = 0, gam = 1, om = r, vertical = FALSE)
alpha.C = C2beta(t(J)%*%J, Y, ml$C, lambda = 0)$alpha
alpha.C

Q1 = randortho(p1, type = "orthonormal")
W1 = Q1[,1:r]%*%(ml1$E[1:r, 1:r]^(1/2))
S1 = t(ml1$V[,1:r])
I1 = W1%*%S1
a1 = ml1$V[,1:r]%*%solve(ml1$E[1:r, 1:r])^(1/2)%*%ml1$Z
# Y = t(I1)%*%(I1)%*%a1%*%alpha1 #+ e
# ml = continuum.ridge.fix(t(I1), Y, 1, lambda = 0, gam = 1, om = r1, vertical = FALSE)
# alpha.C = C2beta(t(I1)%*%I1, Y, ml$a, lambda = 0)$alpha
# alpha.C

Q2 = randortho(p2, type = "orthonormal")
W2 = Q2[,1:r]%*%(ml2$E[1:r, 1:r]^(1/2))
S2 = t(ml2$V[,1:r])
I2 = W2%*%S2
a2 = ml2$V[,1:r]%*%solve(ml2$E[1:r, 1:r])^(1/2)%*%ml2$Z[1:r,]
# Y = t(I2)%*%(I2)%*%a2%*%alpha2
# ml = continuum.ridge.fix(t(I2), Y, 1, lambda = 0, gam = 1, om = r2, vertical = FALSE)
# alpha.C = C2beta(t(I2)%*%I2, Y, ml$a, lambda = 0)$alpha
# alpha.C

E = t(mvrnorm(n, rep(0, p), diag(p))*.01)
X = J + rbind(I1, I2) + E
e = rnorm(n)*.01
Y = t(J)%*%(J)%*%a%*%alpha + t(I1)%*%(I1)%*%a1%*%alpha1 + t(I2)%*%(I2)%*%a2%*%alpha2 + e

X1 = X[1:p1,]
X2 = X[(p1+1):p,]
X.list = list(X1, X2)

# testing
Q = randortho(p, type = "orthonormal")
U = Q[,1:r]%*%(ml$E[1:r, 1:r]^(1/2))
# S = t(ml$V[,1:r])
J.test = U%*%S
# a = ml$V[,1:r]%*%solve(ml$E[1:r, 1:r])^(1/2)%*%ml$Z

Q1 = randortho(p1, type = "orthonormal")
W1 = Q1[,1:r]%*%(ml1$E[1:r, 1:r]^(1/2))
# S1 = t(ml1$V[,1:r])
I1.test = W1%*%S1
# a1 = ml1$V[,1:r]%*%solve(ml1$E[1:r, 1:r])^(1/2)%*%ml1$Z

Q2 = randortho(p2, type = "orthonormal")
W2 = Q2[,1:r]%*%(ml2$E[1:r, 1:r]^(1/2))
# S2 = t(ml2$V[,1:r])
I2.test = W2%*%S2
# a2 = ml2$V[,1:r]%*%solve(ml2$E[1:r, 1:r])^(1/2)%*%ml2$Z[1:r,]

E = t(mvrnorm(n, rep(0, p), diag(p))*.01)
X = J.test + rbind(I1.test, I2.test) + E
e = rnorm(n)*.01
Y.test = t(J.test)%*%(J)%*%a%*%alpha + t(I1.test)%*%(I1)%*%a1%*%alpha1 + t(I2.test)%*%(I2)%*%a2%*%alpha2 + e

X1 = X[1:p1,]
X2 = X[(p1+1):p,]
X.test.list = list(X1, X2)


# run
ml.continuum = continuum.multisource.iter(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2), 
                                          center.X = F, scale.X = F, center.Y = F, scale.Y = F)
mlll = decomposeX(X.test.list, ml.continuum$C, ml.continuum$Cind, ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
Yhat.homo = t(mlll$J)%*%ml.continuum$J%*%ml.continuum$beta.C
Yhat.heter.list = lapply(1:G, function(g) t(mlll$I[[g]])%*%ml.continuum$I[[g]]%*%ml.continuum$beta.Cind[[g]])
#Yhat.heter.list = lapply(1:G, function(g) t(mlll$I[[g]])%*%mlll$I[[g]]%*%ml.continuum$beta.Cind[[g]])
Yhat.heter = do.call("+", Yhat.heter.list)
MSE.intercept.continuum = mean((Y.test - ml.continuum$intercept)^2)
MSE.homo.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo)^2)
MSE.heter.continuum = mean((Y.test - ml.continuum$intercept - Yhat.heter)^2)
MSE.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo - Yhat.heter)^2)

ml.continuum = continuum.multisource.iter(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2), 
                                          center.X = F, scale.X = F, center.Y = F, scale.Y = F)

mlll = decomposeX(X.test.list, ml.continuum$C, ml.continuum$Cind, ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
Yhat.homo = t(mlll$J)%*%ml.continuum$J%*%ml.continuum$beta.C
Yhat.heter.list = lapply(1:G, function(g) t(mlll$I[[g]])%*%ml.continuum$I[[g]]%*%ml.continuum$beta.Cind[[g]])
#Yhat.heter.list = lapply(1:G, function(g) t(mlll$I[[g]])%*%mlll$I[[g]]%*%ml.continuum$beta.Cind[[g]])
Yhat.heter = do.call("+", Yhat.heter.list)
MSE.intercept.continuum = mean((Y.test - ml.continuum$intercept)^2)
MSE.homo.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo)^2)
MSE.heter.continuum = mean((Y.test - ml.continuum$intercept - Yhat.heter)^2)
MSE.continuum = mean((Y.test - ml.continuum$intercept - Yhat.homo - Yhat.heter)^2)





ml = continuum.ridge.fix(t(X), Y, 1, lambda = 0, gam = 1, om = r, vertical = FALSE)
Xa = X%*%ml$C%*%solve(t(ml$C)%*%ml$C)%*%t(ml$C)
# Xa = X%*%ml$V[,1:r]%*%solve(t(ml$V[,1:r])%*%ml$V[,1:r])%*%t(ml$V[,1:r])
C2beta(t(X)%*%(X), Y, ml$a, lambda = 0)$alpha
Y.homo = t(X)%*%(X)%*%C2beta(t(X)%*%(X), Y, ml$a, lambda = 0)$beta
Y.heter = matrix(Y - Y.homo)
# Xx = X%*%(diag(n) - ml$V[,1:r]%*%solve(t(ml$V[,1:r])%*%ml$V[,1:r])%*%t(ml$V[,1:r]))
Xx = X%*%(diag(n) - ml$C%*%solve(t(ml$C)%*%ml$C)%*%t(ml$C))
Xx1 = Xx[1:p1,]
Xx2 = Xx[(1+p1):p,]
ml1 = continuum.ridge.fix(t(Xx1), Y.heter, 1, lambda = 0, gam = 1, om = r, vertical = FALSE)
ml2 = continuum.ridge.fix(t(Xx2), Y.heter, 1, lambda = 0, gam = 1, om = r, vertical = FALSE)
C2beta(t(Xx1)%*%(Xx1), Y.heter, ml1$a, lambda = 0)$alpha
C2beta(t(Xx2)%*%(Xx2), Y.heter, ml2$a, lambda = 0)$alpha





