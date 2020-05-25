m = 10
r = 2
r1 = 2
r2 = 2
n = 50
p1 = 50
p2 = 50
p = p1 + p2

X = t(mvrnorm(n, rep(0, p), diag(p)))
# X.mean = apply(X, 2, mean)
# X = sweep(X, 2, X.mean)
Q = randortho(n, type = "orthonormal")
V = Q[,1:r]
# J = V%*%solve(t(V)%*%V)%*%t(V)%*%X
X1 = X[1:p1,]
X2 = X[(p1+1):p,]
V1 = Q[,(1+r):(r+r1)]
V2 = Q[,(1+r1+r):(r+r1+r2)]
#I1 = V1%*%solve(t(V1)%*%V1)%*%t(V1)%*%X1
#I2 = V2%*%solve(t(V2)%*%V2)%*%t(V2)%*%X2

J = t(X)%*%X%*%V%*%solve(t(V)%*%V)%*%t(V)
I1 = t(X1)%*%X1%*%V1%*%solve(t(V1)%*%V1)%*%t(V1)
I2 = t(X2)%*%X2%*%V2%*%solve(t(V2)%*%V2)%*%t(V2)

beta = rep(1/5, n)*10
beta1 = rep(c(1/100, -1/100), n/2)*10
beta2 = rep(c(0, 1/50), n/2)*10

Y = J%*%beta + I1%*%beta1 + I2%*%beta2 #+ e


ml = continuum.ridge.fixm(J, Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = T)
alpha = C2beta(J, Y, ml$C, lambda = 0)$alpha
alpha

ml1 = continuum.ridge.fixm(I1, Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = T)
alpha1 = C2beta(I1, Y, ml1$C, lambda = 0)$alpha
alpha1

ml2 = continuum.ridge.fixm(I2, Y, 1, lambda = 0, gam = 1, om = 2, m = r, vertical = T)
alpha2 = C2beta(I2, Y, ml2$C, lambda = 0)$alpha
alpha2

J = ml$V%*%(ml$E)^(1/2)%*%t(ml$V)
Y = J%*%ml$C%*%alpha
ml_ = continuum.ridge.fix(J, Y, 1, lambda = 0, gam = 1, om = r, vertical = T)
alpha.C = C2beta(J, Y, ml_$C, lambda = 0)$alpha
alpha.C

# 
# I1 = ml1$V%*%(ml1$E)^(1/2)%*%t(ml1$V)
# Y = I1%*%ml1$C%*%alpha1
# ml1_ = continuum.ridge.fix(I1, Y, 1, lambda = 0, gam = 1, om = r1, vertical = T)
# alpha.C = C2beta(I1, Y, ml1_$C, lambda = 0)$alpha
# alpha.C
# 
# I2 = ml2$V%*%(ml2$E)^(1/2)%*%t(ml2$V)
# Y = I2%*%ml2$C%*%alpha2
# ml2_ = continuum.ridge.fix(I2, Y, 1, lambda = 0, gam = 1, om = r2, vertical = T)
# alpha.C = C2beta(I2, Y, ml2_$C, lambda = 0)$alpha
# alpha.C

K = J%*%C%*%solve(t(C)%*%C)%*%(t(C)) + I1%*%C1%*%solve(t(C1)%*%C1)%*%(t(C1)) + I2%*%C2%*%solve(t(C2)%*%C2)%*%(t(C2))

Q = randortho(p, type = "orthonormal")
U = Q[,1:r]%*%(ml$E^(1/4))
S = t(ml$V)
J = U%*%S
#t(J)%*%J - ml$V%*%(ml$E)^(1/2)%*%t(ml$V)

Q1 = randortho(p1, type = "orthonormal")
W1 = Q1[,1:r1]%*%(ml1$E^(1/4))
S1 = t(ml1$V)
I1 = W1%*%S1
#t(I1)%*%I1 - ml1$V%*%(ml1$E)^(1/2)%*%t(ml1$V)

Q2 = randortho(p2, type = "orthonormal")
W2 = Q2[,1:r2]%*%(ml2$E^(1/4))
S2 = t(ml2$V)
I2 = W2%*%S2
#t(I2)%*%I2 - ml2$V%*%(ml2$E)^(1/2)%*%t(ml2$V)

C = ml$C
C1 = ml1$C
C2 = ml2$C

E = t(mvrnorm(n, rep(0, p), diag(p))*.01)
X = J + rbind(I1, I2) #+ E
e = rnorm(n)*.01
Y = t(J)%*%(J)%*%C%*%alpha + t(I1)%*%(I1)%*%C1%*%alpha1 + t(I2)%*%(I2)%*%C2%*%alpha2 #+ e

K = t(X)%*%X
ml = continuum.ridge.fix(K, Y, 1, lambda = 0, gam = 1, om = r, vertical = F)
alpha_ = C2beta(K, Y, ml$C, lambda = 0)$alpha
alpha_
alpha

X1 = X[1:p1,]
X2 = X[(p1+1):p,]

Yhat = K%*%ml$C%*%alpha_
Y.heter = Y - Yhat

K1 = t(X1)%*%X1%*%(diag(n) - ml$C%*%solve(t(ml$C)%*%ml$C)%*%t(ml$C))
ml1 = continuum.ridge.fix(K1, Y.heter, 1, lambda = 0, gam = 1, om = r1, vertical = T)
C2beta(K1, Y, ml1$C, lambda = 0)$alpha
alpha1

K2 = t(X2)%*%X2%*%(diag(n) - ml$C%*%solve(t(ml$C)%*%ml$C)%*%t(ml$C))
ml2 = continuum.ridge.fix(K2, Y.heter, 1, lambda = 0, gam = 1, om = r2, vertical = T)
C2beta(K2, Y, ml2$C, lambda = 0)$alpha
alpha2




