# m = 5
r = 2
r1 = 2
r2 = 2
n = 40
p1 = 40
p2 = 40
p = p1 + p2

X = t(mvrnorm(n, rep(0, p), diag(p)))
X1 = X[1:p1,]
X2 = X[(p1+1):p,]

beta = rep(1/5, n)*2
beta1 = rep(c(1/100, -1/100), n/2)*1
beta2 = rep(c(0, -1/100, 0, 1/100), n/4)*1

# t(beta)%*%beta1
# t(beta)%*%beta2
# beta1%*%beta2
# Y = t(X)%*%beta #+ t(X1)%*%X1%*%beta1 + t(X2)%*%X2%*%beta2
Y = t(X)%*%X%*%beta #+ t(X1)%*%X1%*%beta1 + t(X2)%*%X2%*%beta2
# ml = continuum.ridge.fixm(t(X)%*%X, Y, 1, lambda = 0, gam = 1, om = r, m = m, vertical = T)
ml = continuum.ridge.fix(t(X)%*%X, Y, 1, lambda = 0, gam = 1, om = r, vertical = T)
alpha = C2beta(t(X)%*%X, Y, ml$C, lambda = 0)$alpha
alpha
C = ml$C
# C%*%alpha



# Y = t(X1)%*%beta1
# t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))
Y = t(X1)%*%X1%*%beta1
ml1 = continuum.ridge.fix(t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C)), Y, 1, lambda = 0, gam = 1, om = r1, vertical = T)
C1 = ml1$C
alpha1 = C2beta(t(X1)%*%X1, Y, ml1$C, lambda = 0)$alpha
alpha1
# C1%*%alpha1

# t(C1)%*%C

Y = t(X2)%*%X2%*%beta2
ml2 = continuum.ridge.fix(t(X2)%*%X2%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C)), Y, 1, lambda = 0, gam = 1, om = r2, vertical = T)
C2 = ml2$C
alpha2 = C2beta(t(X2)%*%X2, Y, ml2$C, lambda = 0)$alpha
alpha2


# orthogonal
# Y = t(X1)%*%beta1
# t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))
Y = t(X1)%*%X1%*%beta1
ml1 = continuum.ridge.fix(t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))%*%(diag(n) - C2%*%solve(t(C2)%*%C2)%*%t(C2)), 
                           Y, 1, lambda = 0, gam = 1, om = r1, vertical = T)
C1 = ml1$C


# Y = t(X2)%*%beta2
# t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))
Y = t(X2)%*%X2%*%beta2
ml2 = continuum.ridge.fix(t(X2)%*%X2%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))%*%(diag(n) - C1%*%solve(t(C1)%*%C1)%*%t(C1)), 
                           Y, 1, lambda = 0, gam = 1, om = r2, vertical = T)
C2 = ml2$C

# Y = t(X1)%*%beta1
Y = t(X1)%*%X1%*%beta1
ml1 = continuum.ridge.fix(t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))%*%(diag(n) - C2%*%solve(t(C2)%*%C2)%*%t(C2)), 
                           Y, 1, lambda = 0, gam = 1, om = r1, vertical = T)
C1 = ml1$C
alpha1 = C2beta(t(X1)%*%X1, Y, C1, lambda = 0)$alpha
#alpha1

# Y = t(X2)%*%beta2
Y = t(X2)%*%X2%*%beta2
ml2 = continuum.ridge.fix(t(X2)%*%X2%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C))%*%(diag(n) - C1%*%solve(t(C1)%*%C1)%*%t(C1)), 
                           Y, 1, lambda = 0, gam = 1, om = r2, vertical = T)
C2 = ml2$C
alpha2 = C2beta(t(X2)%*%X2, Y, C2, lambda = 0)$alpha
#alpha2

alpha
alpha1
alpha2

t(C)%*%C1
t(C)%*%C2
t(C1)%*%C2

Y = t(X)%*%X%*%C%*%alpha + t(X1)%*%X1%*%C1%*%alpha1 + t(X2)%*%X2%*%C2%*%alpha2
ml = continuum.ridge.fix(t(X)%*%X, Y, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_ = C2beta(t(X)%*%X, Y, ml$C, lambda = 0)$alpha
alpha_

C = ml$C
Yhat = t(X)%*%X%*%C%*%alpha_
Y.heter = Y - Yhat

ml1 = continuum.ridge.fix(t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C)), 
                         Y.heter, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_1 = C2beta(t(X1)%*%X1, Y.heter, ml1$C, lambda = 0)$alpha
# ml$E
# ml1$E
# ml2$E

ml2 = continuum.ridge.fix(t(X2)%*%X2%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C)), 
                          Y.heter, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_2 = C2beta(t(X2)%*%X2, Y.heter, ml2$C, lambda = 0)$alpha

Yhat1 = t(X1)%*%X1%*%ml1$C%*%alpha_1
Yhat2 = t(X2)%*%X2%*%ml2$C%*%alpha_2

Y.homo = Y - Yhat1 - Yhat2
ml = continuum.ridge.fix(t(X)%*%X, Y.homo, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_ = C2beta(t(X)%*%X, Y.homo, ml$C, lambda = 0)$alpha
alpha_
alpha

C = ml$C
Yhat = t(X)%*%X%*%C%*%alpha_
Y.heter = Y - Yhat

ml1 = continuum.ridge.fix(t(X1)%*%X1%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C)), 
                          Y.heter, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_1 = C2beta(t(X1)%*%X1, Y.heter, ml1$C, lambda = 0)$alpha
alpha_1
alpha1

ml2 = continuum.ridge.fix(t(X2)%*%X2%*%(diag(n) - C%*%solve(t(C)%*%C)%*%t(C)), 
                          Y.heter, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_2 = C2beta(t(X2)%*%X2, Y.heter, ml2$C, lambda = 0)$alpha
alpha_2
alpha2




Q = randortho(p, type = "orthonormal")
U = Q[,1:m]%*%(ml$E^(1/4))
S = t(ml$V)
J = U%*%S

t(J)%*%J%*%C

Q1 = randortho(p1, type = "orthonormal")
W1 = Q1[,1:m]%*%(ml1$E^(1/4))
S1 = t(ml1$V)
I1 = W1%*%S1
#t(I1)%*%I1 - ml1$V%*%(ml1$E)^(1/2)%*%t(ml1$V)

Q2 = randortho(p2, type = "orthonormal")
W2 = Q2[,1:m]%*%(ml2$E^(1/4))
S2 = t(ml2$V)
I2 = W2%*%S2


E = t(mvrnorm(n, rep(0, p), diag(p))*.01)
X = J + rbind(I1, I2) #+ E
X1 = X[1:p1,]
X2 = X[(p1+1):p,]
e = rnorm(n)*.01
Y = t(X)%*%(J)%*%C%*%alpha + t(X1)%*%(I1)%*%C1%*%alpha1 + t(X2)%*%(I2)%*%C2%*%alpha2 #+ e

I1%*%C2
I2%*%C1

J%*%C1
J%*%C2

K = t(X)%*%X
ml = continuum.ridge.fix(K, Y, 1, lambda = 0, gam = 1, om = 2, vertical = T)
alpha_ = C2beta(K, Y, ml$C, lambda = 0)$alpha
alpha_
alpha

Yhat = K%*%ml$C%*%alpha_
Y.heter = Y - Yhat

K1 = t(X1)%*%X1%*%(diag(n) - ml$C%*%solve(t(ml$C)%*%ml$C)%*%t(ml$C))
ml1 = continuum.ridge.fix(K1, Y.heter, 1, lambda = 0, gam = 1, om = r1, vertical = T)
C2beta(K1, Y.heter, ml1$C, lambda = 0)$alpha
alpha1

K2 = t(X2)%*%X2%*%(diag(n) - ml$C%*%solve(t(ml$C)%*%ml$C)%*%t(ml$C))
ml2 = continuum.ridge.fix(K2, Y.heter, 1, lambda = 0, gam = 1, om = r2, vertical = T)
C2beta(K2, Y.heter, ml2$C, lambda = 0)$alpha
alpha2

