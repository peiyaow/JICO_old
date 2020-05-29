library(MASS)
library(pls)
library(purrr)
library(pracma)

r = 2
r1 = 2
r2 = 2
n1 = 50
n2 = 50
n = n1 + n2
p = 100
G = 2
L = c(10, 10)
L1 = c(2, 1)
L2 = c(2, 1)
alpha = c(100, 100)
alpha1 = rep(50, r1)
alpha2 = rep(50, r2)

Q = randortho(p, type = "orthonormal")

DIFF = list()
CONV = list()
NRUN = list()
RESULT = list()
NCOMP = list()


P = randortho(n, type = "orthonormal")
P1 = randortho(n1, type = "orthonormal")
P2 = randortho(n2, type = "orthonormal")
S = P[,1:r]%*%DIAG(L)
U = matrix(Q[1:r,], ncol = p)
J = S%*%U
J1 = J[1:n1,]
J2 = J[(n1+1):(n1+n2), ]
S1 = P1[,1:r1]%*%DIAG(L1)
S2 = P2[,1:r2]%*%DIAG(L2)
W1 = matrix(Q[(r+1):(r+r1),], ncol = p)
W2 = matrix(Q[(r+r1+1):(r+r1+r2),], ncol = p)
I1 = S1%*%W1
I2 = S2%*%W2
E1 = mvrnorm(n1, rep(0, p), diag(p)*.01)
E2 = mvrnorm(n2, rep(0, p), diag(p)*.01)

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e1 = rnorm(n1)
e2 = rnorm(n2)

Y1 = J1%*%t(U)%*%alpha + I1%*%t(W1)%*%alpha1 + e1
Y2 = J2%*%t(U)%*%alpha + I2%*%t(W2)%*%alpha2 + e2

cor(Y1, J1%*%t(U))
cor(Y1, I1%*%t(W1))

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)
cor(S, Y)
cor(S1, Y1)
cor(S2, Y2)

P = randortho(n, type = "orthonormal")
P1 = randortho(n1, type = "orthonormal")
P2 = randortho(n2, type = "orthonormal")
S = P[,1:r]%*%DIAG(L)
#U = matrix(Q[1:r,], ncol = p)
J = S%*%U
J1 = J[1:n1,]
J2 = J[(n1+1):(n1+n2), ]
S1 = P1[,1:r1]%*%DIAG(L1)
S2 = P2[,1:r2]%*%DIAG(L2)
#W1 = matrix(Q[(r+1):(r+r1),], ncol = p)
#W2 = matrix(Q[(r+r1+1):(r+r1+r2),], ncol = p)
I1 = S1%*%W1
I2 = S2%*%W2
E1 = mvrnorm(n1, rep(0, p), diag(p)*.01)
E2 = mvrnorm(n2, rep(0, p), diag(p)*.01)

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e1 = rnorm(n1)
e2 = rnorm(n2)

Y1 = J1%*%t(U)%*%alpha + I1%*%t(W1)%*%alpha1 + e1
Y2 = J2%*%t(U)%*%alpha + I2%*%t(W2)%*%alpha2 + e2

X.test = rbind(X1, X2)
Y.test = rbind(Y1, Y2)

X.test.list = list(X1, X2)
Y.test.list = list(Y1, Y2)








