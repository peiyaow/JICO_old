# PCR
m = 5
r = 2
L = c(c(10, 3), rep(.1, m-r))

Q = randortho(p, type = "orthonormal")
U = Q[,1:m]%*%DIAG(L)
D = t(U)%*%U
P = randortho(n, type = "orthonormal")
S = matrix(P[1:m,], ncol = n)
E = t(mvrnorm(n, rep(0, p), diag(p))*0.01)

X = U%*%S + E
#Z = randortho(m, type = "orthonormal")[,1:r]
Z = diag(m)[,1:r]
alpha = c(1, 1)
e = rnorm(n)*0.01
a = t(S)%*%DIAG(1/L)%*%Z
Y = t(S)%*%DIAG(L)%*%Z%*%alpha #+ e
t(a)%*%t(X)%*%X%*%t(X)%*%X%*%a
t(a)%*%t(X)%*%X%*%Y

ml = continuum.ridge.fix(t(X), Y, 1, lambda = 0, gam = 1e10, om = 2, vertical = FALSE)
alpha.C = C2beta(t(X)%*%X, Y, ml$a, lambda = 0)$alpha
alpha.C

# PLS
m = 5
r = 2
L = rep(10, m)

Q = randortho(p, type = "orthonormal")
U = Q[,1:m]%*%DIAG(L)
D = t(U)%*%U
P = randortho(n, type = "orthonormal")
S = matrix(P[1:m,], ncol = n)
E = t(mvrnorm(n, rep(0, p), diag(p))*0.01)

X = U%*%S + E

Z = randortho(m, type = "orthonormal")[,1:r]
#Z = diag(m)[,1:r]
alpha = c(10, 5)
e = rnorm(n)*0.01
a = t(S)%*%DIAG(1/L)%*%Z
Y = t(S)%*%DIAG(L)%*%Z%*%alpha #+ e

# beta1 = rep(c(1, -1), n/2)
# beta2 = rep(1/2, n)
# Y = t(U[,1]%*%t(S[1,]))%*%U[,1]%*%t(S[1,])%*%beta1 +t(U[,2]%*%t(S[2,]))%*%U[,2]%*%t(S[2,])%*%beta2

# beta = rep(1, n)
# Y = t(U%*%(S))%*%U%*%(S)%*%beta

# eigen(t(X)%*%X%*%Y%*%t(Y)%*%t(X)%*%X)$value

#a = t(S)%*%DIAG(1/L)%*%Z
#Y = t(S)%*%DIAG(L)%*%Z%*%alpha + e
#t(a)%*%t(X)%*%X%*%t(X)%*%X%*%a
#t(a)%*%t(U%*%S)%*%U%*%S%*%t(U%*%S)%*%U%*%S%*%a
#t(a)%*%t(X)%*%X%*%Y

ml = continuum.ridge.fix(t(X), Y, 1, lambda = 0, gam = 1, om = 2, vertical = FALSE)
t(ml$a)%*%t(X)%*%X%*%t(X)%*%X%*%ml$a
alpha.C = C2beta(t(X)%*%X, Y, ml$a, lambda = 0)$alpha
alpha.C

# -----------------------------------------
m = 10
X = mvrnorm(n, rep(0, p), diag(p))
beta = rep(1/10, n)
e = rnorm(n)
Y = (X)%*%t(X)%*%beta #+ e
ml = continuum.ridge.fixm(X, Y, 1, lambda = 0, gam = 1, om = 2, m = m, vertical = FALSE)
#solve(t(ml$a)%*%(X)%*%t(X)%*%(X)%*%t(X)%*%ml$a)%*%t(ml$a)%*%(X)%*%t(X)%*%Y
alpha = C2beta((X)%*%t(X), Y, ml$a, lambda = 0)$alpha
alpha

ml$a%*%solve(t(ml$a)%*%ml$a)%*%t(ml$a)
ml$Z%*%solve(t(ml$Z)%*%solve(ml$E)%*%ml$Z)%*%t(ml$Z)%*%solve(ml$E)^(1/2)%*%t(ml$V)
# solve(t(ml$a)%*%ml$a)%*%t(ml$a) - solve(t(ml$Z)%*%solve(ml$E)%*%ml$Z)%*%t(ml$Z)%*%solve(ml$E)^(1/2)%*%t(ml$V)

t(ml$Z)%*%ml$E%*%ml$Z

U = Q[,1:m]%*%(ml$E[1:m, 1:m]^(1/2))
S = t(ml$V[,1:m])
X = U%*%S + E
a = t(S)%*%solve(ml$E[1:m, 1:m])^(1/2)%*%ml$Z
t(a)%*%t(X)%*%(X)%*%t(X)%*%(X)%*%a
#Y = t(X)%*%(X)%*%a%*%alpha + e
Y = t(U%*%S)%*%(U%*%S)%*%a%*%alpha #+ e

ml = continuum.ridge.fix(t(X), Y, 1, lambda = 0, gam = 1, om = 2, vertical = FALSE)
alpha.C = C2beta(t(X)%*%X, Y, ml$a, lambda = 0)$alpha
alpha.C

ml.pcr = continuum.ridge.fix(t(X), Y, 1, lambda = 0, gam = 1e10, om = 2, vertical = FALSE)
alpha.C = C2beta(t(X)%*%X, Y, ml.pcr$a, lambda = 0)$alpha
alpha.C

 
# testing 
P = randortho(n, type = "orthonormal")
S = matrix(P[1:m,], ncol = n)
E = t(mvrnorm(n, rep(0, p), diag(p))*0.01)
X.test = U%*%S + E
#Z = randortho(m, type = "orthonormal")[,1:r]
#alpha = c(3, 1)
e = rnorm(n)*0.01
#a = t(S)%*%DIAG(1/L)%*%Z
Y.test = t(S)%*%DIAG(L)%*%Z%*%alpha + e
