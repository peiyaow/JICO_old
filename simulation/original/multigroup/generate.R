library(pracma)

r = 2
r1 = 1
r2 = 1
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 200 # 40
G = 2

s = 100
beta = c(rep(1/25, 50), rep(0, p-50))*s
beta1 = c(rep(0, p-25), rep(1/50, 25))*s
beta2 = c(rep(0, p-50), rep(-1/50, 25), rep(0, 25))*s

S1 = mvrnorm(n1, rep(0, r), diag(r))
S2 = mvrnorm(n2, rep(0, r), diag(r))

U = mvrnorm(r, rep(0, p), diag(p))#%*%V%*%t(V)

T1 = mvrnorm(n1, rep(0, r1), diag(r1))
T2 = mvrnorm(n2, rep(0, r2), diag(r2))

U1 = mvrnorm(r1, rep(0, p), diag(p))#%*%V1%*%t(V1)
U2 = mvrnorm(r2, rep(0, p), diag(p))#%*%V2%*%t(V2)

E1 = mvrnorm(n1, rep(0, p), diag(p))*.1
E2 = mvrnorm(n2, rep(0, p), diag(p))*.1

X1 = S1%*%U + T1%*%U1 #+ E1
X2 = S2%*%U + T2%*%U2 #+ E2

# X1 = T1%*%U1 + E1
# X2 = T2%*%U2 + E2

# e1 = rnorm(n1)*0.1
# e2 = rnorm(n2)*0.1

# Y1 = (T1%*%U1)%*%(beta + beta1) + e1
# Y2 = (T2%*%U2)%*%(beta + beta2) + e2

Y1 = (S1%*%U + T1%*%U1)%*%(beta1) #+ e1
Y2 = (S2%*%U + T2%*%U2)%*%(beta2) #+ e2

# Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) #+ e1
# Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) #+ e2

#Y1 = (S1%*%U + T1%*%U1)%*%beta + (T1%*%U1)%*%beta1 + e1
#Y2 = (S2%*%U + T2%*%U2)%*%beta + (T2%*%U2)%*%beta2 + e2

#Y1 = (S1%*%U)%*%beta + (T1%*%U1)%*%beta1 #+ e1
#Y2 = (S2%*%U)%*%beta + (T2%*%U2)%*%beta2 #+ e2

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200, gam = 1, rankJ = r, rankA = c(r1, r2),
                               center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = T)
                               
R.list = lapply(1:G, function(g) X.list[[g]] - ml$J[[g]] - ml$I[[g]])

R1 = mean(R.list[[1]]^2)
R2 = mean(R.list[[2]]^2)

sum(R1^2)
sum(R2^2)

beta1 = ml$beta.Cind[[1]]
beta2 = ml$beta.Cind[[2]]
beta = ml$beta.C[[1]]
sum(beta1^2)
sum(beta2^2)
sum(beta^2)

sum(((ml$J[[1]] + ml$I[[1]])%*%(beta+beta1) - Y1)^2)
sum(((ml$J[[2]] + ml$I[[2]])%*%(beta+beta2) - Y2)^2)

J = do.call(rbind, ml$J)
I1 = ml$I[[1]]
I2 = ml$I[[2]]

t(ml$Cind[[1]])%*%t(I1)%*%I1%*%ml$Cind[[1]]
t(ml$Cind[[2]])%*%t(I2)%*%I2%*%ml$Cind[[2]]
t(ml$C)%*%t(J)%*%J%*%ml$C

svd.J = svd(J, nu = r, nv = r)
V = svd.J$v
U = svd.J$u
D = DIAG(svd(J, nu = r, nv = r)$d[1:r])


svd.I1 = svd(I1, nu = r1, nv = r1)
V1 = svd.I1$v
U1 = svd.I1$u
D1 = DIAG(svd.I1$d[1:r1])

svd.I2 = svd(I2, nu = r2, nv = r2)
V2 = svd.I2$v
U2 = svd.I2$u
D2 = DIAG(svd.I2$d[1:r2])

save(n1, n2, p, n, r1, r2, r, G, beta1, beta2, beta, V1, V2, V, U1, U2, U, D1, D2, D, file = "pls.RData")

