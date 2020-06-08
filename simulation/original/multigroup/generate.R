library(pracma)
library(MASS)

path = "~/Documents/GitHub/continuum/"
setwd(path)
source("./function/jive_continuum.R")
source("./function/cv_multigroup.R")

r = 2
r1 = 1
r2 = 1
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 200 # 40
G = 2

# Q = randortho(p)
# V = matrix(Q[,1:r], ncol = r)
# V1 = matrix(Q[, r + (1:r1) ], ncol = r1)
# V2 = matrix(Q[, r+r1 + (1:r2) ], ncol = r2)

s = 10
# beta = V%*%t(V)%*%c(rep(1/25, 50), rep(0, p-50))*s
# beta1 = V1%*%t(V1)%*%c(rep(0, p-25), rep(1/50, 25))*s
# beta2 = V2%*%t(V2)%*%c(rep(0, p-25), rep(-1/50, 25))*s
# 
# sum(beta^2)
# sum(beta1^2)
# sum(beta2^2)

beta = c(rep(1/25, 50), rep(0, p-50))*s
beta1 = c(rep(0, p-25), rep(1/50, 25))*s
beta2 = c(rep(0, p-50), rep(-1/50, 25), rep(0, 25))*s

t = 1

S1 = mvrnorm(n1, rep(0, r), diag(r)*t)
S2 = mvrnorm(n2, rep(0, r), diag(r)*t)
U = mvrnorm(r, rep(0, p), diag(p))#%*%V%*%t(V)

t1 = 0.5
t2 = 0.5
T1 = mvrnorm(n1, rep(0, r1), diag(r1)*t1)
T2 = mvrnorm(n2, rep(0, r2), diag(r2)*t2)
U1 = mvrnorm(r1, rep(0, p), diag(p))#%*%V1%*%t(V1)
U2 = mvrnorm(r2, rep(0, p), diag(p))#%*%V2%*%t(V2)

E1 = mvrnorm(n1, rep(0, p), diag(p))*.2
E2 = mvrnorm(n2, rep(0, p), diag(p))*.2

X1 = S1%*%U + T1%*%U1 + E1
X2 = S2%*%U + T2%*%U2 + E2

sum(E1^2)/sum(X1^2)
sum(E2^2)/sum(X2^2)

e1 = rnorm(n1)*1
e2 = rnorm(n2)*1

Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2

sum(e1^2)/sum(Y1^2)
sum(e2^2)/sum(Y2^2)

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 1000, gam = 0.8, rankJ = r, rankA = c(r1, r2),
                               center.X = F, scale.X = T, center.Y = F, scale.Y = F, orthIndiv = F)

ml2 = continuum.multigroup.iter2(X.list, Y.list, lambda = 0, maxiter = 1000, gam = 0.8, rankJ = r, rankA = c(r1, r2),
                                center.X = F, scale.X = T, center.Y = F, scale.Y = F, orthIndiv = F)




R.list = lapply(1:G, function(g) X.list[[g]] - ml$J[[g]]*ml$scaleValues.X[[g]] - ml$I[[g]]*ml$scaleValues.X[[g]])

R1 = R.list[[1]]
R2 = R.list[[2]]

sum(R1^2)
sum(R2^2)

r.list = lapply(1:G, function(g) Y.list[[g]] - ml$J[[g]]%*%ml$beta*ml$scaleValues.Y[[g]] - ml$I[[g]]%*%ml$beta_i[[g]]*ml$scaleValues.Y[[g]])

r1 = mean(r.list[[1]]^2)
r2 = mean(r.list[[2]]^2)

sum(r1^2)
sum(r2^2)

beta1 = ml$beta.Cind[[1]]
beta2 = ml$beta.Cind[[2]]
beta = ml$beta.C[[1]]

sum(beta^2)
sum(beta1^2)
sum(beta2^2)

C = ml$C
C1 = ml$Cind[[1]]
C2 = ml$Cind[[2]]

t(C)%*%C1
t(C)%*%C2
t(C1)%*%C2

# ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = 1, nfolds = 10, m = 10,
#                               center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "scree", plot = T)
# 
# rankJ = ml.step1$rankJ
# ml.step2 = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = 1, nfolds = 10, m = 10,
#                                                 center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "1se", plot = T))
# rankA = sapply(ml.step2, function(ml) ml$rankA)
# rankJ
# rankA
# 
# ml.step2[[1]]$rMSE

# ml = continuum.multigroup.iter2(X.list, Y.list, lambda = 0, maxiter = 200, gam = 1, rankJ = r, rankA = c(r1, r2),
#                                 center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = T)
# 
# R.list = lapply(1:G, function(g) X.list[[g]] - ml$J[[g]]*ml$scaleValues.X[[g]] - ml$I[[g]]*ml$scaleValues.X[[g]])
# 
# R1 = mean(R.list[[1]]^2)
# R2 = mean(R.list[[2]]^2)
# 
# sum(R1^2)
# sum(R2^2)
# 
# beta1 = ml$beta.Cind[[1]]
# beta2 = ml$beta.Cind[[2]]
# beta = ml$beta.C[[1]]
# 
# sum(beta^2)
# sum(beta1^2)
# sum(beta2^2)
# 
# C = ml$C
# C1 = ml$Cind[[1]]
# C2 = ml$Cind[[2]]
# 
# t(C)%*%C1
# t(C)%*%C2
# t(C1)%*%C2

# save data
r = 2
r1 = 1
r2 = 1

J1 = ml$J[[1]]
J2 = ml$J[[2]]
J = do.call(rbind, ml$J)
I1 = ml$I[[1]]
I2 = ml$I[[2]]

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

save(n1, n2, n, p,
     r1, r2, r, G,
     beta1, beta2, beta,
     V1, V2, V,
     U1, U2, U,
     D1, D2, D,
     C1, C2, C,
     file = "pls.RData")

