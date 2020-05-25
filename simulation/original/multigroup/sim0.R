library(MASS)
r = 3
r1 = 1
r2 = 1
n1 = 50
n2 = 50
n = n1 + n2
p = 100 # 40

X1 = mvrnorm(n1, rep(0, p), diag(p))
X2 = mvrnorm(n2, rep(0, p), diag(p))

s = .1

beta = rep(1/5, p)
beta1 = rep(c(1/100, 0, -1/100, 0), p/4)*s
beta2 = rep(c(0, -1/100, 0, 1/100), p/4)*s

e1 = rnorm(n1)*0.1
e2 = rnorm(n2)*0.1

Y1 = X1%*%beta + X1%*%beta1 + e1
Y2 = X2%*%beta + X2%*%beta2 + e2

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 1.1, rankJ = r, rankA = c(r1, r2), 
                          center.X = F, scale.X = F, center.Y = F, scale.Y = F, conv = 1e-6)
t(ml$Cind[[1]])%*%ml$Cind[[2]]
t(ml$C)%*%ml$Cind[[1]]
t(ml$C)%*%ml$Cind[[2]]

X1.test = mvrnorm(n1, rep(0, p), diag(p))
X2.test = mvrnorm(n2, rep(0, p), diag(p))

e1 = rnorm(n1)*0.1
e2 = rnorm(n2)*0.1

Y1.test = X1.test%*%beta + X1.test%*%beta1 + e1
Y2.test = X2.test%*%beta + X2.test%*%beta2 + e2

sum((as.numeric(ml$intercept[[1]])+ X1.test%*%ml$beta.C[[1]] + X1.test%*%ml$beta.Cind[[1]] - Y1.test)^2)
sum((as.numeric(ml$intercept[[2]])+X2.test%*%ml$beta.C[[2]] + X2.test%*%ml$beta.Cind[[2]] - Y2.test)^2)

sum(Y1.test^2)
sum(Y2.test^2)

