library(pracma)
library(MASS)

r = 2
r1 = 2
r2 = 2
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 200 # 40
G = 2

U = matrix(randortho(n, type = "orthonormal")[,1:r], ncol = r)
J = U%*%D%*%t(V)
J1 = J[1:n1,]
J2 = J[n1+(1:n2),]

U1 = matrix(randortho(n1, type = "orthonormal")[,1:r1], ncol = r1)
I1 = U1%*%D1%*%t(V1)

U2 = matrix(randortho(n2, type = "orthonormal")[,1:r2], ncol = r2)
I2 = U2%*%D2%*%t(V2)


E1 = mvrnorm(n1, rep(0, p), diag(p))*.2
E2 = mvrnorm(n2, rep(0, p), diag(p))*.2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

sum(E1^2)/sum(X1^2)
sum(E2^2)/sum(X2^2)

e1 = rnorm(n1)*1
e2 = rnorm(n2)*1

Y1 = J1%*%beta + I1%*%beta1 + e1
Y2 = J2%*%beta + I2%*%beta2 + e2

sum(e1^2)/sum(Y1^2)
sum(e2^2)/sum(Y2^2)

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

L = 10
gam.list = exp(seq(log(0.5), log(1), length.out = L-1))
parameter.set = list()
for (gam in gam.list){
  ml.step1 = cv.continuum.step1(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m = 5, 
                                center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "scree")
  rankJ = ml.step1$rankJ
  ml.step2 = lapply(1:G, function(g) cv.continuum(X.list[[g]], Y.list[[g]], lambda = 0, gam = gam, nfolds = 10, m = 5, 
                                                  center.X = F, scale.X = T, center.Y = F, scale.Y = T, criteria = "1se"))
  rankA = sapply(ml.step2, function(ml) max(ml$rankA - rankJ, 0))
  parameter = list(gam = gam, rankJ = rankJ, rankA = rankA)
  parameter.set = list.append(parameter.set, parameter)
}

ml.2step.best = cv.continnum.2step(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, 
                                   center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                   criteria = "min")

ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, 
                                 center.X = F, scale.X = T, center.Y = F, scale.Y = T, 
                                 maxiter = 200, criteria = "min")



