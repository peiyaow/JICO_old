library(pracma)

r = 2
r1 = 2
r2 = 2
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 200 # 40
G = 2

s = 10

beta = c(rep(1/25, 50), rep(0, p-50))*10
beta1 = c(rep(0, p-25), rep(1/50, 25))*s*10
beta2 = c(rep(0, p-25), rep(-1/50, 25))*s*10

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

e1 = rnorm(n1)*0.1
e2 = rnorm(n2)*0.1

# Y1 = (T1%*%U1)%*%(beta + beta1) + e1
# Y2 = (T2%*%U2)%*%(beta + beta2) + e2

#Y1 = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
#Y2 = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2

#Y1 = (S1%*%U + T1%*%U1)%*%beta + (T1%*%U1)%*%beta1 + e1
#Y2 = (S2%*%U + T2%*%U2)%*%beta + (T2%*%U2)%*%beta2 + e2

Y1 = (S1%*%U)%*%beta + (T1%*%U1)%*%beta1 #+ e1
Y2 = (S2%*%U)%*%beta + (T2%*%U2)%*%beta2 #+ e2

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200, gam = 1, rankJ = r, rankA = c(r1, r2),
                               center.X = F, scale.X = F, center.Y = F, scale.Y = F)
                               
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
D = DIAG(svd(J, nu = 2, nv = 2)$d[1:r])


svd.I1 = svd(I1, nu = r1, nv = r1)
V1 = svd.I1$v
U1 = svd.I1$u
D1 = DIAG(svd.I1$d[1:r1])

svd.I2 = svd(I2, nu = r2, nv = r2)
V2 = svd.I2$v
U2 = svd.I2$u
D2 = DIAG(svd.I2$d[1:r2])

# save(beta1, beta2, beta, V1, V2, V, U1, U2, U, D1, D2, D, file = "pls.RData")


Unew = matrix(randortho(n, type = "orthonormal")[,1:r], ncol = r)
Jnew = Unew%*%D%*%t(V)

U1new = matrix(randortho(n1, type = "orthonormal")[,1:r1], ncol = r1)
I1new = U1new%*%D1%*%t(V1)

U2new = matrix(randortho(n2, type = "orthonormal")[,1:r2], ncol = r2)
I2new = U2new%*%D2%*%t(V2)

t(ml$C)%*%t(Jnew)%*%Jnew%*%ml$C
t(ml$Cind[[1]])%*%t(I1new)%*%I1new%*%ml$Cind[[1]]
t(ml$Cind[[2]])%*%t(I2new)%*%I2new%*%ml$Cind[[2]]








# testing data
S1 = mvrnorm(n1, rep(0, r), diag(r))
S2 = mvrnorm(n2, rep(0, r), diag(r))

#U = mvrnorm(r, rep(0, p), diag(p))

T1 = mvrnorm(n1, rep(0, r1), diag(r1))
T2 = mvrnorm(n2, rep(0, r2), diag(r2))

#U1 = mvrnorm(r1, rep(0, p), diag(p))
#U2 = mvrnorm(r2, rep(0, p), diag(p))

E1 = mvrnorm(n1, rep(0, p), diag(p))*.1
E2 = mvrnorm(n2, rep(0, p), diag(p))*.1

# X1.test = T1%*%U1 + E1
# X2.test = T2%*%U2 + E2

X1.test = S1%*%U + T1%*%U1 + E1
X2.test = S2%*%U + T2%*%U2 + E2

e1 = rnorm(n1)*0.1
e2 = rnorm(n2)*0.1

# Y1.test = (T1%*%U1)%*%(beta + beta1) + e1
# Y2.test = (T2%*%U2)%*%(beta + beta2) + e2

# Y1.test = (S1%*%U + T1%*%U1)%*%(beta + beta1) + e1
# Y2.test = (S2%*%U + T2%*%U2)%*%(beta + beta2) + e2

#Y1.test = (S1%*%U + T1%*%U1)%*%beta + (T1%*%U1)%*%beta1 + e1
#Y2.test = (S2%*%U + T2%*%U2)%*%beta + (T2%*%U2)%*%beta2 + e2

Y1.test = (S1%*%U)%*%beta + (T1%*%U1)%*%beta1 + e1
Y2.test = (S2%*%U)%*%beta + (T2%*%U2)%*%beta2 + e2

X.test.list = list(X1.test, X2.test)
Y.test.list = list(Y1.test, Y2.test)


ml.pls = plsr(Y~X, validation = "CV", center = F)
ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
ncomp.pls

ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F))
ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))
ncomp.pls.list

L = 6
gam.list = exp(seq(log(0.5), log(1), length.out = L-1))
ml.2step.list = lapply(gam.list, function(gam) cv.continnum.2step.separate(X.list, Y.list, lambda = 0, gam = gam, nfolds = 10, m1 = 10, m2 = 5,
                                                                           center.X = F, scale.X = F, center.Y = F, scale.Y = F))

parameter.set = lapply(ml.2step.list, function(ml) ml$parameter) # parameter set

# tune best 2step model
ml.2step.best = cv.continnum.2step(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, criteria = "min", 
                                   center.X = F, scale.X = F, center.Y = F, scale.Y = F)

#tune best iterate model
ml.iter.best = cv.continnum.iter(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 200, criteria = "min",
                                center.X = F, scale.X = F, center.Y = F, scale.Y = F)

# ml.iter.list = list()
# for (parameter in parameter.set){
#   print(parameter)
#   ml.iter.list = list.append(ml.iter.list, 
#                              continuum.multigroup.iter(X.list, Y.list, lambda = 0, maxiter = 200,
#                                                        gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA,
#                                                        center.X = F, scale.X = T, center.Y = F, scale.Y = T))
# }

# ml.list = list()
# for (gam in gam.list){
#   print(gam)
#   ml.list = list.append(ml.list, continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = gam, rankJ = r, rankA = c(r1, r2),
#                                                            center.X = F, scale.X = T, center.Y = F, scale.Y = T))
# }

MSE = list()
MSE.2step = list()
for (ml in ml.2step.list){
  MSE.2step = list.append(MSE.2step, sapply(1:G, function(g) 
    mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
}
MSE.2step = list.append(MSE.2step, MSE.2step[[ml.2step.best$ix]])
MSE[1:L] = MSE.2step


# for (ml in ml.list){
#   MSE = list.append(MSE, sapply(1:G, function(g) 
#     mean((as.numeric(ml$intercept[[g]])+ X.test.list[[g]]%*%ml$beta.C[[g]] + X.test.list[[g]]%*%ml$beta.Cind[[g]] - Y.test.list[[g]])^2)))
# }

ml = ml.pls
# MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pls)[,,1] - Y.test.list[[g]])^2)))

MSE = list.append(MSE, sapply(1:G, function(g) mean((predict(ml, newdata = X.test.list[[g]], ncomp = ncomp.pls)[,,1] - Y.test.list[[g]])^2)))



ml = ml.pls.list
MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = ncomp.pls.list[[g]])[,,1]
)^2)))

MSE
















sum(((ml$J[[1]] + ml$I[[1]])%*%(beta+beta1) - Y1)^2)
sum(((ml$J[[2]] + ml$I[[2]])%*%(beta+beta2) - Y2)^2)

nr_x = 0.1
E1 = mvrnorm(n1, rep(0, p), diag(p))*nr_x
E2 = mvrnorm(n2, rep(0, p), diag(p))*nr_x
X1 = ml$J[[1]] + ml$I[[1]] + E1
X2 = ml$J[[2]] + ml$I[[2]] + E2

mean(E1^2)/mean(X1^2)*100
mean(E2^2)/mean(X2^2)*100



nr_y = 0.05
e1 = rnorm(n1)*nr_y
e2 = rnorm(n2)*nr_y

Y1 = (ml$J[[1]] + ml$I[[1]])%*%(beta+beta1) + e1
Y2 = (ml$J[[2]] + ml$I[[2]])%*%(beta+beta2) + e2

mean(e1^2)/mean(Y1^2)*100
mean(e2^2)/mean(Y2^2)*100

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 0.8, rankJ = r, rankA = c(r1, r2), 
                               center.X = F, scale.X = F, center.Y = F, scale.Y = F)
R.list = lapply(1:G, function(g) X.list[[g]] - ml$J[[g]] - ml$I[[g]])
mean(R.list[[1]]^2)

# ml$J
mean(R.list[[2]]^2)

mean(E1^2)
mean(E2^2)









