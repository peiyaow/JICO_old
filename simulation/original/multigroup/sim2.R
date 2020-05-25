r = 2
r1 = 1
r2 = 1
r.list = list(r1, r2)
n1 = 50
n2 = 50
n = n1 + n2
p = 100 # 40
G = 2

s = 0.1
beta = c(rep(1/25, 50), rep(0, p-50))
beta1 = c(rep(1/50, 25), rep(0, p-25))*s
beta2 = c(rep(-1/50, 25), rep(0, p-25))*s

S1 = mvrnorm(n1, rep(0, r), diag(r))
S2 = mvrnorm(n2, rep(0, r), diag(r))

U = mvrnorm(r, rep(0, p), diag(p))

T1 = mvrnorm(n1, rep(0, r1), diag(r1))
T2 = mvrnorm(n2, rep(0, r2), diag(r2))

U1 = mvrnorm(r1, rep(0, p), diag(p))
U2 = mvrnorm(r2, rep(0, p), diag(p))

E1 = mvrnorm(n1, rep(0, p), diag(p))*0.1
E2 = mvrnorm(n2, rep(0, p), diag(p))*0.1

X1 = S1%*%U + T1%*%U1 #+ E1
X2 = S2%*%U + T2%*%U2 #+ E2

e1 = rnorm(n1)*0.1
e2 = rnorm(n2)*0.1

Y1 = (S1%*%U + T1%*%U1)%*%beta + (T1%*%U1)%*%beta1 #+ e1
Y2 = (S2%*%U + T2%*%U2)%*%beta + (T2%*%U2)%*%beta2 #+ e2

X = rbind(X1, X2)
Y = rbind(Y1, Y2)

sum(X1^2)
sum(X2^2)

sum(Y1^2)
sum(Y2^2)

X.list = list(X1, X2)
Y.list = list(Y1, Y2)

ml = continuum.multigroup.iter(X.list, Y.list, lambda = 0, gam = 0.8, rankJ = r, rankA = c(r1, r2), 
                               center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = T)

R.list = lapply(1:G, function(g) X.list[[g]] - ml$J[[g]] - ml$I[[g]])

R1 = mean(R.list[[1]]^2)
R2 = mean(R.list[[2]]^2)

sum(R1^2)
sum(R2^2)

sum(((ml$J[[1]] + ml$I[[1]])%*%(beta+beta1) - Y1)^2)
sum(((ml$J[[2]] + ml$I[[2]])%*%(beta+beta2) - Y2)^2)

nr_x = 0.1
E1 = mvrnorm(n1, rep(0, p), diag(p))*nr_x
E2 = mvrnorm(n2, rep(0, p), diag(p))*nr_x
X1 = ml$J[[1]] + ml$I[[1]] + E1
X2 = ml$J[[2]] + ml$I[[2]] + E2

mean(E1^2)/mean(X1^2)*100
mean(E2^2)/mean(X2^2)*100

beta1 = ml$beta.Cind[[1]]
beta2 = ml$beta.Cind[[2]]
beta = ml$beta.C[[1]]

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









