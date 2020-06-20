G = 2
result = as.matrix(result)
n = nrow(result)
result.list = list()
for (i in 1:n){
  result.list[[i]] = matrix(as.vector(result[i,-1]), ncol = G)
}

result.mean = matrix(apply(result, 2, mean)[-1], ncol = G)[c(L, (2*L):(2*L+6)), ]
result.sd = matrix(apply(result, 2, sd)[-1], ncol = G)[c(L, (2*L):(2*L+6)), ]
result.mean%*%c(n1, n2)/n

row.names(result.mean) = c("2step", "iter", "global.pls", "global.pcr", "global.ridge", "group.pls", "group.pcr", "group.ridge")

g = 1
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 2
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[1:(L-1),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L,g])

g = 1
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

g = 2
plot(gam.list, matrix(apply(result, 2, mean)[-1], ncol = G)[L + (1:(L-1)),g], xlab = "gam", ylab = "MSE")
abline(h = matrix(apply(result, 2, mean)[-1], ncol = G)[L*2, g])

n = nrow(rank)
rank.list = list()
for (i in 1:n){
  rank.list[[i]] = matrix(rank[i,-1], ncol = G+1)
}


###################################
#PLS

# groupwise
result = as.matrix(result)
n = nrow(result)
L = 50
a = seq(0, 1, length.out = L+1)

MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(n)
round(rbind(MSE[2*(L+1)+26,], MSE[-(1:(5*(L+1))),]), digit = 3)
round(rbind(se[2*(L+1)+26,], se[-(1:(5*(L+1))),]), digit = 3)

# joint
result = as.matrix(result)
n = nrow(result)
L = 50
a = seq(0, 1, length.out = L+1)

MSE.list = list()
for (i in 1:n){
  MSE.list[[i]] = apply(matrix(as.vector(result[i,-1]), ncol = G), 1, mean)
}

MSE = apply(do.call(rbind, MSE.list), 2, mean)
se = apply(do.call(rbind, MSE.list), 2, sd)/sqrt(n)

c(MSE[2*(L+1)+26], MSE[-(1:(5*(L+1)))])
c(se[2*(L+1)+26], se[-(1:(5*(L+1)))])


###################################
# jico
result = as.matrix(result)
n1 = 50
n2 = 50
n = n1 + n2

G = 2
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)

# group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
B = round(cbind(MSE, se)[,ix.MSEse], digits = 3)

A.jico = rbind(A[2*(L+1)+c(1,26,51),], A[-(1:(5*(L+1))),])
B.jico = rbind(B[2*(L+1)+c(1,26,51),], B[-(1:(5*(L+1))),])

# factor
result = as.matrix(result)
n1 = 50
n2 = 50
n = n1 + n2

G = 2
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)

# group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
B = round(cbind(MSE, se)[,ix.MSEse], digits = 3)

A.factor = A[c(1,4),] # only ridge 1,1 result
B.factor = B[c(1,4),] # only ridge 1,1 result

# combine jico and factor
A = rbind(A.jico[1:3,], A.jico[-(1:3),], A.factor)
B = rbind(B.jico[1:3,], B.jico[-(1:3),], B.factor)

C = apply(cbind(B, A), 2, function(x) sprintf(x, fmt = "%0.3f"))
parenth = c(rep(c("\ (", ")\ &\ "), G), c("\ (", ")\ \\\\"))
n.method = nrow(C)
D = cbind(C, matrix(rep(parenth, n.method), nrow = n.method, byrow = T))
ncol.D = ncol(D)
printt = D[, do.call(c, lapply(1:(ncol.D/2), function(i) c(i,(i+ncol.D/2))))]
for (i in 1: nrow(printt)){
  cat(c("& ", printt[i,]), sep = "", fill = T)
}



# OLS1
# jico
result = as.matrix(result)
n1 = 50
n2 = 50
n = n1 + n2

G = 2
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)

# group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
B = round(cbind(MSE, se)[,ix.MSEse], digits = 3)

A.jico = rbind(A[0*(L+1)+c(1,26,51),], A[-(1:(5*(L+1))),])
B.jico = rbind(B[0*(L+1)+c(1,26,51),], B[-(1:(5*(L+1))),])

# factor
result = as.matrix(result)
n1 = 50
n2 = 50
n = n1 + n2

G = 2
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)

# group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
B = round(cbind(MSE, se)[,ix.MSEse], digits = 3)

A.factor = A[c(1,4),] # only ridge 1,1 result
B.factor = B[c(1,4),] # only ridge 1,1 result

# combine jico and factor
A = rbind(A.jico[1:3,], A.jico[-(1:3),], A.factor)
B = rbind(B.jico[1:3,], B.jico[-(1:3),], B.factor)

C = apply(cbind(B, A), 2, function(x) sprintf(x, fmt = "%0.3f"))
parenth = c(rep(c("\ (", ")\ &\ "), G), c("\ (", ")\ \\\\"))
n.method = nrow(C)
D = cbind(C, matrix(rep(parenth, n.method), nrow = n.method, byrow = T))
ncol.D = ncol(D)
printt = D[, do.call(c, lapply(1:(ncol.D/2), function(i) c(i,(i+ncol.D/2))))]
for (i in 1: nrow(printt)){
  cat(c("& ", printt[i,]), sep = "", fill = T)
}

# OLS2
# jico
result = as.matrix(result)
n1 = 50
n2 = 50
n = n1 + n2

G = 2
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)

# group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
B = round(cbind(MSE, se)[,ix.MSEse], digits = 3)

A.jico = rbind(A[3*(L+1)+c(1,26,51),], A[-(1:(5*(L+1))),])
B.jico = rbind(B[3*(L+1)+c(1,26,51),], B[-(1:(5*(L+1))),])

# factor
result = as.matrix(result)
n1 = 50
n2 = 50
n = n1 + n2

G = 2
N = nrow(result)
result.list = list()
MSE.list = list()
for (i in 1:N){
  result.list[[i]] = matrix(result[i,-1], ncol = G)
  MSE.list[[i]] = matrix(result[i,-1], ncol = G)%*%c(n1, n2)/n
}

# overall
MSE = apply(do.call(cbind, MSE.list), 1, mean)
se = apply(do.call(cbind, MSE.list), 1, sd)/sqrt(N)
A = round(cbind(MSE, se), digit = 3)

# group
MSE = matrix(apply(result, 2, mean)[-1], ncol = G)
se = matrix(apply(result, 2, sd)[-1], ncol = G)/sqrt(N)
ix.MSEse = do.call(c, lapply(1:G, function(g) c(1, 1+G) + g-1))
B = round(cbind(MSE, se)[,ix.MSEse], digits = 3)

A.factor = A[c(1,4),] # only ridge 1,1 result
B.factor = B[c(1,4),] # only ridge 1,1 result

# combine jico and factor
A = rbind(A.jico[1:3,], A.jico[-(1:3),], A.factor)
B = rbind(B.jico[1:3,], B.jico[-(1:3),], B.factor)

C = apply(cbind(B, A), 2, function(x) sprintf(x, fmt = "%0.3f"))
parenth = c(rep(c("\ (", ")\ &\ "), G), c("\ (", ")\ \\\\"))
n.method = nrow(C)
D = cbind(C, matrix(rep(parenth, n.method), nrow = n.method, byrow = T))
ncol.D = ncol(D)
printt = D[, do.call(c, lapply(1:(ncol.D/2), function(i) c(i,(i+ncol.D/2))))]
for (i in 1: nrow(printt)){
  cat(c("& ", printt[i,]), sep = "", fill = T)
}

















