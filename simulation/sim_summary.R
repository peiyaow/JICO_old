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

